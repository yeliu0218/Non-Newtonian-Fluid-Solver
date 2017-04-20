/*
 *  Copyright 1995-2010 by IRSN
 *
 *  This software is an application framework, with a set of integrated  
 *  reusable components, whose purpose is to simplify the task of developing 
 *  softwares of numerical mathematics and scientific computing.
 * 
 *  This software is governed by the CeCILL-C license under French law and 
 *  abiding by the rules of distribution of free software. You can use, modify 
 *  and/or redistribute the software under the terms of the CeCILL-C license  
 *  as circulated by CEA, CNRS and INRIA at the following URL 
 *  "http://www.cecill.info". 
 *
 *  As a counterpart to the access to the source code and rights to copy,  
 *  modify and redistribute granted by the license, users are provided only 
 *  with a limited warranty and the software's author, the holder of the  
 *  economic rights, and the successive licensors have only limited liability. 
 *
 *  In this respect, the user's attention is drawn to the risks associated  
 *  with loading, using, modifying and/or developing or reproducing the  
 *  software by the user in light of its specific status of free software,
 *  that may mean that it is complicated to manipulate, and that also  
 *  therefore means that it is reserved for developers and experienced 
 *  professionals having in-depth computer knowledge. Users are therefore 
 *  encouraged to load and test the software's suitability as regards their 
 *  requirements in conditions enabling the security of their systems and/or 
 *  data to be ensured and, more generally, to use and operate it in the same 
 *  conditions as regards security. 
 *
 *  The fact that you are presently reading this means that you have had 
 *  knowledge of the CeCILL-C license and that you accept its terms.
 */

#include <PEL_ObjectReader.hh>

#include <PEL_Communicator.hh>
#include <PEL_Data.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ModuleIterator.hh>
#include <PEL_assertions.hh>

#include <fstream>
#include <sstream>

using std::endl ;
using std::ifstream ;
using std::ostringstream ;
using std::string ;
using std::stack ;

struct PEL_ObjectReader_ERROR
{
   static void n0( string const& fname ) ;
   static void n1( string const& fname ) ;
   static void n2( size_t cycle_number, string const& name ) ;
   static void n3( string const& name ) ;
   static void n4( string const& class_name, string const& nn ) ;
   static void n5( void ) ;
   static void n6( string const& class_name ) ;
   static void n7( size_t stored_nb_rank, size_t nb_ranks ) ;
} ;

//---------------------------------------------------------------------------
PEL_ObjectReader*
PEL_ObjectReader:: create( PEL_Object* a_owner,
			   PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectReader:: create" ) ;

   PEL_ObjectReader* result = new PEL_ObjectReader( a_owner, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->positioned_in_a_valid_cycle() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PEL_ObjectReader:: PEL_ObjectReader( PEL_Object* a_owner,
                                     PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , ROOT_MOD( 0 )
   , NB_CYCLES( 0 )
   , LAST_CYCLE( 0 )
   , iOBJECT( 0 )
{
   string fname = exp->string_data( "file_name" ) ;
   PEL_Communicator const* com = PEL_Exec::communicator() ;
   if( com->nb_ranks()>1 )
   {
      ostringstream rank ;
      rank << "." << com->rank() ;
      fname  += rank.str() ;
   }

   ifstream file( fname.c_str(), std::ios::in  ) ;
   if( !file ) PEL_ObjectReader_ERROR::n0( fname ) ;
   file.close() ;

   ROOT_MOD = PEL_Module::create( this, "ROOT", fname ) ;

   // Check fname structure :
   PEL_ModuleIterator* it = ROOT_MOD->create_module_iterator( 0 ) ;
   if( !it->is_valid() ) 
      PEL_ObjectReader_ERROR::n1( fname ) ;
   if( !( it->item()->name()=="communicator") ) 
      PEL_ObjectReader_ERROR::n1( fname ) ;
   it->go_next() ;
   if( !it->is_valid() ) 
      PEL_ObjectReader_ERROR::n1( fname ) ;
   if( !( it->item()->name()=="header") ) 
      PEL_ObjectReader_ERROR::n1( fname ) ;
   for( it->go_next() ; it->is_valid() ; it->go_next() )
   {
      PEL_Module const* mm = it->item() ;
      if( !(mm->name().substr(0,6)=="cycle#") )
         PEL_ObjectReader_ERROR::n1( fname ) ;
      if( !mm->has_entry( "cycle_number" ) )
         PEL_ObjectReader_ERROR::n3( mm->name() ) ;
      int i = mm->data_of_entry( "cycle_number" )->to_int() ;
      if( i<0 ) PEL_ObjectReader_ERROR::n3( mm->name() ) ;
      LAST_CYCLE = (size_t) i ;
      ++NB_CYCLES ;
   }
   it->destroy() ;

   // Check communicator :
   PEL_ModuleExplorer const* com_exp =
      PEL_ModuleExplorer::create( 0, ROOT_MOD->module( "communicator" ) ) ;
   size_t const nb_ranks = com_exp->int_data( "nb_ranks" ) ;
   size_t const rank = com_exp->int_data( "rank" ) ;
   if( nb_ranks!=com->nb_ranks() )
   {
      PEL_ObjectReader_ERROR::n7( nb_ranks, com->nb_ranks() ) ;
   }
   if( rank!=com->rank() )
   {
      PEL_ObjectReader_ERROR::n1( fname ) ;
   }
   com_exp->destroy() ;
}

//---------------------------------------------------------------------------
PEL_ObjectReader:: ~PEL_ObjectReader( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
PEL_Module*
PEL_ObjectReader:: header_module( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectReader:: header_module" ) ;

   // the existence of such a module has been tested in the constructor
   PEL_Module* result = ROOT_MOD->module( "header" ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
PEL_ObjectReader:: nb_cycles( void ) const
//---------------------------------------------------------------------------
{
   size_t result = NB_CYCLES ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
PEL_ObjectReader:: seek_cycle( size_t cycle_number )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectReader:: seek_cycle" ) ;

   if( cycle_number==0 ) cycle_number=LAST_CYCLE ;

   ostringstream name ;
   name << "cycle#" << cycle_number ;
   if( !ROOT_MOD->has_module( name.str() ) ) 
      PEL_ObjectReader_ERROR::n2( cycle_number, name.str() ) ;

   PEL_Module* mm = ROOT_MOD->module( name.str() ) ;
   MODS.push( mm ) ;
   MOD_ITS.push( mm->create_module_iterator( this ) ) ; 

   if( mm->data_of_entry( "cycle_number" )->to_int()!=(int)cycle_number )
      PEL_ObjectReader_ERROR::n3( name.str() ) ;

   iOBJECT = 0 ;

   PEL_CHECK_POST( current_object_number() == 0 ) ;
   PEL_CHECK_POST( positioned_in_a_valid_cycle() ) ;
}

//---------------------------------------------------------------------------
bool
PEL_ObjectReader:: positioned_in_a_valid_cycle( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectReader:: positioned_in_a_valid_cycle" ) ;

   bool result = !MOD_ITS.empty() ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
PEL_ObjectReader:: close_cycle( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectReader:: close_cycle" ) ;
   PEL_CHECK_PRE( positioned_in_a_valid_cycle() ) ;

   iOBJECT = 0 ;

   PEL_ModuleIterator* it = MOD_ITS.top() ;
   MOD_ITS.pop() ;
   destroy_possession( it ) ;

   MODS.pop() ;

   if( !MODS.empty() ) PEL_ObjectReader_ERROR::n5() ;

   PEL_CHECK_POST( !positioned_in_a_valid_cycle() ) ;
}

//---------------------------------------------------------------------------
void
PEL_ObjectReader:: start_object_retrieval( std::string const& class_name )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectReader:: start_object_retrieval" ) ;
   PEL_CHECK_PRE( positioned_in_a_valid_cycle() ) ;

   if( !MOD_ITS.top()->is_valid() ) PEL_ObjectReader_ERROR::n6( class_name ) ;

   PEL_Module* mm = MOD_ITS.top()->item() ;
   MODS.push( mm ) ;
   MOD_ITS.push( mm->create_module_iterator( this ) ) ; 

   PEL_ModuleExplorer const* exp = PEL_ModuleExplorer::create( 0, mm ) ;
   string const& nn = exp->string_data( "class" ) ;
   if( nn != class_name ) PEL_ObjectReader_ERROR::n4( class_name, nn ) ;

   iOBJECT = exp->int_data( "object_number" ) ;
   if( ! ( iOBJECT > 0 ) ) 
      PEL_Error::object()->raise_bad_data_value( exp, 
                                                 "object_number", 
                                                 "greater or equal to 1" ) ;
   exp->destroy() ;

   PEL_CHECK_POST( current_object_number() != 0 ) ;
}

//---------------------------------------------------------------------------
size_t
PEL_ObjectReader:: current_object_number( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectReader:: current_object_number" ) ;

   return( (size_t)iOBJECT ) ;
}

//---------------------------------------------------------------------------
bool
PEL_ObjectReader:: has_entry( std::string const& keyword ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectReader:: has_entry" ) ;
   PEL_CHECK_PRE( positioned_in_a_valid_cycle() ) ;
   PEL_CHECK_PRE( current_object_number() != 0 ) ;
   
   return( MODS.top()->has_entry( keyword ) ) ;
}

//---------------------------------------------------------------------------
PEL_Data const*
PEL_ObjectReader:: data_of_entry( std::string const& keyword ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectReader:: data_of_entry" ) ;
   PEL_CHECK_PRE( positioned_in_a_valid_cycle() ) ;
   PEL_CHECK_PRE( current_object_number() != 0 ) ;

   PEL_Data const* result = MODS.top()->data_of_entry( keyword ) ;

   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
PEL_ObjectReader:: end_object_retrieval( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ObjectReader:: end_object_retrieval" ) ;

   PEL_ModuleIterator* it = MOD_ITS.top() ;
   MOD_ITS.pop() ;
   destroy_possession( it ) ;

   MODS.pop() ;

   MOD_ITS.top()->go_next() ;
}

//---------------------------------------------------------------------------
bool
PEL_ObjectReader:: invariant( void ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   PEL_ASSERT( MODS.size() == MOD_ITS.size() ) ;
   return( true ) ;
}


//internal--------------------------------------------------------------
void 
PEL_ObjectReader_ERROR:: n0( string const& fname )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "PEL_ObjectReader :" << endl
        << "   unable to open file \"" << fname << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PEL_ObjectReader_ERROR:: n1( string const& fname )
//internal--------------------------------------------------------------
{
   //????? faire un message plus explicite qui expose la structure
   //????? attendue du fichier telle que décrite dans la doc de 
   //????? PEL_ObjectWriter
   ostringstream mesg ;
   mesg << "PEL_ObjectReader :" << endl
        << "   file \"" << fname << "\" has an invalid structure" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PEL_ObjectReader_ERROR:: n2( size_t cycle_number, string const& name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "PEL_ObjectReader :" << endl
        << "   object retrieval from cycle " << cycle_number 
        << " is impossible " << endl 
        << "   since the underlying file has no module called \"" 
        << name << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PEL_ObjectReader_ERROR:: n3( string const& name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "PEL_ObjectReader :" << endl
        << "   module \"" << name << "\" has a missing or invalid" << endl
        << "   entry of keyword \"cycle_number\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PEL_ObjectReader_ERROR:: n4( string const& class_name, string const& nn )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "PEL_ObjectReader :" << endl
        << "   attempt to retrieve an object of" << endl 
        << "   class \"" << class_name << "\" from a module whose " << endl
        << "   entry of keyword \"class_name\" is \"" << nn << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PEL_ObjectReader_ERROR:: n5( void )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "PEL_ObjectReader :" << endl
        << "   impossible to close a cycle when an " << endl 
        << "   object retrieval is in progress" << endl
        << "   (\"start_object_retrieval\" and \"end_object_retrieval\"" 
        << endl
        << "    should be called the same number of times between two calls to"
	<< endl
	<< "    \"seek_cycle\" and \"close_cycle\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PEL_ObjectReader_ERROR:: n6( string const& class_name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "PEL_ObjectReader :" << endl
        << "   attempt to retrieve an object of" << endl 
        << "   class \"" << class_name << "\" from a module who " << endl
        << "   does not contain any more object storage"  ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PEL_ObjectReader_ERROR:: n7( size_t stored_nb_rank, size_t nb_ranks )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "PEL_ObjectReader :" << endl
        << "   object retrieval is impossible on "
        << nb_ranks << " processes" << endl 
        << "   because the storage has been done on "
        << stored_nb_rank << endl ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}
