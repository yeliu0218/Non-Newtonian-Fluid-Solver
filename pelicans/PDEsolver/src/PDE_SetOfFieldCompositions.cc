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

#include <PDE_SetOfFieldCompositions.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_FieldComposition.hh>
#include <PDE_FieldCompositionExp.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <PEL.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>
#include <PEL_assertions.hh>

#include <iostream>

//----------------------------------------------------------------------
PDE_SetOfFieldCompositions*
PDE_SetOfFieldCompositions:: create( PEL_Object* a_owner )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfFieldCompositions:: create" ) ;

   PDE_SetOfFieldCompositions* result =
                             new PDE_SetOfFieldCompositions( a_owner ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_compositions() == 0 ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_SetOfFieldCompositions:: PDE_SetOfFieldCompositions( PEL_Object* a_owner )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , COMPO_TABLE( PEL_Vector::create( this, 0 ) )
   , INDEX( PEL::bad_index() )
{
}

//----------------------------------------------------------------------
void
PDE_SetOfFieldCompositions:: build_compositions(
                                     PEL_ModuleExplorer* exp,
                                     size_t nb_sp_dimensions,
                                     PDE_SetOfDiscreteFields const* fields )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfFieldCompositions:: build_compositions" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   PEL_CHECK_PRE( nb_sp_dimensions == 0 || nb_sp_dimensions == 1 ||
                  nb_sp_dimensions == 2 || nb_sp_dimensions == 3 ) ;
   PEL_CHECK_PRE( fields!= 0 ) ;
      
   // Composition from expression :
   
   exp->start_entry_iterator() ;
   for( ; exp->is_valid_entry() ; exp->go_next_entry() )
   {
      PDE_FieldComposition* compo =
         PDE_FieldCompositionExp::create( 0,
                                          exp->keyword(), exp->data(0),
                                          fields ) ;
      if( has( compo->name() ) )
      {
         std::string msg = "*** Error while field composition building\n" ;
         msg += "    \""+compo->name()+"\" is defined twice." ;
         PEL_Error::object()->raise_plain( msg ) ;
      }
      append( compo ) ;
   }

   // Plug in compositions :
   exp->start_module_iterator() ;
   for( ; exp->is_valid_module() ;  exp->go_next_module() )
   {
      PEL_ModuleExplorer* e = exp->create_subexplorer( 0 ) ;

      PDE_FieldComposition* compo =
         PDE_FieldComposition::make( 0, e,  nb_sp_dimensions, fields ) ;
      if( has( compo->name() ) )
      {
         std::string msg = "*** Error while field composition building\n" ;
         msg += "    \""+compo->name()+"\" is defined twice." ;
         PEL_Error::object()->raise_plain( msg ) ;
      }
      append( compo ) ;
      e->destroy() ; e = 0 ;
   }

   // Do the links between compositions :
   for( start() ; is_valid() ; go_next() )
   {
      item()->do_the_links( this ) ;
   }

   // Finalize compositions :
   for( start() ; is_valid() ; go_next() )
   {
      item()->complete_internal_dependencies() ;
   }
   
   PEL_CHECK_POST( FORALL( ( start() ; is_valid() ; go_next() ),
                           item()->owner() == this ) ) ;
   PEL_CHECK_POST(
      FORALL(
         ( start() ; is_valid() ; go_next() ),
         FORALL(
            ( item()->start_variable_iterator() ;
              item()->valid_variable() ;
              item()->go_next_variable() ),
            fields->has( item()->variable()->name() ) ) ) ) ;
   PEL_CHECK_POST( !is_valid() ) ;
}

//----------------------------------------------------------------------
PDE_SetOfFieldCompositions:: ~PDE_SetOfFieldCompositions( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
PDE_FieldComposition*
PDE_SetOfFieldCompositions:: item( std::string const& name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfFieldCompositions:: item( std::string const& )" ) ;
   PEL_CHECK_PRE( !name.empty() ) ;

   size_t const idx = find_index( name ) ;
   if( idx>=COMPO_TABLE->count() )
   {
      std::string msg ;
      msg += "*** access to : \""+name+"\"\n" ;
      msg += "    unknown composition of fields" ;
      PEL_Error::object()->raise_plain( msg ) ;
   }

   PDE_FieldComposition* result =
                  static_cast<PDE_FieldComposition*>( COMPO_TABLE->at(idx) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->name() == name ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_SetOfFieldCompositions:: has( std::string const& name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfFieldCompositions:: has" ) ;
   PEL_CHECK_PRE( !name.empty() ) ;

   return( find_index( name )<nb_compositions() ) ;
}

//----------------------------------------------------------------------
size_t
PDE_SetOfFieldCompositions:: nb_compositions( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfFieldCompositions:: nb_compositions" ) ;
   return( COMPO_TABLE->count() ) ;
}

//----------------------------------------------------------------------
void
PDE_SetOfFieldCompositions:: start( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfFieldCompositions:: start" ) ;
   INDEX = 0 ;
}

//----------------------------------------------------------------------
void
PDE_SetOfFieldCompositions:: go_next( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfFieldCompositions:: go_next" ) ;
   INDEX++ ;
}

//----------------------------------------------------------------------
bool
PDE_SetOfFieldCompositions:: is_valid( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfFieldCompositions:: is_valid" ) ;
   return( INDEX<nb_compositions() ) ;
}

//----------------------------------------------------------------------
PDE_FieldComposition*
PDE_SetOfFieldCompositions:: item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfFieldCompositions:: item()" ) ;
   PEL_CHECK_PRE( is_valid() ) ;

   PDE_FieldComposition* result = 
          static_cast<PDE_FieldComposition*>( COMPO_TABLE->at( INDEX ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_SetOfFieldCompositions:: print( std::ostream& os,
                                    size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string space( indent_width, ' ' ) ;
   os << space << nb_compositions() << " composition of fields : "
      << std::endl ;
   for( start() ; is_valid() ; go_next() )
   {
      os << space << "   - \"" << item()->name() << "\"" << std::endl ;
   }
}

//----------------------------------------------------------------------
void
PDE_SetOfFieldCompositions:: append( PDE_FieldComposition* composition )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfFieldCompositions:: append" ) ;
   PEL_CHECK( composition!= 0 ) ;
   PEL_CHECK( composition->owner()==0 ) ;
   PEL_CHECK( !has( composition->name() ) ) ;
   PEL_SAVEOLD( size_t, nb_compositions, nb_compositions() ) ;

   COMPO_TABLE->append( composition ) ;
   composition->set_owner( this ) ;
   
   PEL_CHECK_POST( nb_compositions() == OLD( nb_compositions )+1 ) ;
   PEL_CHECK_POST( has( composition->name() ) ) ;
   PEL_CHECK_POST( item( composition->name() ) == composition ) ;
   PEL_CHECK_POST( composition->is_under_ownership_of( this ) ) ;
}

//----------------------------------------------------------------------
size_t
PDE_SetOfFieldCompositions:: find_index( std::string const& name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfFieldCompositions:: find_index" ) ;

   size_t result = PEL::bad_index() ;
   for( size_t i=0 ; result==PEL::bad_index() && i<COMPO_TABLE->count() ; ++i )
   {
      PDE_FieldComposition* compo =
                static_cast<PDE_FieldComposition*>( COMPO_TABLE->at(i) ) ;
      if( compo->name()==name )
      {
         result = i ;
      }
   }
   return( result ) ;
}
