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

#include <FE_AdaptationIndicator.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_assertions.hh>

#include <FE_SetOfParameters.hh>

#include <iostream>
#include <string>
#include <sstream>

using std::string ;
using std::cout ; using std::endl ;
using std::ostringstream ;

struct FE_AdaptationIndicator_ERROR
{
   static void n0( void ) ;
} ;

//---------------------------------------------------------------------------
FE_AdaptationIndicator*
FE_AdaptationIndicator:: make( PEL_Object* a_owner,
                               PDE_DomainAndFields const* dom,
                               FE_SetOfParameters const* prms,
                               PEL_ModuleExplorer const* exp,
                               size_t a_verbose_level  )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationIndicator:: make" ) ;
   PEL_CHECK_PRE( dom != 0 ) ;
   PEL_CHECK_PRE( prms != 0 ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   
   string name = exp->string_data( "concrete_name" ) ;
   FE_AdaptationIndicator const* proto =
   static_cast<FE_AdaptationIndicator const*>( plugins_map()->item( name ) ) ;
      
   FE_AdaptationIndicator* result = proto->create_replica( a_owner, 
                                                           dom,
                                                           prms,
                                                           exp,
                                                           a_verbose_level ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
FE_AdaptationIndicator:: FE_AdaptationIndicator( PEL_Object* a_owner,
                                                 size_t a_verbose_level  )
//---------------------------------------------------------------------------
   : PDE_AdaptationIndicator( a_owner, a_verbose_level )
{
}

//-------------------------------------------------------------------------
FE_AdaptationIndicator:: FE_AdaptationIndicator( std::string const& a_name )
//-------------------------------------------------------------------------
   : PDE_AdaptationIndicator( "" )
{
   PEL_LABEL( "FE_AdaptationIndicator:: FE_AdaptationIndicator" ) ;

   plugins_map()->register_item( a_name, this ) ;

   PEL_CHECK_POST( is_a_prototype() ) ;
}

//---------------------------------------------------------------------------
FE_AdaptationIndicator:: ~FE_AdaptationIndicator( void  )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
FE_AdaptationIndicator:: set_time_iterator( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationIndicator:: set_time_iterator" ) ;
   PEL_CHECK_PRE( t_it != 0 ) ;
   
   T_IT = t_it ;
}

//---------------------------------------------------------------------------
void
FE_AdaptationIndicator:: reset( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationIndicator:: reset" ) ;   
   
   // to be able to check that set_time_iterator() is called
   // before build()
   T_IT = 0 ;
   
   reset_self() ;
}

//---------------------------------------------------------------------------
void
FE_AdaptationIndicator:: build( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationIndicator:: build" ) ;   
   
   if( T_IT == 0 ) FE_AdaptationIndicator_ERROR::n0() ;
   
   build_self( T_IT ) ;
}

//---------------------------------------------------------------------------
PDE_AdaptationIndicator*
FE_AdaptationIndicator:: create_replica( 
                                  PEL_Object* a_owner,
                                  PDE_DomainAndFields const* dom,
                                  PEL_ModuleExplorer const* exp,
                                  size_t a_verbose_level ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationIndicator:: create_replica" ) ;
   PEL_ASSERT( false ) ;
   PDE_AdaptationIndicator* result = 0 ;
   return( result ) ;
}

//---------------------------------------------------------------------------
bool
FE_AdaptationIndicator:: create_replica_PRE( PEL_Object* a_owner,
                                             PDE_DomainAndFields const* dom,
                                             FE_SetOfParameters const* prms,
                                             PEL_ModuleExplorer const* exp,
                                             size_t a_verbose_level ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( dom != 0 ) ;
   PEL_ASSERT( prms != 0 ) ;
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
FE_AdaptationIndicator:: create_replica_POST(  
                                       FE_AdaptationIndicator const* result,
                                       PEL_Object* a_owner,
                                       PDE_DomainAndFields const* dom,
                                       FE_SetOfParameters const* prms,
                                       PEL_ModuleExplorer const* exp,
                                       size_t a_verbose_level ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( result->verbose_level() == a_verbose_level ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
PEL_ObjectRegister*
FE_AdaptationIndicator:: plugins_map( void )
//---------------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
          PEL_ObjectRegister::create( PEL_Root::object(),
                                      "FE_AdaptationIndicator descendant" ) ;
   return( result ) ;
}

//internal-------------------------------------------------------------------
void
FE_AdaptationIndicator_ERROR:: n0( void )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "*** FE_AdaptationIndicator error:" << endl << endl ;
   mesg << "    \"set_time_iterator\" must be called before \"build\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

