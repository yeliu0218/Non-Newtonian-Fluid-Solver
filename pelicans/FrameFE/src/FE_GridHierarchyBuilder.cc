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

#include <FE_GridHierarchyBuilder.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Timer.hh>

#include <PDE_AdapterCHARMS.hh>
#include <PDE_DomainAndFields.hh>

#include <FE_TimeIterator.hh>

#include <iostream>
#include <string>

using std::endl ;

FE_GridHierarchyBuilder const*
FE_GridHierarchyBuilder::PROTOTYPE = new FE_GridHierarchyBuilder() ;

//----------------------------------------------------------------------
FE_GridHierarchyBuilder:: FE_GridHierarchyBuilder( void )
//----------------------------------------------------------------------
   : FE_OneStepIteration( "FE_GridHierarchyBuilder" ) 
{
}

//----------------------------------------------------------------------
FE_GridHierarchyBuilder*
FE_GridHierarchyBuilder:: create_replica( PEL_Object* a_owner,
                                          PDE_DomainAndFields const* dom,
                                          FE_SetOfParameters const* prms,
                                          PEL_ModuleExplorer* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_GridHierarchyBuilder:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_GridHierarchyBuilder* result = 
                new FE_GridHierarchyBuilder( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_GridHierarchyBuilder:: FE_GridHierarchyBuilder( 
                                          PEL_Object* a_owner,
                                          PDE_DomainAndFields const* dom,
                                          FE_SetOfParameters const* prms,
                                          PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , DOM( dom )
   , DA( 0 )
{
   PEL_LABEL( "FE_GridHierarchyBuilder:: FE_GridHierarchyBuilder" ) ;

   DA = dom->adapter_CHARMS() ;
   if( DA == 0 ) 
      PEL_Error::object()->raise_plain( "discretization adapter requested" ) ;
}

//----------------------------------------------------------------------
FE_GridHierarchyBuilder:: ~FE_GridHierarchyBuilder( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
FE_GridHierarchyBuilder:: do_before_time_stepping( 
                                         FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_GridHierarchyBuilder:: do_before_time_stepping" ) ;

   start_total_timer( "FE_GridHierarchyBuilder:: do_before_time_stepping" ) ;
   // --------------

   DA->reset() ;

   bool keep_adapting = false ;
   do 
   {      
      if( verbose_level() > 1 ) 
         PEL::out() << indent() << "adaptation..." << endl ;
      
      DA->adapt() ;

      keep_adapting = DA->something_changed() ;

      if( keep_adapting )
      {
         DOM->apply_requests_of_DOFs_values_modules( true ) ;
      }
   } while( keep_adapting ) ;

   stop_total_timer() ;
   // -------------
}

//---------------------------------------------------------------------------
void 
FE_GridHierarchyBuilder:: print_additional_times( std::ostream& os,
                                                  size_t indent_width ) const 
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
FE_GridHierarchyBuilder:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_GridHierarchyBuilder:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;
}
