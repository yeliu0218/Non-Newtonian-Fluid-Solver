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

#include <FE_AdaptationStepCHARMStest.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <PDE_DomainAndFields.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_AdapterCHARMS.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <ios>
#include <iostream>
#include <iomanip>

using std::cout ; using std::endl ;
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;

FE_AdaptationStepCHARMStest const* 
FE_AdaptationStepCHARMStest::PROTOTYPE = new FE_AdaptationStepCHARMStest() ;

//----------------------------------------------------------------------
FE_AdaptationStepCHARMStest:: FE_AdaptationStepCHARMStest( void )
//----------------------------------------------------------------------
   : FE_OneStepIteration( "FE_AdaptationStepCHARMStest" )
{
}

//----------------------------------------------------------------------
FE_AdaptationStepCHARMStest*
FE_AdaptationStepCHARMStest:: create_replica( PEL_Object* a_owner,
                                              PDE_DomainAndFields const* dom,
                                              FE_SetOfParameters const* prms,
                                              PEL_ModuleExplorer* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationStepCHARMStest:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_AdaptationStepCHARMStest* result = 
                new FE_AdaptationStepCHARMStest( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_AdaptationStepCHARMStest:: FE_AdaptationStepCHARMStest( 
                                           PEL_Object* a_owner,
                                           PDE_DomainAndFields const* dom,
                                           FE_SetOfParameters const* prms,
                                           PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , DFS( dom->set_of_discrete_fields() )
   , EXP( exp->create_subexplorer( this, "list_of_PDE_DiscreteField" ) )
   , DA( dom->adapter_CHARMS() )
{
}

//----------------------------------------------------------------------
FE_AdaptationStepCHARMStest:: ~FE_AdaptationStepCHARMStest( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
FE_AdaptationStepCHARMStest:: do_one_inner_iteration( 
                                           FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationStepCHARMStest:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "FE_AdaptationStepCHARMStest:: do_one_inner_iteration" );
   
   EXP->start_module_iterator() ;
   for( ; EXP->is_valid_module() ; EXP->go_next_module() )
   {
      PEL_ModuleExplorer* se = EXP->create_subexplorer( 0 ) ;
      PDE_DiscreteField* sf = DFS->item( se->string_data( "current" ) ) ;
         
      for( size_t n=0 ; n<sf->nb_nodes() ; ++n )
      {
         if( sf->DOF_value( 0, n ) == PEL::bad_double() )
         {
            sf->set_DOF_value( 0, n, -3.0 ) ;
         }
      }
         
      se->destroy() ; se = 0 ;
   }

   stop_total_timer() ;
}
