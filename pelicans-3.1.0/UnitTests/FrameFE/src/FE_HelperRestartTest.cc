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

#include <FE_HelperRestartTest.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <FE_TimeIterator.hh>

#include <iostream>

FE_HelperRestartTest const*
FE_HelperRestartTest::PROTOTYPE = new FE_HelperRestartTest() ;

//---------------------------------------------------------------------------
FE_HelperRestartTest:: FE_HelperRestartTest( void )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( "FE_HelperRestartTest" )
   , FIELD( 0 )
   , FIELD_LEVEL( PEL::bad_index() )
   , COEF( PEL::bad_double() )
{
}

//---------------------------------------------------------------------------
FE_HelperRestartTest*
FE_HelperRestartTest:: create_replica(
                             PEL_Object* a_owner,
                             PDE_DomainAndFields const* dom,
                             FE_SetOfParameters const* prms,
                             PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_HelperRestartTest:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_HelperRestartTest* result =
      new FE_HelperRestartTest( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
FE_HelperRestartTest:: FE_HelperRestartTest(
                             PEL_Object* a_owner,
                             PDE_DomainAndFields const* dom,
                             FE_SetOfParameters const* prms,
                             PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , FIELD(
      dom->set_of_discrete_fields()->item( exp->string_data( "field_name" ) ) )
   , FIELD_LEVEL( (size_t) exp->int_data( "field_level" ) )
   , COEF( -3.14 ) 
{
   PEL_LABEL( "FE_HelperRestartTest:: FE_HelperRestartTest" ) ;

   check_field_storage_depth( FIELD, FIELD_LEVEL ) ;
}

//---------------------------------------------------------------------------
FE_HelperRestartTest:: ~FE_HelperRestartTest( void )
//---------------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//---------------------------------------------------------------------------
void
FE_HelperRestartTest:: do_before_time_stepping( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_HelperRestartTest:: do_before_time_stepping" ) ;
   
   // executed even when there is a restart
   COEF = 2. ;
   
   bool restart = ( t_it->time() != t_it->initial_time() ) ;
   if( !restart )
   {
      modify_field( t_it->time_step()/2.0 ) ;
   }
}

//---------------------------------------------------------------------------
void
FE_HelperRestartTest:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_HelperRestartTest:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   modify_field( COEF * t_it->time_step() ) ;
}

//----------------------------------------------------------------------------
void
FE_HelperRestartTest:: modify_field( double xx )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_HelperRestartTest:: modify_field" ) ;
   
   MAX_VAL = PEL::bad_double() ;
   MIN_VAL = PEL::bad_double() ;
   for( size_t n=0 ; n<FIELD->nb_nodes() ; ++n )
   {
      for( size_t ic=0 ; ic<FIELD->nb_components() ; ++ic )
      {
         double new_val = FIELD->DOF_value( FIELD_LEVEL, n, ic ) + xx ;
         FIELD->set_DOF_value( FIELD_LEVEL, n, new_val, ic ) ;
         if( ( MAX_VAL == PEL::bad_double() ) || ( new_val > MAX_VAL ) )
         {
            MAX_VAL = new_val ;
         }
         if( ( MIN_VAL == PEL::bad_double() ) || ( new_val < MIN_VAL ) )
         {
            MIN_VAL = new_val ;
         }
      }
   }
   
   PEL::out() << indent() << "   max(" << FIELD->name() << ") = " 
              << MAX_VAL << std::endl ;;
   PEL::out() << indent() << "   min(" << FIELD->name() << ") = " 
              << MIN_VAL << std::endl ;;
}

