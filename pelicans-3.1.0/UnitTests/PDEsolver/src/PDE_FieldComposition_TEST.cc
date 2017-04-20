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

#include <PDE_FieldComposition_TEST.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_FieldComposition.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SetOfFieldCompositions.hh>

#include <iostream>

PDE_FieldComposition_TEST*
PDE_FieldComposition_TEST:: UNIQUE_INSTANCE = new PDE_FieldComposition_TEST() ;

//---------------------------------------------------------------------------
PDE_FieldComposition_TEST:: PDE_FieldComposition_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "PDE_FieldComposition", "PDE_FieldComposition_TEST" )
{
}

//---------------------------------------------------------------------------
PDE_FieldComposition_TEST:: ~PDE_FieldComposition_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_FieldComposition_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FieldComposition_TEST:: process_one_test" ) ;
   
   PEL_ModuleExplorer* se =
                        exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   PDE_DomainAndFields* dom = PDE_DomainAndFields::create( 0, se ) ;
   se->destroy() ; se = 0 ;

   PDE_SetOfDiscreteFields const* dfs =  dom->set_of_discrete_fields() ;
   PDE_DiscreteField const* fresult = dfs->item( "result" ) ;

   PDE_SetOfFieldCompositions const* cs = dom->set_of_field_compositions() ;
   PDE_FieldComposition* result = cs->item( "result" ) ;
   
   size_t nb_val = exp->int_data( "nb_values" ) ;
   notify_one_test_result( exp->name() + " : nb_values",
                           nb_val <= fresult->nb_nodes() ) ;
   
   double const d_eps = exp->double_data( "dbl_epsilon" ) ;
   double const d_min = exp->double_data( "dbl_minimum" ) ;
   
   bool ok = true ;
   for( size_t i=0 ; i<nb_val ; i++ )
   {
      result->start_variable_iterator() ;
      for( ; result->valid_variable() ; result->go_next_variable() )
      {
         PDE_DiscreteField const* f = result->variable() ;
         for( size_t c=0 ; c<f->nb_components() ; ++c )
         {
            result->set_variable_value( f, c, f->DOF_value( 0, i, c ) ) ;
         }
      }
      
      result->compute() ;
      
      for( size_t c=0 ; c<fresult->nb_components() ; ++c )
      {
         double const v_ref = fresult->DOF_value( 0, i, c ) ;
         double const v = result->value( c ) ;
         bool okt = PEL::double_equality( v, v_ref, d_eps, d_min ) ;
         ok = ok && okt ;
         if( !okt )
         {
            out() << v_ref << " <> " << v << std::endl ;
            if( PEL::abs( v_ref )>d_min )
            {
               out() << "   epsilon : "
                     << PEL::abs( PEL::min( 1.-v/v_ref, 1.+v/v_ref ) )
                     << std::endl ;
            }
         }
      }
   }
   notify_one_test_result( exp->name() + " : values", ok) ;
   
   dom->destroy() ; dom = 0 ;
}

// //---------------------------------------------------------------------------
// void
// PDE_FieldComposition_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
// //---------------------------------------------------------------------------
// {
//    PEL_LABEL( "PDE_FieldComposition_TEST:: process_one_test" ) ;

//    PEL_ModuleExplorer* se = 0 ;

//    se = exp->create_subexplorer( this, "PDE_DomainAndFields" ) ;
//    PDE_DomainAndFields* dom = PDE_DomainAndFields::create( 0, se, 0 ) ;

//    PDE_SetOfDiscreteFields const* dfs =  dom->set_of_discrete_fields() ;
//    PDE_DiscreteField const* fresult = dfs->item( "result" ) ;
   
//    PEL_List * list_of_comps = PEL_List::create( this ) ;
//    PDE_FieldCompositionExp* result = 0 ;
   
//    se = exp->create_subexplorer( this, "field_compositions" ) ;
//    se->start_entry_iterator() ;
//    for( ; se->is_valid_entry() ; se->go_next_entry() )
//    {
//       PDE_FieldCompositionExp* l = 
//                             PDE_FieldCompositionExp::create( list_of_comps, 
//                                                              se->keyword(), 
//                                                              se->data( this ), 
//                                                              dfs ) ;
//       list_of_comps->append( l ) ;
//       if( l->name() == "result" ) result = l ;
//    }
//    PEL_ASSERT( result != 0 ) ;
   
//    PEL_Iterator* it = list_of_comps->create_iterator( 0 ) ;
//    for( it->start() ; it->is_valid() ; it->go_next() )
//    {
//       PDE_FieldCompositionExp* l = 
//                          static_cast<PDE_FieldCompositionExp*>( it->item() ) ;
//       l->do_the_link( list_of_comps ) ;
//    }
//    for( it->start() ; it->is_valid() ; it->go_next() )
//    {
//       PDE_FieldCompositionExp* l = 
//                          static_cast<PDE_FieldCompositionExp*>( it->item() ) ;
//       l->complete_internal_dependencies() ;
//    }
//    it->destroy() ;

//    size_t nb_val = exp->int_data( "nb_values" ) ;
//    notify_one_test_result( exp->name() + " : nb_values",
//                            nb_val <= fresult->nb_nodes() ) ;
   
//    bool ok = true ;
//    for( size_t i=0 ; i<nb_val ; i++ )
//    {
//       result->start_variable_iterator() ;
//       for( ; result->valid_variable() ; result->go_next_variable() )
//       {
//          PDE_DiscreteField const* f = result->variable() ;
//          for( size_t c=0 ; c<f->nb_components() ; ++c )
//          {
//             result->set_variable_value( f, c, f->DOF_value( 0, i, c ) ) ;
//          }
//       }
//       result->compute() ;
      
//       for( size_t c=0 ; c<fresult->nb_components() ; ++c )
//       {
//          bool okt = PEL::equal( fresult->DOF_value( 0, i, c ), 
//                                 result->value( c ) ) ;
//          ok = ok && okt ;
//          if( !okt ) std::cout << fresult->DOF_value( 0, i, c ) << " <> " 
//                               << result->value( c ) << std::endl ;
//       }
//    }
//    notify_one_test_result( exp->name() + " : value", ok) ;
//    dom->destroy() ;
// }

