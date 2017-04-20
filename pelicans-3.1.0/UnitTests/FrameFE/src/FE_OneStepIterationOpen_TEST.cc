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

#include <FE_OneStepIterationOpen_TEST.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>

#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_SetOfDomains.hh>

#include <FE.hh>
#include <FE_OneStepIteration.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <FE_OneStepIterationOpen.hh>

#include <ios>
#include <iomanip>
#include <iostream>
#include <sstream>

using std::cout ; using std::endl ;
using std::ios_base ;
using std::string ; using std::ostringstream ;

FE_OneStepIterationOpen_TEST const* 
FE_OneStepIterationOpen_TEST::PROTOTYPE = new FE_OneStepIterationOpen_TEST() ;

//---------------------------------------------------------------------------
FE_OneStepIterationOpen_TEST:: FE_OneStepIterationOpen_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "FE_OneStepIterationOpen", 
                     "FE_OneStepIterationOpen_TEST" )
{
}

//---------------------------------------------------------------------------
FE_OneStepIterationOpen_TEST:: ~FE_OneStepIterationOpen_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void 
FE_OneStepIterationOpen_TEST:: process_one_test( 
                                              PEL_ModuleExplorer const* exp ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIterationOpen_TEST:: process_one_test" ) ;

   // -- similar to instantiations in FE_StepByStepProgression
   FE::set_geometry( FE::cartesian ) ;

   PEL_ModuleExplorer* ee = 
                       exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   PDE_DomainAndFields* dom = PDE_DomainAndFields::create( 0, ee, PEL_Exec::communicator() ) ;
   size_t nb_dims = dom->nb_space_dimensions() ;
   ee->destroy() ;

   ee = exp->create_subexplorer( 0, "FE_SetOfParameters" ) ;
   FE_SetOfParameters* prms = FE_SetOfParameters::create( dom, dom, ee ) ;
   ee->destroy() ;

   ee = exp->create_subexplorer( 0, "FE_OneStepIteration" ) ;
   FE_OneStepIteration* one_it = 
                        FE_OneStepIteration::make( dom, dom, prms, ee ) ;
   ee->destroy() ;

   ee = exp->create_subexplorer( 0, "FE_TimeIterator" ) ;
   FE_TimeIterator* t_it = FE_TimeIterator::create( 0, ee ) ;
   ee->destroy() ; ee = 0 ;

   ee = exp->create_subexplorer( 0, "jacobian_test" ) ;
   FE_OneStepIterationOpen* one_it_jac = FE_OneStepIterationOpen::object( 
                    ee->string_data( "one_step_iteration_with_jacobian" ) ) ;
   double d_eps = ee->double_data( "dbl_epsilon" ) ;
   double d_min = ee->double_data( "dbl_minimum" ) ;
   double hh = ee->double_data( "hh" ) ;
   ee->destroy() ; ee = 0 ;

   one_it->do_before_time_stepping( t_it ) ;
   t_it->start() ;

   for( size_t i=0 ; i<one_it_jac->nb_unknowns() ; ++i )
   {
      for( size_t j=0 ; j<one_it_jac->nb_unknowns() ; ++j )
      {
         PDE_DiscreteField* uu = one_it_jac->field( j ) ;
         size_t l_uu = one_it_jac->level_of_field( j ) ;
         PDE_LinkDOF2Unknown const* uu_link =
                                    one_it_jac->link_DOF_2_unknown( j ) ;

         one_it_jac->build_function_and_jacobian( t_it ) ;

         LA_SeqMatrix const* jac  = one_it_jac->create_jacobian( 0, i, j ) ;
         LA_SeqVector const* func = one_it_jac->create_function( 0, i ) ;

         bool eq = false ;
         bool ok = true ;
         for( size_t n=0 ; n<uu->nb_nodes() ; ++n )
         {
            for( size_t ic=0 ; ic<uu->nb_components() ; ++ic )
            {
               if( uu_link->DOF_is_unknown( n, ic ) )
               {
                  size_t i_unk = uu_link->unknown_linked_to_DOF( n, ic ) ;
                  double old_val = uu->DOF_value( l_uu, n, ic ) ;
                  uu->set_DOF_value( l_uu, n, old_val + hh , ic ) ;

                  one_it_jac->build_function_and_jacobian( t_it ) ;
                  LA_SeqVector const* new_func = 
                        one_it_jac->create_function( 0, i) ;
                  PEL_ASSERT( new_func->nb_rows() == func->nb_rows() ) ;
                  for( size_t ii=0 ; ii<func->nb_rows() ; ++ii )
                  {
                     double newf = new_func->item( ii ) ;
                     double oldf = func->item( ii ) ;
                     double xx = - ( newf - oldf )/ hh ;
                     eq = PEL::double_equality( xx, jac->item( ii, i_unk ),
                                                d_eps, d_min ) ;
                     if( !eq ) 
                     {
                        display_error( i, j, 
                                       xx, jac->item( ii, i_unk ), 
                                       newf, oldf ) ;
                     }
                     ok = ok && eq ;
                  }

                  new_func->destroy() ;
                  uu->set_DOF_value( l_uu, n, old_val, ic ) ;
               }
            }
         }
         {
            ostringstream mesg ;
            mesg << " block(" << i << "," << j << ") (" << nb_dims << "D)" ;
            notify_one_test_result( one_it_jac->name() + mesg.str(), ok ) ;
         }
         jac->destroy() ;
         func->destroy() ;
      }
   }
   dom->destroy() ;
   t_it->destroy() ;
}

//----------------------------------------------------------------------------
void
FE_OneStepIterationOpen_TEST:: display_error( size_t i, size_t j,
                                              double xx_1, double xx_2,
                                              double newf, double oldf ) const
//----------------------------------------------------------------------------
{
   ios_base::fmtflags original_flags = out().flags() ;
   out().setf( ios_base::uppercase | ios_base::scientific ) ;
   out() << std::setprecision( 10 ) ;

   out() << "block (" << i << "," << j << ")" << endl ;
   out() << std::setw( 20 ) << xx_1
         << std::setw( 20 ) << xx_2
         << std::setw( 20 ) << newf
         << std::setw( 20 ) << oldf
         << endl ;

   out().flags( original_flags ) ;
}
