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

#include <PEL_DoubleComparator_TEST.hh>

#include <PEL_DoubleComparator.hh>
#include <PEL_ModuleExplorer.hh>
#include <stringVector.hh>

PEL_ObjectTest const*
PEL_DoubleComparator_TEST:: PROTOTYPE = new PEL_DoubleComparator_TEST() ;

//-------------------------------------------------------------------------
PEL_DoubleComparator_TEST:: PEL_DoubleComparator_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "PEL_DoubleComparator", "PEL_DoubleComparator_TEST" )
{
}

//-------------------------------------------------------------------------
PEL_DoubleComparator_TEST:: ~PEL_DoubleComparator_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
PEL_DoubleComparator_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DoubleComparator_TEST:: process_one_test" ) ;

   // Building the comparator:
   PEL_ModuleExplorer const* sexp =
                       exp->create_subexplorer( 0, "PEL_DoubleComparator" ) ;
   PEL_DoubleComparator const* comp = PEL_DoubleComparator::make( 0, sexp ) ;
   sexp->destroy() ; sexp = 0 ;

   // Processing tests:
   PEL_ModuleExplorer* texp = exp->create_subexplorer( 0, "Tests" ) ;
   for( texp->start_module_iterator() ;
        texp->is_valid_module() ;
        texp->go_next_module() )
   {
      PEL_ModuleExplorer const* ee = texp->create_subexplorer( 0 ) ;
      double x = ee->double_data( "x" ) ;
      double y = ee->double_data( "y" ) ;
      int res = ee->int_data( "three_way_comparison" ) ;
      ee->test_data_in( "three_way_comparison", "0,-1,1" ) ;
      notify_one_test_result( exp->name()+"/"+ee->name(),
                              res == comp->three_way_comparison( x, y ) ) ;
      ee->destroy() ; ee = 0 ;
   }
   texp->destroy() ; texp = 0 ;

   comp->destroy() ; comp = 0 ;
}
