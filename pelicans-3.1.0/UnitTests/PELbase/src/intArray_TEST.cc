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

#include <intArray_TEST.hh>

#include <PEL.hh>
#include <PEL_System.hh>
#include <intArray3D.hh>

#include <iostream>

using std::cout ;
using std::endl ;

intArray_TEST const* 
intArray_TEST:: REGISTRATOR = new intArray_TEST() ;

//-------------------------------------------------------------------------
intArray_TEST:: intArray_TEST( void )
//-------------------------------------------------------------------------
   :   PEL_ObjectTest( "intArray", "intArray_TEST" )
{
}

//-------------------------------------------------------------------------
intArray_TEST:: ~intArray_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
intArray_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   do_tests_3D() ;
}

//-------------------------------------------------------------------------
void
intArray_TEST:: do_tests_3D( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "intArray_TEST:: do_tests_3D" ) ;

   bool ok = true ;
   
   intArray3D xxx( 3, 1, 2 ) ;
   int val = -20 ;
   for( size_t i=0 ; i<3 ; ++i )
      for( size_t j=0 ; j<1 ; ++j )
         for( size_t k=0 ; k<2 ; ++k )
            xxx( i, j, k ) = val++ ; 
   
   intArray3D yyy( 12, 8, 9 ) ;
   
   yyy = xxx ;
   ok &= ( xxx.index_bound( 0 ) == 3 ) ;
   ok &= ( xxx.index_bound( 1 ) == 1 ) ;
   ok &= ( xxx.index_bound( 2 ) == 2 ) ;
   for( size_t i=0 ; i<3 ; ++i )
      for( size_t j=0 ; j<1 ; ++j )
         for( size_t k=0 ; k<2 ; ++k )
            ok &= ( yyy( i, j, k ) == xxx( i, j, k ) ) ;
   
   notify_one_test_result( "intArray3D::operator=", ok ) ;
}

