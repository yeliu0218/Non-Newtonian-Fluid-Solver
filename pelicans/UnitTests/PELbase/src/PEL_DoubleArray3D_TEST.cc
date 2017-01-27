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

#include <PEL_DoubleArray3D_TEST.hh>

#include <PEL.hh>
#include <PEL_DoubleArray3D.hh>
#include <PEL_Map.hh>
#include <doubleArray3D.hh>

//-------------------------------------------------------------------------
PEL_DoubleArray3D_TEST* PEL_DoubleArray3D_TEST::registered_test =
new PEL_DoubleArray3D_TEST() ;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void
PEL_DoubleArray3D_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   size_t const n = 2, m = 3 , p = 5 ;
   
   doubleArray3D d(n,m,p) ;
   for( size_t i=0 ; i<n ; i++ )
      for( size_t j=0 ; j<n ; j++ )
         for( size_t k=0 ; k<p ; k++ )
            d(i,j,k) = i*1.0+j*2.0+k*3.0 ;
   
   PEL_DoubleArray3D const* p_d = PEL_DoubleArray3D::create( this, d ) ;
   bool ok = true ;
   for( size_t i=0 ; i<n ; i++ )
      for( size_t j=0 ; j<n ; j++ )
         for( size_t k=0 ; k<p ; k++ )
            ok = ok && d(i,j,k) == p_d->to_double_array3D()(i,j,k) ;
   notify_one_test_result("  to_double_array3D", ok ) ;
   PEL_DoubleArray3D* another_d = p_d->create_clone( this ) ;
   ok = true ;
   for( size_t i=0 ; i<n ; i++ )
      for( size_t j=0 ; j<n ; j++ )
         for( size_t k=0 ; k<p ; k++ )
            ok = ok && d(i,j,k) == another_d->to_double_array3D()(i,j,k) ;
   notify_one_test_result("  clone", ok ) ;
}



//-------------------------------------------------------------------------
PEL_DoubleArray3D_TEST:: PEL_DoubleArray3D_TEST( void )
//-------------------------------------------------------------------------
   :   PEL_ObjectTest( "PEL_DoubleArray3D", "PEL_DoubleArray3D_TEST" )
{
}



//-------------------------------------------------------------------------
PEL_DoubleArray3D_TEST:: ~PEL_DoubleArray3D_TEST( void )
//-------------------------------------------------------------------------
{
}



