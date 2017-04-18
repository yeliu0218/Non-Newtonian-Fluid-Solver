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

#include <PEL_DoubleComparatorFloat.hh>

#include <PEL_assertions.hh>
#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>

PEL_DoubleComparator const*
PEL_DoubleComparatorFloat:: PROTOTYPE = new PEL_DoubleComparatorFloat() ;

//----------------------------------------------------------------------
PEL_DoubleComparator const*
PEL_DoubleComparatorFloat:: create( PEL_Object* a_owner,
                                    double a_dbl_min )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DoubleComparatorFloat:: create" ) ;
   PEL_CHECK_PRE( a_dbl_min >= 0. ) ;

   PEL_DoubleComparator const* result =
                  new PEL_DoubleComparatorFloat( a_owner, a_dbl_min ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_DoubleComparator const*
PEL_DoubleComparatorFloat:: create_replica(
              PEL_Object* a_owner, PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DoubleComparatorFloat:: create_replica" ) ;
   PEL_CHECK_PRE( create_replica_PRE( a_owner, exp ) ) ;

   double a_dbl_min = exp->double_data( "dbl_minimum" ) ;
   exp->test_data( "dbl_minimum", "dbl_minimum>=0." ) ;
   
   PEL_DoubleComparator const* result =
                   new PEL_DoubleComparatorFloat( a_owner, a_dbl_min ) ;
   
   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_DoubleComparatorFloat:: PEL_DoubleComparatorFloat(
                                PEL_Object* a_owner, double a_dbl_min  )
//----------------------------------------------------------------------
   : PEL_DoubleComparator( a_owner )
   , EPSILON( a_dbl_min )
{
}

//----------------------------------------------------------------------
PEL_DoubleComparatorFloat:: PEL_DoubleComparatorFloat( void )
//----------------------------------------------------------------------
   : PEL_DoubleComparator( "PEL_DoubleComparatorFloat" )
   , EPSILON( PEL::bad_double() )
{
}

//----------------------------------------------------------------------
PEL_DoubleComparatorFloat:: ~PEL_DoubleComparatorFloat( void )
//----------------------------------------------------------------------
{
   if( is_a_prototype() ) PROTOTYPE = 0 ;
}

//----------------------------------------------------------------------
int
PEL_DoubleComparatorFloat:: three_way_comparison( double x,
                                                  double y ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DoubleComparatorFloat:: three_way_comparison" ) ;

   int result = 0 ;
   float const xx = ( PEL::abs(x)>EPSILON ? (float) x : 0. ) ;
   float const yy = ( PEL::abs(y)>EPSILON ? (float) y : 0. ) ;
   float diff =  xx - yy ;
   if( diff > 0. )
   {
      result = 1 ;
   }
   else if( diff < 0. )
   {
      result = -1 ;
   }
   return( result ) ;
}

