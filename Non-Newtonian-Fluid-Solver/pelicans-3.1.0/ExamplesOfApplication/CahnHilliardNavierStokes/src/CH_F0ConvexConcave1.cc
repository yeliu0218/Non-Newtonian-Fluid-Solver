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

#include <CH_F0ConvexConcave1.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

CH_F0ConvexConcave1 const* 
CH_F0ConvexConcave1:: PROTOTYPE = new CH_F0ConvexConcave1() ;

//----------------------------------------------------------------------
CH_F0ConvexConcave1:: CH_F0ConvexConcave1( void )
//----------------------------------------------------------------------
   : CH_BulkEnergy( "CH_F0ConvexConcave1" )
{
}

//----------------------------------------------------------------------
CH_F0ConvexConcave1*
CH_F0ConvexConcave1:: create_replica( PEL_Object* a_owner,
                                      double sigma_1,
                                      double sigma_2,
                                      double sigma_3,
                                      PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_F0ConvexConcave1:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, 
                                  sigma_1, sigma_2, sigma_3, exp ) ) ;

   CH_F0ConvexConcave1* result = new CH_F0ConvexConcave1( a_owner, 
                                                          sigma_1, 
                                                          sigma_2, 
                                                          sigma_3,
                                                          exp ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, 
                                   sigma_1, sigma_2, sigma_3, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
CH_F0ConvexConcave1:: CH_F0ConvexConcave1( PEL_Object* a_owner,
                                           double sigma_1,
                                           double sigma_2,
                                           double sigma_3,
					                            PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : CH_BulkEnergy( a_owner, sigma_1, sigma_2, sigma_3 )
{
}

//----------------------------------------------------------------------
CH_F0ConvexConcave1:: ~CH_F0ConvexConcave1( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
double 
CH_F0ConvexConcave1:: F( double c1, double c2, double c3 ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_F0ConvexConcave1:: F" ) ;

   double result = Sigma1()*c1*c1*(1.0-c1)*(1.0-c1)/2.0
                 + Sigma2()*c2*c2*(1.0-c2)*(1.0-c2)/2.0
                 + Sigma3()*c3*c3*(1.0-c3)*(1.0-c3)/2.0 ;

   return( result ) ;
}

//----------------------------------------------------------------------
double 
CH_F0ConvexConcave1:: diF( double c1, double c2, double c3,
                           double c1_exp, double c2_exp, double c3_exp,
                           size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_F0ConvexConcave1:: diF" ) ;
   PEL_CHECK_PRE( diF_PRE( c1, c2, c3, c1_exp, c2_exp, c3_exp, i ) ) ;

   double result = PEL::bad_double() ;
   if( i==1 )
   { 
      result = 2.*Sigma1()*PEL::pow(c1-0.5,3) - Sigma1()*(2.*c1_exp-1.0)/4. ;
   }
   else if( i==2 )
   {
      result = 2.*Sigma2()*PEL::pow(c2-0.5,3) - Sigma2()*(2.*c2_exp-1.0)/4.  ;

   }
   else if( i==3 )
   {
      result = 2.*Sigma3()*PEL::pow(c3-0.5,3) -Sigma3()*(2.*c3_exp-1.0)/4.  ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double 
CH_F0ConvexConcave1:: DDiF( double c1, double c2, 
                            double c1_exp, double c2_exp, 
                            size_t i, double eps ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_F0ConvexConcave1:: DDiF" ) ;
   PEL_CHECK_PRE( DDiF_PRE( c1, c2, c1_exp, c2_exp, i, eps ) ) ;
   
   double result = PEL::bad_double() ;
   
   double c3 = 1.0 - c1 - c2 ;
   if( i== 1 )
   {
      result = 6./eps * ( Sigma1() * ( 4.*PEL::pow(c1-0.5,3) - c1_exp+0.5 ) -
                          4.*SigmaT() * c1*c2*c3 ) ;
   }
   else if( i == 2 )
   {
      result = 6./eps * ( Sigma2() * ( 4.*PEL::pow(c2-0.5,3) - c2_exp+0.5 ) -
                          4.*SigmaT() * c1*c2*c3 ) ;
   }
   
   return( result ) ;
}

//----------------------------------------------------------------------
double 
CH_F0ConvexConcave1:: dj_DDiF( double c1, double c2, 
                               double c1_exp, double c2_exp, 
                               size_t i, size_t j, double eps ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_F0ConvexConcave1:: dj_DDiF" ) ;
   PEL_CHECK_PRE( dj_DDiF_PRE( c1, c2, c1_exp, c2_exp, i, j, eps ) ) ;
   
   double result = PEL::bad_double() ;
   
   double c3 = 1.0 - c1 - c2 ;
   if( i == 1 )
   {
      if( j == 1 )
      {
         result = 24./eps * ( Sigma1() * 3. * (c1-0.5)*(c1-0.5)
                              - SigmaT() * ( c2*c3-c1*c2 ) ) ;
      }
      else if( j == 2 )
      {
         result = -24./eps * SigmaT() * ( c1*c3 - c1*c2 ) ;
      }
   }
   else if( i == 2 )
   {
      if( j == 2 )
      {
         result = 24./eps * ( Sigma2() * 3. * (c2-0.5)*(c2-0.5)
                              - SigmaT() * ( c1*c3-c1*c2 ) ) ;
      }
      else if( j == 1 )
      {
         result = -24./eps * SigmaT() * ( c2*c3 - c1*c2 ) ;
      }
   }
   
   return( result ) ;
}
