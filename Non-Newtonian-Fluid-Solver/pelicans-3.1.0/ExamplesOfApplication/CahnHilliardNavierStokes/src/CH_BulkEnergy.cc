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

#include <CH_BulkEnergy.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Object.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_assertions.hh>

#include <iostream>

using std::string ;

//----------------------------------------------------------------------
CH_BulkEnergy*
CH_BulkEnergy:: make( PEL_Object* a_owner,
                      double sigma_1,
                      double sigma_2,
                      double sigma_3,
                      PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_BulkEnergy:: make" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   string name = exp->string_data( "concrete_name" ) ;
   CH_BulkEnergy const* proto =
      static_cast<CH_BulkEnergy const*>( plugins_map()->item( name ) ) ;
      
   CH_BulkEnergy* result = proto->create_replica( a_owner, 
                                                  sigma_1,
                                                  sigma_2,
                                                  sigma_3,
                                                  exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->Sigma1() == sigma_1 ) ;
   PEL_CHECK_POST( result->Sigma2() == sigma_2 ) ;
   PEL_CHECK_POST( result->Sigma3() == sigma_3 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
CH_BulkEnergy:: ~CH_BulkEnergy( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
CH_BulkEnergy:: CH_BulkEnergy( std::string const& a_name )
//----------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
{
   PEL_LABEL( "CH_BulkEnergy:: CH_BulkEnergy" ) ;

   plugins_map()->register_item( a_name, this ) ;

   PEL_CHECK_POST( is_a_prototype() ) ;
}

//----------------------------------------------------------------------
CH_BulkEnergy:: CH_BulkEnergy( PEL_Object* a_owner,
                               double sigma_1,
                               double sigma_2,
                               double sigma_3 )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
   , S1( sigma_1 )
   , S2( sigma_2 )
   , S3( sigma_3 )
{
   ST = 3.0*S1*S2*S3/(S1*S2+S1*S3+S2*S3) ;
}

//----------------------------------------------------------------------
double 
CH_BulkEnergy:: Sigma1( void ) const
//----------------------------------------------------------------------
{
   return( S1 ) ; 
}

//----------------------------------------------------------------------
double 
CH_BulkEnergy:: Sigma2( void ) const
//----------------------------------------------------------------------
{
   return( S2 ) ; 
}

//----------------------------------------------------------------------
double 
CH_BulkEnergy:: Sigma3( void ) const
//----------------------------------------------------------------------
{
   return( S3 ) ; 
}

//----------------------------------------------------------------------
double 
CH_BulkEnergy:: SigmaT( void ) const
//----------------------------------------------------------------------
{
   return( ST ) ; 
}

//----------------------------------------------------------------------
double 
CH_BulkEnergy:: dj_ddiF( double c1, double c2, 
                         double c1_exp, double c2_exp, 
                         size_t i, size_t j ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_BulkEnergy:: dj_ddiF" ) ;
   PEL_CHECK_PRE( dj_ddiF_PRE( c1, c2, c1_exp, c2_exp, i, j ) ) ;
   
   PEL_Error::object()->raise_not_implemented( this, "dj_ddiF" ) ;
   return( PEL::bad_double() ) ;
}

//----------------------------------------------------------------------
double 
CH_BulkEnergy:: DDiF( double c1, double c2, 
                      double c1_exp, double c2_exp, 
                      size_t i, double eps ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_BulkEnergy:: DDiF" ) ;
   PEL_CHECK_PRE( DDiF_PRE( c1, c2, c1_exp, c2_exp, i, eps ) ) ;

   double c3 = 1.0 - c1 - c2 ;
   double c3_exp = 1.0 - c1_exp - c2_exp ;
   
   double d1F = diF( c1, c2, c3, c1_exp, c2_exp, c3_exp, 1 ) ;
   double d2F = diF( c1, c2, c3, c1_exp, c2_exp, c3_exp, 2 ) ;
   double d3F = diF( c1, c2, c3, c1_exp, c2_exp, c3_exp, 3 ) ;
   
   double result =  PEL::bad_double() ;
   if( i==1 ) 
   {
      result = 4.*ST*( ( d1F - d2F )/S2 + ( d1F - d3F )/S3 )/eps ;
   }
   else if( i==2 )
   {
      result = 4.*ST*( ( d2F - d1F )/S1 + ( d2F - d3F )/S3 )/eps ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
double
CH_BulkEnergy:: dj_DDiF( double c1, double c2,
                         double c1_exp, double c2_exp, 
                         size_t i, size_t j, double eps ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_BulkEnergy:: dj_DDiF" ) ;
   PEL_CHECK_PRE( dj_DDiF_PRE( c1, c2, c1_exp, c2_exp, i, j, eps ) ) ;

   double dj_dd1F = dj_ddiF( c1, c2, c1_exp, c2_exp, 1, j ) ;
   double dj_dd2F = dj_ddiF( c1, c2, c1_exp, c2_exp, 2, j ) ;
   double dj_dd3F = dj_ddiF( c1, c2, c1_exp, c2_exp, 3, j ) ;
   
   double result =  PEL::bad_double() ;

   if( i==1 ) 
   {
      result = 4.*ST/eps*( ( dj_dd1F - dj_dd2F )/S2 + 
                           ( dj_dd1F - dj_dd3F )/S3 ) ;
   }
   else if( i==2 )
   {
      result = 4.*ST/eps*( ( dj_dd2F - dj_dd1F )/S1 + 
                           ( dj_dd2F - dj_dd3F )/S3 ) ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
bool
CH_BulkEnergy:: is_a_prototype( void ) const
//----------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//----------------------------------------------------------------------
bool
CH_BulkEnergy:: diF_PRE( double c1, double c2, double c3,
                         double c1_exp, double c2_exp, double c3_exp,
                         size_t i ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( i==1 || i==2 || i==3 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
CH_BulkEnergy:: DDiF_PRE( double c1, double c2, 
                          double c1_exp, double c2_exp,
                          size_t i, double eps ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( i==1 || i==2 ) ;
   PEL_ASSERT( eps > 0.0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
CH_BulkEnergy:: dj_DDiF_PRE( double c1, double c2, 
                             double c1_exp, double c2_exp,
                             size_t i, size_t j, double eps ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( i==1 || i==2 ) ;
   PEL_ASSERT( j==1 || j==2 ) ;
   PEL_ASSERT( eps > 0.0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
CH_BulkEnergy:: dj_ddiF_PRE( double c1, double c2, 
                             double c1_exp, double c2_exp,
                             size_t i, size_t j ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( i==1 || i==2 || i==3 ) ;
   PEL_ASSERT( j==1 || j==2 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
CH_BulkEnergy:: create_replica_PRE( PEL_Object* a_owner,
                                    double sigma_1,
                                    double sigma_2,
                                    double sigma_3,
                                    PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
CH_BulkEnergy:: create_replica_POST( CH_BulkEnergy const* result,
                                     PEL_Object* a_owner,
                                     double sigma_1,
                                     double sigma_2,
                                     double sigma_3,
                                     PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   return( true ) ;  
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
CH_BulkEnergy:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
          PEL_ObjectRegister::create( PEL_Root::object(),
                                      "CH_BulkEnergy descendant" ) ;
   return( result ) ;
}
