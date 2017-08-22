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

#include <CH_BulkChemicalPotential.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Object.hh>
#include <PEL_assertions.hh>

#include <CH_BulkEnergy.hh>

#include <iostream>

//----------------------------------------------------------------------
CH_BulkChemicalPotential*
CH_BulkChemicalPotential:: create( PEL_Object* a_owner, 
                                   PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_BulkChemicalPotential:: create" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   CH_BulkChemicalPotential* result = new CH_BulkChemicalPotential( a_owner, exp ) ;
   
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
CH_BulkChemicalPotential:: CH_BulkChemicalPotential( 
                                          PEL_Object* a_owner,
                                          PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , S1( exp->double_data( "coef_Sigma_1" ) )
   , S2( exp->double_data( "coef_Sigma_2" ) )
   , S3( exp->double_data( "coef_Sigma_3" ) )
   , F0( 0 )
   , P( 0 )
{
   PEL_LABEL( "CH_BulkChemicalPotential:: CH_BulkChemicalPotential" ) ;
   
   ST = 3.0*S1*S2*S3/(S1*S2+S1*S3+S2*S3) ;
   
   PEL_ModuleExplorer* e = exp->create_subexplorer( 0, "CH_BulkEnergy#F0" ) ;
   F0 = CH_BulkEnergy::make( this, S1, S2, S3, e ) ;
   e->destroy() ; e = 0 ;
   
   if( exp->has_module( "CH_BulkEnergy#P" ) )
   {
      e = exp->create_subexplorer( 0, "CH_BulkEnergy#P" ) ;
      P = CH_BulkEnergy::make( this, S1, S2, S3, e ) ;
      e->destroy() ; e = 0 ;
   }
}

//----------------------------------------------------------------------
CH_BulkChemicalPotential:: ~CH_BulkChemicalPotential( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
double 
CH_BulkChemicalPotential:: Sigma1( void ) const
//----------------------------------------------------------------------
{
   return( S1 ) ; 
}

//----------------------------------------------------------------------
double 
CH_BulkChemicalPotential:: Sigma2( void ) const
//----------------------------------------------------------------------
{
   return( S2 ) ; 
}

//----------------------------------------------------------------------
double 
CH_BulkChemicalPotential:: Sigma3( void ) const
//----------------------------------------------------------------------
{
   return( S3 ) ; 
}

//----------------------------------------------------------------------
double 
CH_BulkChemicalPotential:: F( double c1, double c2, double c3 ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_BulkChemicalPotential:: F" ) ;
   
   double result = F0->F( c1, c2, c3 ) ;
   if( P != 0 )
   {
      result += P->F( c1, c2, c3 ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double 
CH_BulkChemicalPotential:: DDiF( double c1, double c2, 
                                double c1_exp, double c2_exp, 
                                size_t i, double eps ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_BulkChemicalPotential:: DDiF" ) ;
   PEL_CHECK_PRE( i==1 || i==2 ) ;
   PEL_CHECK_PRE( eps > 0.0 ) ;

   double result = F0->DDiF( c1, c2, c1_exp, c2_exp, i, eps ) ;
   if( P != 0 ) 
   {
      result+= P->DDiF( c1, c2, c1_exp, c2_exp, i, eps ) ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
double
CH_BulkChemicalPotential:: dj_DDiF( double c1, double c2,
                                   double c1_exp, double c2_exp, 
                                   size_t i, size_t j, double eps ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "CH_BulkChemicalPotential:: dj_DDiF" ) ;
   PEL_CHECK_PRE( i==1 || i==2 ) ;
   PEL_CHECK_PRE( j==1 || j==2 ) ;
   PEL_CHECK_PRE( eps > 0.0 ) ;

   double result = F0->dj_DDiF( c1, c2, c1_exp, c2_exp, i, j, eps ) ;
   if( P != 0 ) 
   {
      result+= P->dj_DDiF( c1, c2, c1_exp, c2_exp, i, j, eps ) ;
   }

   return( result ) ;
}

