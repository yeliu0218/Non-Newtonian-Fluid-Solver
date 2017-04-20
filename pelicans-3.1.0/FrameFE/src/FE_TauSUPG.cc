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

#include <FE_TauSUPG.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <iostream>

using std::string ;

//-----------------------------------------------------------------------------
FE_TauSUPG:: FE_TauSUPG( PEL_Object* a_owner, 
                         PEL_ModuleExplorer const* exp )
//-----------------------------------------------------------------------------
   : FE_TauStab( a_owner, exp )
   , COEF( 1.0 )
{
   PEL_LABEL( "FE_TauSUPG:: FE_TauSUPG" ) ;
   PEL_CHECK( exp != 0 ) ;

   COEF = exp->double_data( "upwind_factor" ) ;
}

//-----------------------------------------------------------------------------
FE_TauSUPG:: ~FE_TauSUPG( void )
//-----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
double
FE_TauSUPG:: tau( double h, double alpha, double normv, double mu ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_TauSUPG:: tau" ) ;
   PEL_CHECK_PRE( tau_PRE( h, alpha, normv, mu ) ) ;

   double result = 0.0 ;

   if( normv > 1.0e-20 )
   {
      double Pe = normv * h / ( 2. * mu ) ;
      result = ( h / ( 2. * normv ) ) * dzeta( Pe ) ;
   }
   else
   {
      result = h * h / ( 4.* mu ) ;
   }

   result *= COEF ;
   return( result ) ;
}
