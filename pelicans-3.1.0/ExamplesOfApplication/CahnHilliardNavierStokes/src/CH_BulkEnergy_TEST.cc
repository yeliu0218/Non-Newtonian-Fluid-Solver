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

#include <CH_BulkEnergy_TEST.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_SetOfDomains.hh>
#include <PDE_LocalFEcell.hh>

#include <FE.hh>
#include <FE_OneStepIteration.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <GE_QRprovider.hh>

#include <CH_BulkEnergy.hh>

#include <ios>
#include <iomanip>
#include <iostream>
#include <sstream>

using std::cout ; using std::endl ;
using std::ios_base ;
using std::string ; using std::ostringstream ;

CH_BulkEnergy_TEST const* 
CH_BulkEnergy_TEST::PROTOTYPE = new CH_BulkEnergy_TEST() ;

//---------------------------------------------------------------------------
CH_BulkEnergy_TEST:: CH_BulkEnergy_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "CH_BulkEnergy", 
                     "CH_BulkEnergy_TEST" )
{
}

//---------------------------------------------------------------------------
CH_BulkEnergy_TEST:: ~CH_BulkEnergy_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void 
CH_BulkEnergy_TEST:: process_one_test( PEL_ModuleExplorer const* exp ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_BulkEnergy_TEST:: process_one_test" ) ;
   
   out() << "| ... " << exp->name() << std::endl ;

   D_EPS = exp->double_data( "dbl_epsilon" ) ;
   D_MIN = exp->double_data( "dbl_minimum" ) ;

   THICKNESS = exp->double_data( "thickness" ) ;

   PEL_ModuleExplorer* e = exp->create_subexplorer( 0, "CH_BulkEnergy" ) ;
   S1 = exp->double_data( "coef_Sigma_1" ) ;
   S2 = exp->double_data( "coef_Sigma_2" ) ;
   S3 = exp->double_data( "coef_Sigma_3" ) ;
   CH_BulkEnergy const* energy = CH_BulkEnergy::make( this, S1, S2, S3, e ) ;
   e->destroy() ; e = 0 ;
   
   ST = energy->SigmaT() ;
   
   bool has_dj_ddiF = exp->bool_data( "has_dj_ddiF" ) ;

   HH = exp->double_data( "hh" ) ;

   bool ok_diF = true ;
   bool ok_DDiF = true ;
   bool ok_dj_DDiF = true ;
   bool ok_dj_ddiF = true ;

   size_t nb_pts = 100 ;
   for( size_t i=0 ; i!=nb_pts ; ++i )
   {
      double c1 = (double) i/(nb_pts-1) ;
      for( size_t j=0 ; j!=nb_pts ; ++j )
      {
         double c2 = (double) j/(nb_pts-1) ;
         
         if( has_dj_ddiF ) check_dj_ddiF( energy, c1, c2, ok_dj_ddiF ) ;
         
         check_F_DDiF( energy, c1, c2, ok_DDiF ) ;
         check_dj_DDiF( energy, c1, c2, ok_dj_DDiF ) ;

         double c3 = 1.0 - c1 - c2 ;
         check_diF( energy, c1, c2, c3, ok_diF ) ;
         c3 += PEL::pi() ; // why not ?    
         check_diF( energy, c1, c2, c3, ok_diF ) ;
      }
   }
   notify_one_test_result( energy->type_name() + "  diF/F", ok_diF ) ;
   notify_one_test_result( energy->type_name() + "  DDiF/F", ok_DDiF) ;
   notify_one_test_result( energy->type_name() + "  dj_DDiF/DDiF", ok_dj_DDiF) ;
   if( has_dj_ddiF ) 
      notify_one_test_result( energy->type_name() + " dj_ddiF/diF", 
                              ok_dj_ddiF) ;
}

//----------------------------------------------------------------------------
void
CH_BulkEnergy_TEST:: check_diF( CH_BulkEnergy const* energy,
                                   double c1, double c2, double c3,
                                   bool& ok ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "CH_BulkEnergy_TEST:: check_diF" ) ;

   double FF = energy->F( c1, c2, c3 ) ;
   double d1F = ( energy->F( c1+HH, c2, c3 ) - FF ) / HH ;
   double d2F = ( energy->F( c1, c2+HH, c3 ) - FF ) / HH ;
   double d3F = ( energy->F( c1, c2, c3+HH ) - FF ) / HH ;

   double d1F_theo = energy->diF( c1, c2, c3, c1, c2, c3, 1 ) ;
   double d2F_theo = energy->diF( c1, c2, c3, c1, c2, c3, 2 ) ;
   double d3F_theo = energy->diF( c1, c2, c3, c1, c2, c3, 3 ) ;

   bool eq = PEL::double_equality( d1F_theo, d1F, D_EPS, D_MIN ) ;
   if( !eq ) display_error( "d1F", d1F_theo, d1F ) ;
   ok = ok && eq ;
   eq = PEL::double_equality( d2F_theo, d2F, D_EPS, D_MIN ) ;
   if( !eq ) display_error( "d2F", d2F_theo, d2F ) ;
   ok = ok && eq ;
   eq = PEL::double_equality( d3F_theo, d3F, D_EPS, D_MIN ) ;
   if( !eq ) display_error( "d3F", d3F_theo, d3F ) ;
   ok = ok && eq ;
}
 
//----------------------------------------------------------------------------
void
CH_BulkEnergy_TEST:: check_F_DDiF( CH_BulkEnergy const* energy,
                                      double c1, double c2, bool& ok ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "CH_BulkEnergy_TEST:: check_F_DDiF" ) ;

   double c3 = 1.0 - c1 - c2 ;

   double d1F = energy->diF( c1, c2, c3, c1, c2, c3, 1 ) ; 
   double d2F = energy->diF( c1, c2, c3, c1, c2, c3, 2 ) ; 
   double d3F = energy->diF( c1, c2, c3, c1, c2, c3, 3 ) ; 

   double f1_diff = 4.0 * ST * ( (d1F-d2F)/S2 + (d1F-d3F)/S3 ) / THICKNESS ;
   double f1_theo = energy->DDiF( c1, c2, c1, c2, 1, THICKNESS ) ;

   bool eq = PEL::double_equality( f1_theo, f1_diff, D_EPS, D_MIN ) ;
   if( !eq ) display_error( "DDiF(i=1)/F", f1_theo, f1_diff ) ;
   ok = ok && eq ;

   double f2_diff = 4.0 * ST * ( (d2F-d1F)/S1 + (d2F-d3F)/S3 ) / THICKNESS ;
   double f2_theo = energy->DDiF( c1, c2, c1, c2, 2, THICKNESS ) ;

   eq = PEL::double_equality( f2_theo, f2_diff, D_EPS, D_MIN ) ;
   if( !eq ) display_error( "DDiF(i=2)/F", f2_theo, f2_diff ) ;
   ok = ok && eq ;
}

//----------------------------------------------------------------------------
void
CH_BulkEnergy_TEST:: check_dj_DDiF( CH_BulkEnergy const* energy, 
                                       double c1, double c2, bool& ok ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "CH_BulkEnergy_TEST:: check_dj_DDiF" ) ;

   double c1_exp = 2.*c1 ;
   double c2_exp = 2.*c2 ;

   for( size_t i=1 ; i!=3 ; ++i )
   {
      size_t j = 1 ;
      double d1_DiF_theo = energy->dj_DDiF( c1, c2, 
                                            c1_exp, c2_exp, i, j, THICKNESS ) ;
      double d1_DiF_diff = ( energy->DDiF( c1+HH, c2, 
                                           c1_exp, c2_exp, i, THICKNESS ) - 
                             energy->DDiF( c1, c2, 
                                           c1_exp, c2_exp, i, THICKNESS )  
                           ) / HH ;

      bool eq = PEL::double_equality( d1_DiF_theo, d1_DiF_diff, D_EPS, D_MIN ) ;
      if( !eq )
      {
         std::ostringstream mesg ;
         mesg <<  "d1_DD" << i << "F" ;
         display_error( mesg.str(), d1_DiF_theo, d1_DiF_diff) ;
      }
      ok = ok && eq ;

      j = 2 ;
      double d2_DiF_theo = energy->dj_DDiF( c1, c2, 
                                            c1_exp, c2_exp, i, j, THICKNESS ) ;
      double d2_DiF_diff = ( energy->DDiF( c1, c2+HH, 
                                           c1_exp, c2_exp, i, THICKNESS ) -
                             energy->DDiF( c1, c2, 
                                           c1_exp, c2_exp, i, THICKNESS ) 
                           ) / HH ;
      eq = PEL::double_equality( d2_DiF_theo , d2_DiF_diff , D_EPS, D_MIN ) ;
      if( !eq )
      {
         std::ostringstream mesg ;
         mesg << "d2_DD" << i << "F" ;
         display_error( mesg.str(), d2_DiF_theo, d2_DiF_diff ) ;
      }
      ok = ok && eq ;
   }
}

//----------------------------------------------------------------------------
void
CH_BulkEnergy_TEST:: check_dj_ddiF( CH_BulkEnergy const* energy,
                                       double c1, double c2, bool& ok ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "CH_BulkEnergy_TEST:: check_dj_ddiF" ) ;
   double c1_exp = 2.*c1 ;
   double c2_exp = 2.*c2 ;
   for( size_t i=1 ; i!=4 ; i++)
   {
      double c3 = 1.0 - c1 - c2 ;
      double c3_exp = 1.0 - c1_exp - c2_exp ;
      double diF= energy->diF( c1, c2, c3, c1_exp, c2_exp, c3_exp, i ) ;
      
      double c1HH = c1+HH ;
      double c3HH = 1.0 - c1HH - c2 ;
      double d1_ddiF = 
         ( energy->diF( c1HH, c2, c3HH, c1_exp, c2_exp, c3_exp, i ) - diF ) 
         / HH ;
      
      double c2HH = c2+HH ;
      c3HH = 1.0 - c1 - c2HH ;
      double d2_ddiF =
         ( energy->diF( c1, c2HH, c3HH, c1_exp, c2_exp, c3_exp, i ) - diF ) 
         / HH ;
      
      double d1_ddiF_theo = energy->dj_ddiF( c1, c2, 
                                                  c1_exp, c2_exp, i, 1 ) ;
      double d2_ddiF_theo = energy->dj_ddiF( c1, c2, 
                                                  c1_exp, c2_exp, i, 2 ) ;
      bool eq = PEL::double_equality(  d1_ddiF, d1_ddiF_theo, D_EPS, D_MIN ) ;
      if( !eq )
      {
         std::ostringstream mesg ;
         mesg <<  "d1_dd" << i << "F" ;
         display_error( mesg.str(), d1_ddiF, d1_ddiF_theo ) ;
      }
      ok = ok && eq ;
      eq = PEL::double_equality( d2_ddiF, d2_ddiF_theo, D_EPS, D_MIN ) ;
      if( !eq ) 
      {
         std::ostringstream mesg ;
         mesg <<  "d2_dd" << i << "F" ;
         display_error( mesg.str(), d2_ddiF, d2_ddiF_theo ) ;
      }
      ok = ok && eq ;
   }
}

//----------------------------------------------------------------------------
void
CH_BulkEnergy_TEST:: display_error( std::string const& mesg,
                                       double xx_1, 
                                       double xx_2 ) const
//----------------------------------------------------------------------------
{
   ios_base::fmtflags original_flags = cout.flags() ;
   cout.setf( ios_base::uppercase | ios_base::scientific ) ;

   cout << std::setw( 20 ) << mesg << " "
        << std::setprecision( 10 ) << std::setw( 20 ) << xx_1
        << std::setprecision( 10 ) << std::setw( 20 ) << xx_2
        << endl ;

   cout.flags( original_flags ) ;
}
