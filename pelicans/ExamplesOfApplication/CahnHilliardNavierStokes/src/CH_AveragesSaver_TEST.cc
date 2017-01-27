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

#include <CH_AveragesSaver_TEST.hh>

#include <CH_AveragesSaver.hh>

#include <PEL.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>

#include <PDE_DomainAndFields.hh>

#include <FE.hh>
#include <FE_SetOfParameters.hh>

#include <doubleVector.hh>
#include <doubleArray2D.hh>
#include <size_t_vector.hh>

#include <ios>
#include <iomanip>
#include <iostream>
#include <sstream>

using std::cout ; using std::endl ;
using std::ios_base ;
using std::string ; using std::ostringstream ;

CH_AveragesSaver_TEST const* 
CH_AveragesSaver_TEST::PROTOTYPE = new CH_AveragesSaver_TEST() ;

//---------------------------------------------------------------------------
CH_AveragesSaver_TEST:: CH_AveragesSaver_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "CH_AveragesSaver", "CH_AveragesSaver_TEST" )
   , D_EPS( PEL::bad_double() )
   , D_MIN( PEL::bad_double() ) 
   , EXVOL( PEL::bad_double() ) 
   , EXCENTER( doubleVector( 0 ) ) 
   , EXVELOCITY( doubleVector( 0 ) ) 
   , EXPER( PEL::bad_double() )
{
}

//---------------------------------------------------------------------------
CH_AveragesSaver_TEST:: ~CH_AveragesSaver_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void 
CH_AveragesSaver_TEST:: process_one_test( PEL_ModuleExplorer const* exp ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_AveragesSaver_TEST:: process_one_test" ) ;
   
   D_EPS = exp->double_data( "dbl_epsilon" ) ;
   D_MIN = exp->double_data( "dbl_minimum" ) ;
   
   PEL_ModuleExplorer* se = 0 ;
   
   se = exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   PDE_DomainAndFields* dom =
             PDE_DomainAndFields::create( 0, se, PEL_Exec::communicator() ) ;
   se->destroy() ; se = 0 ;

   se = exp->create_subexplorer( 0, "FE_SetOfParameters" ) ;
   FE_SetOfParameters* prms = FE_SetOfParameters::create( dom, dom, se ) ;
   se->destroy() ; se = 0 ;

   se = exp->create_subexplorer( 0, "CH_AveragesSaver" ) ;
   CH_AveragesSaver* sv = 
      static_cast<CH_AveragesSaver*>( CH_AveragesSaver::make( dom, dom, 
                                                              prms, se ) ) ;
   se->destroy() ; se = 0 ;

   se = exp->create_subexplorer( 0, "exact_averages" ) ;
   EXVOL = se->double_data( "exact_volume" ) ;
   EXCENTER = se->doubleVector_data( "exact_center_coordinates" ) ;
   EXVELOCITY = se->doubleVector_data( "exact_center_velocity" ) ;
   EXPER = se->double_data( "exact_perimeter" ) ;
   se->destroy() ; se = 0 ;

   double per = sv->perimeter() ;
   doubleVector volume( 0 ) ;
   doubleArray2D center( 0, 0 ) ;
   doubleArray2D velocity( 0, 0 ) ;
   size_t_vector nb_cells( 0 ) ;
   sv->center_volume( volume, center, velocity, nb_cells ) ;
   
   bool ok = true ;
   for( size_t i = 0 ; i < volume.size() ; ++i )
   {
      bool eq = PEL::double_equality( EXVOL, volume( i ), D_EPS, D_MIN ) ;
      std::ostringstream mesg ;
      mesg <<  "Volume(" << i << ")" ;
      if( !eq ) display_error( mesg.str(), EXVOL, 
                                         volume( i ) ) ;
      ok &= eq ;
   }
   notify_one_test_result( "Volume", ok ) ;
   
   size_t nb_dims = dom->nb_space_dimensions() ;
   PEL_ASSERT( nb_dims == center.index_bound( 1 ) ) ;
   bool ok_center = true ;
   bool ok_velocity = true ;
   for( size_t i = 0 ; i < center.index_bound( 0 ) ; ++i )
   {
      for( size_t d = 0 ; d < nb_dims ; ++d ) 
      {
         bool eq = PEL::double_equality( EXCENTER( d ), center( i, d ), 
                                         D_EPS, D_MIN ) ;
         std::ostringstream mesg ;
         mesg <<  "Center coordinates(" << i <<"," << d << ")" ;
         if( !eq ) display_error( mesg.str(), EXCENTER( d ), center( i, d ) ) ;
         ok_center &= eq ;
         
         eq = PEL::double_equality( EXVELOCITY( d ), velocity( i, d ), 
                                    D_EPS, D_MIN ) ;
         std::ostringstream mesg2 ;
         mesg2 <<  "Velocity(" << i <<"," << d << ")" ;
         if( !eq ) display_error( mesg2.str(), EXVELOCITY( d ), 
                                               velocity( i, d ) ) ;
         ok_velocity &= eq ;
      }
   }
   notify_one_test_result( "Center coordinates", ok_center ) ;
   notify_one_test_result( "Center velocity", ok_velocity ) ;

   ok = PEL::double_equality( EXPER, per, D_EPS, D_MIN ) ;
   if( !ok ) display_error( "Perimeter", EXPER , per ) ;
   notify_one_test_result( "Perimeter", ok) ;

   dom->destroy() ;
}

//----------------------------------------------------------------------------
void
CH_AveragesSaver_TEST:: display_error( std::string const& mesg,
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


