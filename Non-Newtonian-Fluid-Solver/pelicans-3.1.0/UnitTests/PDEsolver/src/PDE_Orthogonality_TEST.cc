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

#include <PDE_Orthogonality_TEST.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_Vector.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEcell.hh>

#include <iostream>

PDE_Orthogonality_TEST*
PDE_Orthogonality_TEST:: REGISTRATOR = new PDE_Orthogonality_TEST() ;

//---------------------------------------------------------------------------
PDE_Orthogonality_TEST:: PDE_Orthogonality_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "PDE_CursorFEside", 
                     "PDE_Orthogonality_TEST" )
{
}

//---------------------------------------------------------------------------
PDE_Orthogonality_TEST:: ~PDE_Orthogonality_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_Orthogonality_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_Orthogonality_TEST:: process_one_test" ) ;

   double my_eps = exp->double_data( "dbl_epsilon" ) ;
   double my_min = exp->double_data( "dbl_minimum" ) ;

   PEL_ModuleExplorer* se = 
                       exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   PDE_DomainAndFields* dom = PDE_DomainAndFields::create( 0, se ) ;
   se->destroy() ; se = 0 ;

   GE_Vector* cc = GE_Vector::create( 0, dom->nb_space_dimensions() ) ;
   GE_Vector* vv = GE_Vector::create( 0, dom->nb_space_dimensions() ) ;

   bool eq ;
   bool ok_exist = true ;
   bool ok_dist = true ;
   bool ok_normal = true ;

   PDE_CursorFEside* sFE = dom->create_CursorFEside( 0 ) ;

   PDE_LocalFEcell const* cFE0 = sFE->adjacent_localFEcell( 0 ) ;
   PDE_LocalFEcell const* cFE1 = sFE->adjacent_localFEcell( 1 ) ;
   
   for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
   {
      GE_Mpolyhedron const* s_poly = sFE->polyhedron() ;

      GE_Point const* c0 = cFE0->polyhedron()->finite_volume_center() ;
      GE_Point const* c1 = cFE1->polyhedron()->finite_volume_center() ;
      ok_exist = ( c0!=0 && c1!=0 ) ;

      if( ok_exist )
      {
         double h1 = c0->distance( c1 ) ;
         double h2 = 
            PEL::abs( sFE->distance_to_adjacent_finite_volume_center( 0 ) ) +
            PEL::abs( sFE->distance_to_adjacent_finite_volume_center( 1 ) ) ;
         eq = PEL::double_equality( h1, h2, my_eps, my_min ) ;
         ok_dist = ok_dist && eq ;
         if( !eq ) display_error( "distance", h1, h2 ) ;

         cc->re_initialize( c0, c1 ) ;
         for( size_t iv=1 ; iv<s_poly->nb_vertices() ; ++iv )
         {
            vv->re_initialize( s_poly->vertex( 0 ), s_poly->vertex( 1 ) ) ;
            double xx = cc->dot_product( vv ) ;
            eq = PEL::double_equality( xx, 0.0, my_eps, my_min ) ;
            ok_normal = ok_normal && eq ;
            if( !eq ) display_error( "orthogonality", xx, 0.0 ) ;
         }
      }
      else
      {
         break ;
      }
   }

   cc->destroy() ;
   vv->destroy() ;

   notify_one_test_result( "existence of FV centers", ok_exist ) ;
   if( ok_exist )
   {
      notify_one_test_result( "distances", ok_dist ) ;
      notify_one_test_result( "orthogonality", ok_normal ) ;
   }

   sFE->destroy() ;
   dom->destroy() ;
}

//-----------------------------------------------------------------------
void 
PDE_Orthogonality_TEST:: display_error( std::string const& test,
                                        double xx1, double xx2 ) const
//-----------------------------------------------------------------------
{
   std::cout << test << "  xx1 " << xx1 << " xx2 " << xx2 << std::endl ;
}
