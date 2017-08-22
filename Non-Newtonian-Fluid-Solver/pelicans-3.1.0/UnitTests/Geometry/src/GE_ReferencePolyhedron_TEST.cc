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

#include <GE_ReferencePolyhedron_TEST.hh>

#include <doubleArray2D.hh>
#include <stringVector.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_Vector.hh>

#include <LA_DenseMatrix.hh>
#include <LA_GaussLU_DS.hh>
#include <LA_SeqVector.hh>

#include <PEL.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>

#include <iostream>
#include <sstream>
using std::ostringstream ;
using std::endl ;

//---------------------------------------------------------------------------
GE_ReferencePolyhedron_TEST*
GE_ReferencePolyhedron_TEST:: REGISTRATOR = new GE_ReferencePolyhedron_TEST() ;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
GE_ReferencePolyhedron_TEST:: GE_ReferencePolyhedron_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "GE_ReferencePolyhedron", "GE_ReferencePolyhedron_TEST" )
{
}

//---------------------------------------------------------------------------
GE_ReferencePolyhedron_TEST:: ~GE_ReferencePolyhedron_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
GE_ReferencePolyhedron_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron_TEST:: process_one_test" ) ;
   PEL_CHECK( exp!=0 ) ;

   std::string const& tname = exp->string_data( "polyhedron" ) ;
   GE_ReferencePolyhedron const* ref =
                         GE_Mpolyhedron::reference_polyhedron( tname ) ;

   out() << "| ... " << tname << "/" << ref->name() << " : " << std::endl ;

   PEL_ModuleExplorer const* e = exp->create_subexplorer( 0, "geometry" ) ;
   process_geometry( ref, e ) ;
   e->destroy() ; e = 0 ;

   process_face_containment( ref ) ;
}

//---------------------------------------------------------------------------
void
GE_ReferencePolyhedron_TEST:: process_geometry(
                                            GE_ReferencePolyhedron const* ref,
                                            PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron_TEST:: process_geometry" ) ;
   PEL_CHECK( ref!=0 ) ;
   PEL_CHECK( exp!=0 ) ;

   // Name :
   bool const n = ref->name()==exp->string_data( "reference_polyhedron" ) ;
   notify_one_test_result( "name", n ) ;

   // Dimension :
   bool const dim = ref->dimension()==(size_t) exp->int_data( "dimension" ) ;
   notify_one_test_result( "dimension", dim  ) ;

   // Number of vertices :
   bool const nb_verts =
                 ref->nb_vertices()==(size_t) exp->int_data( "nb_vertices" ) ;
   notify_one_test_result( "nb_vertices", nb_verts ) ;

   // Data checker:
   if( exp->has_entry( "vertex_coordinates" ) )
   {
      if( ref->dimension()==0 )
      {
         PEL_Error::object()->raise_data_error(
            exp, "vertex_coordinates", "test is not expected in 0D" ) ;
      }
   }
   
   // Vertices coordinates :
   if( ref->dimension() != 0 )
   {
      doubleArray2D const& v_coords =
         exp->doubleArray2D_data( "vertex_coordinates" ) ;
      PEL_ASSERT( v_coords.index_bound(0) ==
                  (size_t) exp->int_data( "nb_vertices" ) ) ;
      PEL_ASSERT( v_coords.index_bound(1) ==
                  (size_t) exp->int_data( "dimension" ) ) ;
      bool vcoord = true ;
      for( size_t i=0 ; i<v_coords.index_bound(0) ; ++i )
      {
         GE_Point const* v = ref->vertex(i) ;
         bool ok = true ;
         for( size_t j=0 ; ok && j<v_coords.index_bound(1) ; ++j )
         {
            ok &= ( v_coords(i,j) == v->coordinate(j) ) ;
         }
         vcoord &= ok ;
         if( !ok )
         {
            out() << "Vertex: " << i << std::endl ;
            out() << "   expected:" ;
            for( size_t j=0 ; j<v_coords.index_bound(1) ; ++j )
            {
               out() << " " << v_coords(i,j) ;
            }
            out() << std::endl ;
            out() << "   computed:" ;
            for( size_t j=0 ; j<v_coords.index_bound(1) ; ++j )
            {
               out() << " " << v->coordinate(j) ;
            }
            out() << std::endl ;
         }
      }
      notify_one_test_result( "vertex", vcoord ) ;
   }

   // Number of faces :
   bool const nb_faces = ref->nb_faces()==(size_t) exp->int_data( "nb_faces" ) ;
   notify_one_test_result( "nb_faces", nb_faces ) ;


   // Data check:
   if( exp->has_entry( "face_outward_normals" ) )
   {
      if( ref->dimension()==0 )
      {
         PEL_Error::object()->raise_bad_data_value(
            exp, "face_outward_normals",
            "Not allowed in 0D" ) ;
      }
   }
   
   // Faces connectivity and normal :
   if( ref->dimension() != 0 )
   {
      doubleArray2D const& no =
         exp->doubleArray2D_data( "face_outward_normals" ) ;
      PEL_ASSERT( no.index_bound(0) ==
                  (size_t) exp->int_data( "nb_faces" ) ) ;
      PEL_ASSERT( no.index_bound(1) ==
                  (size_t) exp->int_data( "dimension" ) ) ;
      bool normal = true ;
      for( size_t i=0 ; i<no.index_bound(0) ; ++i )
      {
         GE_Vector const* nn = ref->face_outward_normal(i) ;
         bool ok = true ;
         for( size_t j=0 ; ok && j<no.index_bound(1) ; ++j )
         {
            ok &= PEL::double_equality( no(i,j), nn->component(j),
                                        1.e-12, 1.e-30 ) ;
         }
         normal &= ok ;
         if( !ok )
         {
            out() << "Face " << i << std::endl ;
            out() << "   connectivity: " ;
            for( size_t j=0 ; j<ref->nb_face_vertices(i) ; ++j )
            {
               out() << " " << ref->face_vertex(i,j) ;
            }
            out() << std::endl ;
            out() << "   normal: "<< std::endl ;
            out() << "       expected:" ;
            for( size_t j=0 ; j<no.index_bound(1) ; ++j )
            {
               out() << " " << no(i,j) ;
            }
            out() << std::endl ;
            out() << "       computed:" ;
            for( size_t j=0 ; j<no.index_bound(1) ; ++j )
            {
               out() << " " << nn->component(j) ;
            }
            out() << std::endl ;
         }
      }
      notify_one_test_result( "face_outward_normal", normal ) ;
   }
}

//--------------------------------------------------------------------------
void
GE_ReferencePolyhedron_TEST:: process_face_containment(
                                         GE_ReferencePolyhedron const* ref )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron_TEST:: process_face_containment" ) ;

   for( size_t i_face=0 ; i_face<ref->nb_faces() ; ++i_face )
   {
      bool ok_face = true ;
      if( ref->nb_face_vertices( i_face ) == 1 )
      {
         GE_Point const* V0 = ref->vertex( ref->face_vertex( i_face, 0 ) ) ;

         check_face_containment_1( ref, V0, i_face, ok_face ) ;
      }
      if( ref->nb_face_vertices( i_face ) == 2 )
      {
         GE_Point const* V0 = ref->vertex( ref->face_vertex( i_face, 0 ) ) ;
         GE_Point const* V1 = ref->vertex( ref->face_vertex( i_face, 1 ) ) ;

         check_face_containment_1( ref, V0, i_face, ok_face ) ;
         check_face_containment_1( ref, V1, i_face, ok_face ) ;

         check_face_containment_2( ref, V0, V1, i_face, ok_face ) ;
      }
      else if( ref->nb_face_vertices( i_face ) == 3 )
      {
         GE_Point const* V0 = ref->vertex( ref->face_vertex( i_face, 0 ) ) ;
         GE_Point const* V1 = ref->vertex( ref->face_vertex( i_face, 1 ) ) ;
         GE_Point const* V2 = ref->vertex( ref->face_vertex( i_face, 2 ) ) ;

         check_face_containment_1( ref, V0, i_face, ok_face ) ;
         check_face_containment_1( ref, V1, i_face, ok_face ) ;
         check_face_containment_1( ref, V2, i_face, ok_face ) ;

         check_face_containment_2( ref, V0, V1, i_face, ok_face ) ;
         check_face_containment_2( ref, V1, V2, i_face, ok_face ) ;
         check_face_containment_2( ref, V2, V0, i_face, ok_face ) ;

         check_face_containment_3( ref, V0, V1, V2, i_face, ok_face ) ;
      }
      else if( ref->nb_face_vertices( i_face ) == 4 )
      {
         GE_Point const* V0 = ref->vertex( ref->face_vertex( i_face, 0 ) ) ;
         GE_Point const* V1 = ref->vertex( ref->face_vertex( i_face, 1 ) ) ;
         GE_Point const* V2 = ref->vertex( ref->face_vertex( i_face, 2 ) ) ;
         GE_Point const* V3 = ref->vertex( ref->face_vertex( i_face, 3 ) ) ;

         check_face_containment_1( ref, V0, i_face, ok_face ) ;
         check_face_containment_1( ref, V1, i_face, ok_face ) ;
         check_face_containment_1( ref, V2, i_face, ok_face ) ;
         check_face_containment_1( ref, V3, i_face, ok_face ) ;

         check_face_containment_2( ref, V0, V1, i_face, ok_face ) ;
         check_face_containment_2( ref, V1, V2, i_face, ok_face ) ;
         check_face_containment_2( ref, V2, V3, i_face, ok_face ) ;
         check_face_containment_2( ref, V3, V0, i_face, ok_face ) ;

         check_face_containment_4( ref, V0, V1, V2, V3, i_face, ok_face ) ;
      }
      ostringstream mesg ;
      mesg << "containment in face " << i_face ;
      notify_one_test_result( mesg.str(), ok_face ) ;
   }
}

//---------------------------------------------------------------------------
void
GE_ReferencePolyhedron_TEST:: check_face_containment_1(
                                    GE_ReferencePolyhedron const* ref,
                                    GE_Point const* V0,
                                    size_t i_face,
                                    bool& ok_face )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron_TEST:: check_face_containment_2" ) ;

   bool fc_class = ref->face_contains( i_face, V0 ) ;
   if( !fc_class ) display_face_containment_error( ref, i_face, V0, true ) ;
   ok_face = ok_face && fc_class ;
}

//---------------------------------------------------------------------------
void
GE_ReferencePolyhedron_TEST:: check_face_containment_2(
                                    GE_ReferencePolyhedron const* ref,
                                    GE_Point const* V0,
                                    GE_Point const* V1,
                                    size_t i_face,
                                    bool& ok_face )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron_TEST:: check_face_containment_2" ) ;

   double tol = ref->epsilon() ;

   size_t nb_pts = 100 ;
   double a_min = -2.0 ;
   double a_max =  2.0 ;
   GE_Point* pt = GE_Point::create( 0, ref->dimension() ) ;

   for( size_t i0=0 ; i0<nb_pts ; ++i0 )
   {
      double a0 = a_min + i0*(a_max-a_min)/(nb_pts-1) ;
      for( size_t i1=0 ; i1<nb_pts ; ++i1 )
      {
         double a1 = a_min + i1*(a_max-a_min)/(nb_pts-1) ;
         double sum = a0+a1 ;
         if( PEL::abs( sum ) < 1.e-8 ) continue ;
         a0 /= sum ;
         a1 /= sum ;
         for( size_t ic=0 ; ic<ref->dimension() ; ++ic )
         {
            double xx = a0 * V0->coordinate( ic ) +
                        a1 * V1->coordinate( ic ) ;
            pt->set_coordinate( ic, xx ) ;
         }
         bool fc_class = ref->face_contains( i_face, pt ) ;
         bool fc_theo ;
         if( a0+tol>0. && a0-tol<1. && a1+tol>0. && a1-tol<1. )
         {
            fc_theo = true ;
         }
         else
         {
            fc_theo = false ;
         }
         ok_face = ok_face && ( fc_theo == fc_class ) ;
         if( fc_theo != fc_class )
         {
            out() << "a0=" << a0 << endl ;
            out() << "a1=" << a1 << endl ;
            display_face_containment_error( ref, i_face, pt, fc_class ) ;
         }
      }
   }
   pt->destroy() ;
}

//---------------------------------------------------------------------------
void
GE_ReferencePolyhedron_TEST:: check_face_containment_3(
                                    GE_ReferencePolyhedron const* ref,
                                    GE_Point const* V0,
                                    GE_Point const* V1,
                                    GE_Point const* V2,
                                    size_t i_face,
                                    bool& ok_face )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron_TEST:: check_face_containment_3" ) ;

   double tol = ref->epsilon() ;

   size_t nb_pts = 50 ;
   double a_min = -2.0 ;
   double a_max =  2.0 ;
   GE_Point* pt = GE_Point::create( 0, ref->dimension() ) ;

   for( size_t i0=0 ; i0<nb_pts ; ++i0 )
   {
      double a0 = a_min + i0*(a_max-a_min)/(nb_pts-1) ;
      for( size_t i1=0 ; i1<nb_pts ; ++i1 )
      {
         double a1 = a_min + i1*(a_max-a_min)/(nb_pts-1) ;
         for( size_t i2=0 ; i2<nb_pts ; ++i2 )
         {
            double a2 = a_min + i2*(a_max-a_min)/(nb_pts-1) ;
            double sum = a0+a1+a2 ;
            if( PEL::abs( sum ) < 1.e-8 ) continue ;
            a0 /= sum ;
            a1 /= sum ;
            a2 /= sum ;
            for( size_t ic=0 ; ic<ref->dimension() ; ++ic )
            {
               double xx = a0 * V0->coordinate( ic ) +
                           a1 * V1->coordinate( ic ) +
                           a2 * V2->coordinate( ic ) ;
               pt->set_coordinate( ic, xx ) ;
            }
            bool fc_class = ref->face_contains( i_face, pt ) ;
            bool fc_theo ;
            if( a0+tol>0. && a0-tol<1. &&
                a1+tol>0. && a1-tol<1. &&
                a2+tol>0. && a2-tol<1. )
            {
               fc_theo = true ;
            }
            else
            {
               fc_theo = false ;
            }
            ok_face = ok_face && ( fc_theo == fc_class ) ;
            if( fc_theo != fc_class )
            {
               out() << "a0=" << a0 << endl ;
               out() << "a1=" << a1 << endl ;
               out() << "a2=" << a2 << endl ;
               display_face_containment_error( ref, i_face, pt, fc_class ) ;
            }
         }
      }
   }
   pt->destroy() ;
}

//---------------------------------------------------------------------------
void
GE_ReferencePolyhedron_TEST:: check_face_containment_4(
                                    GE_ReferencePolyhedron const* ref,
                                    GE_Point const* V0,
                                    GE_Point const* V1,
                                    GE_Point const* V2,
                                    GE_Point const* V3,
                                    size_t i_face,
                                    bool& ok_face )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_ReferencePolyhedron_TEST:: check_face_containment_3" ) ;

   size_t nb_pts = 20 ;
   double a_min = -2.0 ;
   double a_max =  2.0 ;
   GE_Point* pt = GE_Point::create( 0, ref->dimension() ) ;

   for( size_t i0=0 ; i0<nb_pts ; ++i0 )
   {
      double a0 = a_min + i0*(a_max-a_min)/(nb_pts-1) ;
      for( size_t i1=0 ; i1<nb_pts ; ++i1 )
      {
         double a1 = a_min + i1*(a_max-a_min)/(nb_pts-1) ;
         for( size_t i2=0 ; i2<nb_pts ; ++i2 )
         {
            double a2 = a_min + i2*(a_max-a_min)/(nb_pts-1) ;
            for( size_t i3=0 ; i3<nb_pts ; ++i3 )
            {
               double a3 = a_min + i3*(a_max-a_min)/(nb_pts-1) ;
               double sum = a0+a1+a2+a3 ;
               if( PEL::abs( sum ) < 1.e-8 ) continue ;
               a0 /= sum ;
               a1 /= sum ;
               a2 /= sum ;
               a3 /= sum ;
               for( size_t ic=0 ; ic<ref->dimension() ; ++ic )
               {
                  double xx = a0 * V0->coordinate( ic ) +
                              a1 * V1->coordinate( ic ) +
                              a2 * V2->coordinate( ic ) +
                              a3 * V3->coordinate( ic ) ;
                  pt->set_coordinate( ic, xx ) ;
               }
               bool fc_class = ref->face_contains( i_face, pt ) ;
               bool ok = true ;
               if( a0>0. && a0<1. &&
                   a1>0. && a1<1. &&
                   a2>0. && a2<1. &&
                   a3>0. && a3<1. )
               {
                  ok = ok && fc_class;
               }
               if( a0+a1+a2+a3 < 0. )
               {
                  ok = ok && !fc_class ;
               }
               ok_face = ok_face && ok ;
               if( !ok )
               {
                  out() << "a0=" << a0 << endl ;
                  out() << "a1=" << a1 << endl ;
                  out() << "a2=" << a2 << endl ;
                  out() << "a3=" << a3 << endl ;
                  display_face_containment_error( ref, i_face, pt, fc_class ) ;
               }
            }
         }
      }
   }
   pt->destroy() ;
}

//---------------------------------------------------------------------------
void
GE_ReferencePolyhedron_TEST:: display_face_containment_error(
                                      GE_ReferencePolyhedron const* ref,
                                      size_t i_face,
                                      GE_Point const* pt,
                                      bool fc_class )
//---------------------------------------------------------------------------
{
   out() << "----------------------" << endl ;
   ref->print( out() , 0 ) ;
   out() << endl << "point : " ; pt->print( out(), 3 ) ; out() << endl ;
   if( fc_class )
      out() << "is given contained in face : " << i_face ;
   else
      out() << "is not given contained in face : " << i_face ;
   out() << endl << "----------------------" << endl ;
   PEL_Error::object()->raise_plain( "exit" ) ;
}
