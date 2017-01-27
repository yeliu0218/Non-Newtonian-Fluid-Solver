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

#include <GE_SegPolyIntersector_TEST.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh> 

#include <doubleVector.hh>
#include <size_t_vector.hh>
#include <stringVector.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_Segment.hh>
#include <GE_SegmentPolyhedron_INT.hh>
#include <GE_SetOfPoints.hh>

#include <string>

//-------------------------------------------------------------------------
GE_SegPolyIntersector_TEST*
GE_SegPolyIntersector_TEST::unique_instance = new GE_SegPolyIntersector_TEST() ;
//-------------------------------------------------------------------------

//----------------------------------------------------------------------------
GE_SegPolyIntersector_TEST:: GE_SegPolyIntersector_TEST( void ) 
//----------------------------------------------------------------------------
   : PEL_ObjectTest( "GE_SegmentPolyhedron_INT", "GE_SegPolyIntersector_TEST" )
   , ITS( 0 )
   , M( 0 )
   , S( 0 )
{
}

//----------------------------------------------------------------------------
GE_SegPolyIntersector_TEST:: ~GE_SegPolyIntersector_TEST( void )
//----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
void
GE_SegPolyIntersector_TEST:: process_one_test( PEL_ModuleExplorer const* texp )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegPolyIntersector_TEST:: process_one_test" ) ;

   // Segment and polyhedron creation
   create_polyhedra( texp ) ;

   // Build intersector
   std::string int_name = "" ;
   size_t const dim = M->nb_space_dimensions() ;
   if( dim == 2 )
   {
      int_name = "GE_SegmentPolyhedron1D_INT" ;
   }
   else 
   {
      int_name = "GE_SegmentPolyhedron2D_INT" ;
   }
   PEL_ModuleExplorer* se = 0 ;
   if( texp->has_module( "GE_SegmentPolyhedron_INT" ) )
   {
      se = texp->create_subexplorer( 0, "GE_SegmentPolyhedron_INT" ) ;
      int_name = se->string_data( "concrete_name" ) ;
   }
   ITS = GE_SegmentPolyhedron_INT::make( 0, int_name, dim, se ) ;

   test_intersection( texp, int_name ) ;

   ITS->destroy() ; ITS = 0 ;
   if( se != 0 ) se->destroy() ;
   S = 0 ; M = 0 ;
}

//-----------------------------------------------------------------------------
void
GE_SegPolyIntersector_TEST:: test_intersection( PEL_ModuleExplorer const* texp,
                                                std::string const& outp )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegPolyIntersector_TEST:: test_intersection" ) ;

   // Check intersection
   ITS->check_intersection( S->vertex(0), S->vertex(1), M ) ;
   PEL_ASSERT( ITS->intersection_checked() ) ;

    // Recover solution
   PEL_ModuleExplorer* sexp = texp->create_subexplorer( 0, "solution" ) ;
   GE_Point* SOL_PT = GE_Point::create( 0, ITS->nb_space_dimensions() ) ;
   GE_Point* INT_PT = GE_Point::create( 0, ITS->nb_space_dimensions() ) ;
   bool result = false ;
   std::string const& type = sexp->string_data( "type" ) ;
   if( type == "one_single_intersection" )
   {
      result = true ;
      SOL_PT->set_coordinates( 
                       sexp->doubleVector_data( "intersection_point" ) ) ;
   }
   else if( type == "no_single_intersection" )
   {
      result = false ;
   }
   else
   {
      PEL_Error::object()->raise_bad_data_value(
         sexp, "type",
         "  - \"one_single_intersection\"\n"
         "  - \"no_single_intersection\"" ) ;
   }
   sexp->destroy() ; sexp = 0 ;

   // Test to provided solution
   bool ok = true ;
   if( ITS->one_single_intersection() )
   {
       ok = result ;
       if( ok )
       {
          ITS->intersection_point( INT_PT ) ;
          ok = ( SOL_PT->distance( INT_PT )<1.E-08 ) ;
       }
   }
    else
    {
       ok = !result ;
    }

    if( !ok )
    {
       out() << "first pt: " ;
       S->vertex(0)->print( out(), 0 ) ;
       out() << std::endl ;
       out() << "second pt: " ;
       S->vertex(1)->print( out(), 0 ) ;
       out() << std::endl ;
       out() << "polyhedron: " ;
       M->print( out(), 0 ) ;
       if( ITS->one_single_intersection() )
       {
	  out() << "inter: " ;
          ITS->intersection_point( INT_PT ) ;
          INT_PT->print( out(), 0 ) ;
          out() << std::endl ;
       }
       else
       { 
	  out() << "No intersection found" << std::endl ;
       }
    }
    notify_one_test_result( outp, ok ) ;

    SOL_PT->destroy() ; SOL_PT = 0 ;
    INT_PT->destroy() ; INT_PT = 0 ;
}

//-----------------------------------------------------------------------------
void
GE_SegPolyIntersector_TEST:: create_polyhedra( PEL_ModuleExplorer const* texp )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_SegPolyIntersector_TEST:: create_polyhedra" ) ;

   size_t const dim = texp->int_data( "nb_space_dimensions" ) ;
   texp->test_data_in( "nb_space_dimensions", "2,3" ) ;

   // Polyhedron
   PEL_ModuleExplorer* sexp = texp->create_subexplorer( 0, "polyhedron" ) ;
   doubleVector const& vec = sexp->doubleVector_data( "coordinates_of_vertices" ) ;
   doubleVector coord( dim ) ;
   size_t poly_nb_points  = vec.size()/dim ;
   GE_SetOfPoints* vm = GE_SetOfPoints:: create( this, dim ) ;
   size_t_vector vertices( poly_nb_points ) ;
   for( size_t i = 0 ; i<poly_nb_points; ++i )
   {
      for( size_t j=0; j<dim; j++ )
         coord( j ) = vec(dim*i+j) ;
      GE_Point* pt = GE_Point::create( 0, coord ) ;
      if( vm->has( pt ) )
         vertices(i) =  vm->index( pt ) ;
      else
      {
         vm->append( pt ) ;
         vertices(i) = i ;
      }
      pt->destroy() ; pt = 0 ;
   }

   M = GE_Mpolyhedron:: create( sexp->string_data( "name" ),
                                vm, vertices ) ;
   sexp->destroy(); sexp = 0 ;

   // Segment
   sexp = texp->create_subexplorer( 0, "segment" ) ;
   doubleVector const& vec2 
                   = sexp->doubleVector_data( "coordinates_of_vertices" ) ;
   PEL_ASSERT( vec2.size()==2*dim ) ;
   vertices.re_initialize( 2 ) ;
   for( size_t i=0; i<2; ++i )
   {
      for( size_t j=0; j<dim; j++ )
         coord( j ) = vec2(dim*i+j) ;
      GE_Point* pt = GE_Point::create( 0, coord ) ;
      if( vm->has( pt ) )
         vertices(i) =  vm->index( pt ) ;
      else
      {
         vm->append( pt ) ;
         vertices(i) = vm->index( pt ) ;
      }
      pt->destroy() ; pt = 0 ;
   }
   GE_Mpolyhedron* poly 
             = GE_Mpolyhedron::create( "GE_Segment", vm, vertices ) ;
   PEL_ASSERT( dynamic_cast<GE_Segment*>( poly ) ) ;
   S = static_cast<GE_Segment*>( poly ) ;
   sexp->destroy(); sexp = 0 ;
}
