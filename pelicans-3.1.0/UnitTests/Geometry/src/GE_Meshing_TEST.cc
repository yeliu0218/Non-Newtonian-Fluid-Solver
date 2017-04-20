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

#include <GE_Meshing_TEST.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_IndexSet.hh>
#include <PEL_Int.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>
#include <PEL_VectorIterator.hh>

#include <doubleVector.hh>
#include <size_t_vector.hh>

#include <GE_Color.hh>
#include <GE_Meshing.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_SetOfPoints.hh>

#include <PDE_CellFE.hh>
#include <PDE_DomainBuilder.hh>
#include <PDE_FaceFE.hh>
#include <PDE_GridFE.hh>

#include <fstream>
#include <ios>
#include <iomanip>
#include <iostream>

using std::cout ;
using std::endl ;
using std::ios_base ;
using std::string ;

//-------------------------------------------------------------------------
GE_Meshing_TEST*
GE_Meshing_TEST::unique_instance = new GE_Meshing_TEST() ;
//-------------------------------------------------------------------------

//----------------------------------------------------------------------------
GE_Meshing_TEST:: GE_Meshing_TEST( void ) 
//----------------------------------------------------------------------------
   : PEL_ObjectTest( "GE_Meshing", "GE_Meshing_TEST" )
{
}

//----------------------------------------------------------------------------
GE_Meshing_TEST:: ~GE_Meshing_TEST( void )
//----------------------------------------------------------------------------
{
}

//-----------------------------------------------------------------------------
void
GE_Meshing_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing_TEST:: process_one_test" ) ;

   if( exp->has_module( "traces" ) )
   {
      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "traces" ) ;
      OFS.open( se->string_data( "output_file" ).c_str() ) ;
      PEL_ASSERT( OFS ) ;
      se->destroy() ;
   }

   check_face_numbering( exp ) ;

   int space_dim = exp->int_data( "nb_space_dimensions" ) ;
   PEL_ModuleExplorer const* t =
                          exp->create_subexplorer( 0, "GE_Meshing" ) ;
   GE_Meshing* tested_meshing = GE_Meshing::create( 0, t, space_dim ) ;
   t->destroy() ; t = 0 ;

   double d_eps = 1.E-8 ;
   if( exp->has_entry( "dbl_epsilon" ) )
   {
      d_eps = exp->double_data( "dbl_epsilon" ) ;
      exp->test_data( "dbl_epsilon", "dbl_epsilon>0." ) ;
   }
   double d_min = 1.e-12 ;
   if( exp->has_entry( "dbl_minimum" ) )
   {
      d_min = exp->double_data( "dbl_minimum" ) ;
      exp->test_data( "dbl_minimum", "dbl_minimum>0." ) ;
   }
   
   
   if( exp->has_module( "GE_Meshing#VERIFY" ) )
   {
      PEL_ModuleExplorer const* v =
                         exp->create_subexplorer( 0, "GE_Meshing#VERIFY" ) ;
      GE_Meshing* verified_meshing = GE_Meshing::create( 0, v, space_dim  ) ;
      v->destroy() ; v = 0 ;
   
      notify_one_test_result(
         "nb_space_dimensions",
         tested_meshing->nb_space_dimensions()
                                 ==verified_meshing->nb_space_dimensions() ) ;
      
      bool ok_total = true ;
      check_vertices( tested_meshing, verified_meshing, d_eps, d_min, ok_total ) ;
      check_faces( tested_meshing, verified_meshing, ok_total ) ;
      check_cells( tested_meshing, verified_meshing, ok_total ) ;
      if( !ok_total )
      {
         std::string nn = exp->name() + "_TESTED" ;
         std::ofstream ff0( nn.c_str() ) ;
         // std::cout << "Tested mesh : " << std::endl ;
         tested_meshing->print_as_an_explicit_meshing( ff0 ) ;
         ff0.close() ;
         
         nn = exp->name() + "_VERIFY" ;
         std::ofstream ff1( nn.c_str() ) ;
         verified_meshing->print_as_an_explicit_meshing( ff1 ) ;
         ff1.close() ;
      }         

      verified_meshing->destroy() ; verified_meshing = 0 ;
   }

   tested_meshing->destroy() ; tested_meshing = 0 ;
   if( OFS ) OFS.close() ;
}

//-------------------------------------------------------------------------
void
GE_Meshing_TEST:: check_vertices( GE_Meshing* testM, GE_Meshing* verifM,
                                  double d_eps, double d_min,
                                  bool& ok_total )
//-------------------------------------------------------------------------
{
   bool ok = ( testM->nb_vertices() == verifM->nb_vertices() ) ;
   notify_one_test_result( "nb_vertices", ok ) ;
   
   if( ok )
   {
      bool ok_vert = true ;
      bool ok_col = true ;
      // To avoid vertex coordinates depend on order used to loop on vertex
      for( size_t i=testM->nb_vertices()/2 ; i>0  ; i-- )
      {
         testM->start_vertex_iterator() ;
         for( size_t j=0 ; j<i ; j++ )
         {
            testM->go_next_vertex() ;
         }
         testM->vertex_coordinates() ;
      }
      
      testM->start_vertex_iterator() ;
      verifM->start_vertex_iterator() ;
      for( ; testM->valid_vertex() ;
           testM->go_next_vertex(), verifM->go_next_vertex() )
      {
         PEL_ASSERT( verifM->valid_vertex() ) ;
  
         doubleVector vctest = testM->vertex_coordinates() ;
         doubleVector vcverif = verifM->vertex_coordinates() ;
         ok_vert &= ( vctest.size() == vcverif.size() ) ;
         
         for( size_t i=0 ; ok_vert && i<vctest.size() ; i++ )
         {
            ok_vert &= PEL::double_equality(
                                  vctest(i), vcverif(i), d_eps, d_min ) ;
         }
         ok_col &= verifM->vertex_color()->is_equal( testM->vertex_color() ) ;
      }
      notify_one_test_result( "vertex_coordinates", ok_vert ) ;
      notify_one_test_result( "vertex_color", ok_col ) ;
      ok = ( ok_vert && ok_col ) ;
   }
   
   ok_total &= ok ;
}

//-------------------------------------------------------------------------
void
GE_Meshing_TEST:: check_faces( GE_Meshing* testM, GE_Meshing* verifM,
                               bool& ok_total )
//-------------------------------------------------------------------------
{
   bool ok = testM->nb_faces() == verifM->nb_faces() ;
   notify_one_test_result( "nb_faces", ok ) ;
   
   if( ok )
   {
      bool ok_vert = true ;
      bool ok_poly = true ;
      bool ok_color = true ;
      verifM->start_face_iterator() ;
      testM->start_face_iterator() ;
      for( ; testM->valid_face() ;
           testM->go_next_face(), verifM->go_next_face() )
      {
         PEL_ASSERT( verifM->valid_face() ) ;
         ok_vert = ok_vert && 
              testM->face_vertices() == verifM->face_vertices() ;
         ok_poly = ok_poly &&
              testM->face_polyhedron_name()==verifM->face_polyhedron_name() ;
         ok_color = ok_color &&
              testM->face_color()==verifM->face_color() ;
            
      }
      notify_one_test_result( "face_vertices", ok_vert ) ;
      notify_one_test_result( "face_polyhedron_name", ok_poly ) ;
      notify_one_test_result( "face_color", ok_color ) ;
      ok = ( ok_vert && ok_poly && ok_color ) ;   
   }
   
   ok_total &= ok ;
}

//-------------------------------------------------------------------------
void
GE_Meshing_TEST:: check_cells( GE_Meshing* testM, GE_Meshing* verifM,
                               bool& ok_total )
//-------------------------------------------------------------------------
{
   bool ok = testM->nb_cells() == verifM->nb_cells() ;
   notify_one_test_result( "nb_cells", ok ) ;
   
   if( ok )
   {
      bool ok_vert = true ;
      bool ok_poly = true ;
      bool ok_color = true ;
      bool ok_conn = true ;
      verifM->start_cell_iterator() ;
      testM->start_cell_iterator() ;
      for( ; testM->valid_cell() ;
           testM->go_next_cell(), verifM->go_next_cell() )
      {
         PEL_ASSERT( verifM->valid_cell() ) ;
         ok_vert = ok_vert &&
            testM->cell_vertices() == verifM->cell_vertices() ;
     	 ok_poly = ok_poly &&
	    testM->cell_polyhedron_name()==verifM->cell_polyhedron_name() ;
	 ok_color = ok_color &&
	    testM->cell_color()==verifM->cell_color() ;       
         ok_conn = ok_conn &&
	    testM->cell_faces()==verifM->cell_faces() ;       
      }
      notify_one_test_result( "cell_vertices", ok_vert ) ;
      notify_one_test_result( "cell_polyhedron_name", ok_poly ) ;
      notify_one_test_result( "cell_color", ok_color ) ;
      notify_one_test_result( "cell_faces", ok_conn ) ;
      ok = ( ok_vert && ok_poly && ok_color && ok_conn ) ;
   }
   
   ok_total &= ok ;
}

//-------------------------------------------------------------------------
void
GE_Meshing_TEST:: check_face_numbering( PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Meshing_TEST:: check_face_numbering" ) ;

   PEL_Module* m = PEL_Module::create( 0, "PDE_DomainAndFields" ) ;
   PEL_ModuleExplorer const* e = exp->create_subexplorer( 0, "GE_Meshing" ) ;
   m->add_module( e->create_clone_of_attached_module( m ) ) ;
   m->add_entry( "verbose_level", PEL_Int::create( m, 0 ) ) ;
   m->add_entry( "nb_space_dimensions", 
                 PEL_Int::create( m, exp->int_data( "nb_space_dimensions" ))) ;
   m->add_module( PEL_Module::create( m, "interior_fields" ) ) ;
   e->destroy() ; e=0 ;

   PEL_ModuleExplorer const* ee = PEL_ModuleExplorer::create( 0, m ) ;
   PDE_DomainBuilder* domain = PDE_DomainBuilder::create( 0, ee, "me" ) ;
   ee->destroy() ;

   GE_SetOfPoints const* verts = domain->set_of_vertices() ;
   PDE_GridFE const* grid = domain->finite_element_grid() ;
   PEL_VectorIterator* it = PEL_VectorIterator::create( domain,
                                                        grid->cells() ) ;
   bool ok = true ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PDE_CellFE const* cell = dynamic_cast<PDE_CellFE*>( it->item() ) ;
      if( OFS ) OFS << "cell : " << cell->id_number() << endl ;
      PEL_ASSERT( cell != 0 ) ;
      GE_Mpolyhedron const* cpoly = cell->polyhedron() ;
      GE_ReferencePolyhedron const* cpoly_ref = cpoly->reference_polyhedron() ;
      PEL_VectorIterator* itf = PEL_VectorIterator::create( domain,
                                                            cell->faces() ) ;
      size_t i_face = 0 ;
      for( itf->start() ; itf->is_valid() ; itf->go_next() )
      {
         PDE_FaceFE const* face = dynamic_cast<PDE_FaceFE*>( itf->item() ) ;
         GE_Mpolyhedron const* fpoly = face->polyhedron() ;
         size_t_vector idxset_expected( fpoly->nb_vertices() ) ;
         size_t_vector idxset_observed( fpoly->nb_vertices() ) ;
         for( size_t iv=0 ; iv<fpoly->nb_vertices() ; ++iv )
         {
            idxset_expected( iv ) = cpoly_ref->face_vertex( i_face, iv ) ;
            GE_Point const* pt_1 = fpoly->vertex( iv ) ;
            size_t idx_1 = verts->index( pt_1 ) ;
            bool found = false ;
            size_t ii = 0 ;
            for( ; ii<cpoly->nb_vertices() ; ++ii )
            {
               GE_Point const* pt = cpoly->vertex( ii ) ;
               size_t idx = verts->index( pt ) ;
               if( idx == idx_1 )
               {
                  found = true ;
                  break ;
               }
            }
            idxset_observed( iv ) = ii ; 
         }
         PEL_IndexSet* s1 = PEL_IndexSet::create( 0, idxset_expected, 0 ) ;
         PEL_IndexSet* s2 = PEL_IndexSet::create( 0, idxset_observed, 1 ) ;
         bool eq = s1->is_equal( s2 ) ;
         ok = ok && eq ;
         if( OFS )
         {
            OFS << "   face : " << i_face ;
            if( eq ) 
               OFS << "  ok" << endl ;
            else 
               OFS << "   ERROR" << endl ;
            s1->print( OFS, 6 ) ; 
            OFS << "(" << idxset_expected << ")  (expected)" << endl ;
            s2->print( OFS, 6 ) ; 
            OFS << "(" << idxset_observed << ")  (observed)" << endl ;
         }
         s1->destroy() ;
         s2->destroy() ;
         i_face++ ;
      }
   }
   notify_one_test_result( exp->name() + " : face numbers", ok ) ;

   domain->destroy() ;
   m->destroy() ;
}
