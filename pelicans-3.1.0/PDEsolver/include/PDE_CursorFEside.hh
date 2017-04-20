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

#ifndef PDE_CURSOR_FE_SIDE_HH
#define PDE_CURSOR_FE_SIDE_HH

#include <PEL_Object.hh>

#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFE.hh>
#include <doubleArray2D.hh>
#include <size_t_vector.hh>

class GE_Color ;
class GE_Mpolyhedron ;
class GE_Point ;
class GE_QRprovider ;
class GE_CustomizedQR_provider ;
class GE_Transform ;
class GE_Vector ;

class PEL_Vector ;
class PEL_VectorIterator ;

class PDE_DiscreteField ;
class PDE_LocalFEcell ;
class PDE_GridFE ;
class PDE_FaceFE ;
class PDE_HelperFIC ;

/*
Iterators to visit the internal sides of a finite element meshing
(ie the sides that are not part of the boundary).

INSTANCE DELIVERY : `PDE_DomainsAndFields::create_CursorFEside'
*/

class PEL_EXPORT PDE_CursorFEside : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Update the specific attributes of `self' in case of modification
      // of the underlying grid.
      virtual void update( void ) ;
      
   //-- Geometry(90)

      // number of space dimensions
      size_t nb_space_dimensions( void ) const ;

   //-- Adjacent cells(92)
      
      // server whose current position is the `i_adj'-th adjacent cell
      PDE_LocalFEcell const* adjacent_localFEcell( size_t i_adj ) const ;

   //-- Handled fields(95)

      // Notify that the `order'-th spatial derivatives of  
      // `ff' basis functions will be requested from
      // `::adjacent_localFEcell'(0) and `::adjacent_localFEcell'(1).
      void require_field_calculation( PDE_DiscreteField const* ff,
                                      int order ) ;
      
      // Will the `order'-th spatial derivatives of `ff' basis function be
      // available on request from `::adjacent_localFEcell'(0) and 
      // `::adjacent_localFEcell'(1) ?
      bool field_calculation_is_handled( PDE_DiscreteField const* ff,
                                         int order ) const ;

      // Is `ff' handled by `::adjacent_localFEcell'(0) and 
      // `::adjacent_localFEcell'(1) ? 
      bool field_is_handled( PDE_DiscreteField const* ff ) const ;

      // number of fields handled by `::adjacent_localFEcell'(0) and 
      // `::adjacent_localFEcell'(1)
      size_t nb_handled_fields( void ) const ;

      // `iF'-th field handled by `::adjacent_localFEcell'(0) and 
      // `::adjacent_localFEcell'(1)
      PDE_DiscreteField const* handled_field( size_t iF ) const ;
      
   //-- Side-iterator movement(100)

      // number of sides that the side-iterator can traverse
      size_t nb_meshes( void ) const ;

      // Exclude sides for which color is `a_color'.
      void exclude_color( GE_Color const* a_color ) ;

      // Include excluded sides for which color is `a_color'.
      void include_color( GE_Color const* a_color ) ;

      // Is `a_color' excluded ?
      bool is_excluded( GE_Color const* a_color ) const ;
      
      // Move side-iterator to the `an_id_mesh'-th position and
      // configure `::adjacent_localFEcell'(0) and `::adjacent_localFEcell'(1)
      // accordingly.
      void go_i_th( size_t an_id_mesh ) ;

      // Move side-iterator to the first position and
      // configure `::adjacent_localFEcell'(0) and `::adjacent_localFEcell'(1)
      // accordingly.
      void start( void ) ;
      
      // Is side-iterator position valid ? 
      bool is_valid( void ) const ;

      // Move side-iterator one position and
      // configure `::adjacent_localFEcell'(0) and `::adjacent_localFEcell'(1)
      // accordingly.
      void go_next( void ) ;
      
      bool is_periodic( void ) const ;

      // geometry of the current side
      GE_Mpolyhedron const* polyhedron( void ) const ;
      
      /* `i_adj'-th polyhedron of the current side, which is periodic 
          (a periodic side is a notion that gathers two polyhedra related by 
          a geometric transform); the meaning of `i_adj' is such that:  
             `::polyhedron'(`i_adj') is a guaranteed to be a face of 
             `::adjacent_localFEcell'(`i_adj')->polyhedron()
      */
      GE_Mpolyhedron const* polyhedron( size_t i_adj ) const ;
      
      /* geometric transform 
           `::polyhedron'(`i_adj_source') ---> `::polyhedron'(1-`i_adj_source')
      */
      GE_Transform const* periodic_transform( size_t i_adj_source ) const ;
      
      // level of refinement of the current mesh (the coarsest level is 0)
      size_t refinement_level( void ) const ;

      // color of the current side
      GE_Color const* color( void ) const  ;
      
      // color associated to the i_adj-th polyhedron of the current side, 
      // which is periodic
      GE_Color const* color( size_t i_adj ) const ;

      // identifier of the current side
      size_t mesh_id( void ) const ;

      // unit normal to side,
      // oriented from `::adjacent_localFEcell'(0)->polyhedron()
      // to            `::adjacent_localFEcell'(1)->polyhedron()
      GE_Vector const* normal( void ) const ;
      
      /* 
      unit normal to the `i_adj'-th polyhedron of the current side, 
      which is periodic; 
       `::normal'(0) being oriented from the inside to the outside 
                     relatively to `::adjacent_localFEcell'(0)->polyhedron(); 
       `::normal'(1) being oriented from the outside to the inside 
                     relatively to `::adjacent_localFEcell'(1)->polyhedron()
      thus making these two normals conceptually
      oriented from `::adjacent_localFEcell'(0)->polyhedron()
      to            `::adjacent_localFEcell'(1)->polyhedron()
      */
      GE_Vector const* normal( size_t i_adj ) const ;
      
      // signed distance from the current side to the finite volume center
      // of the `i_adj'-th adjacent cell (positive if this finite volume 
      // center is inside the half-plane delimited by `self' which contains
      // the `i_adj'-th adjacent cell, negative otherwise)
      double distance_to_adjacent_finite_volume_center( size_t i_adj ) const ;

      PDE_HelperFIC* helper_FIC( void ) const ;
      
   //-- Computations at a given point of the current mesh(101.5)

      // Set `pt' as the calculation point of `::adjacent_localFEcell(0)'
      // and `::adjacent_localFEcell(1)'.
      void set_calculation_point( GE_Point const* pt ) ;

      GE_Point const* calculation_point( void ) const ;
      
      void set_calculation_points( GE_Point const* pt_0,
                                   GE_Point const* pt_1 ) ;

      GE_Point const* calculation_point( size_t i_adj ) const ;
      
   //-- Local discrete variational problem definition(102)

      // Identify the space of trial solutions with `row_field' and the
      // space of weighting functions with `col_field' and
      // configure `::adjacent_localFEcell'(0) and `::adjacent_localFEcell'(1)
      // accordingly.
      void set_row_and_col_fields( PDE_DiscreteField const* row_field, 
                                   PDE_DiscreteField const* col_field  ) ;

      // field identifying the space of weighting functions if `sf' is 0 
      // or the space of trial solutions if `sf' is 1
      PDE_DiscreteField const* field( PDE_LocalFE::field_id sf ) const ;
      
      // for the space of weighting functions if `sf' is 0 and the space 
      // of trial solutions if `sf' is 1, number of basis functions whose 
      // support intersects the two cells adjacent to the current side
      size_t nb_basis_functions( PDE_LocalFE::field_id sf ) const ;
      
      // for the space of weighting functions, global numbers of the basis 
      // functions whose support intersects the two cells adjacent
      // to the current side
      size_t_vector const& row_field_node_connectivity( void ) const ;

      // for the space of weighting functions, local numbers in the `i_adj'
      // cell adjacent to the current side of the basis functions relevant
      // to the current  side (`PEL::bad_index()' for weighting functions whose 
      // support does not intersect this cell)
      size_t_vector const& row_field_node_side_2_cell_index( 
                                                         size_t i_adj ) const ;

      // for the space of trial solutions, global numbers of the basis 
      // functions whose support intersects the two cells adjacent
      // to the current side
      size_t_vector const& col_field_node_connectivity( void ) const ;

      // for the space of trial functions, local numbers in the `i_adj'
      // cell adjacent to the current side of the basis functions relevant
      // to the current side (`PEL::bad_index()' for weighting functions whose 
      // support does not intersect this cell)
      size_t_vector const& col_field_node_side_2_cell_index( 
                                                        size_t i_adj ) const ;
      
   //-- IP-iterator movement (IP standing for Integration Points)(103)

      // Move IP-iterator to the first IP of the quadrature rule provided
      // by `qrp' on the current side and configure 
      // `::adjacent_localFEcell'(0) and `::adjacent_localFEcell'(1)
      // accordingly.
      void start_IP_iterator( GE_QRprovider const* qrp ) ;

      // Move IP-iterator one position (within the current side) and configure 
      // `::adjacent_localFEcell'(0) and `::adjacent_localFEcell'(1)
      // accordingly. 
      void go_next_IP( void ) ;

      // Is IP-iterator position valid ?
      bool valid_IP( void ) const ;

      // surfacic weight of current IP  
      double weight_of_IP( void ) const ;
      
      double weight_of_IP( size_t i_adj ) const ;

      // absolute coordinates of current IP
      GE_Point const* coordinates_of_IP( void ) const ;
      
      GE_Point const* coordinates_of_IP( size_t i_adj ) const ;
      
   //-- Input - Output

      void print_current_mesh( std::ostream& os, size_t indent_width ) const ;

   protected: //--------------------------------------------------------
      
   private: //----------------------------------------------------------

      PDE_CursorFEside( void ) ;
     ~PDE_CursorFEside( void ) ;
      PDE_CursorFEside( PDE_CursorFEside const& other ) ;
      PDE_CursorFEside& operator=( PDE_CursorFEside const& other ) ;

      PDE_CursorFEside( PEL_Object* a_owner,
                        PDE_GridFE const* a_grid,
                        PDE_LocalFEcell* cfe_0,
                        PDE_LocalFEcell* cfe_1 ) ;

      // to simulate selective export
      friend PDE_CursorFEside* 
      PDE_DomainAndFields::create_CursorFEside( PEL_Object* a_owner ) const ;

   //-- Internals on the current mesh

      void compute_current_distances( void ) const ;

      void do_internals_on_current_mesh( void ) ;

      bool rejected_side( void ) const ;

      static PDE_FaceFE const* the_side( PEL_VectorIterator const* it ) ;

      void customize_providers( GE_QRprovider const* qrp ) ;
      
      static void determine_local( size_t_vector& n,
                                   size_t_vector const& n0,
                                   size_t_vector const& n1,
                                   size_t_vector& loc0,
                                   size_t_vector& loc1 ) ;
      
   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Attributes

      PDE_GridFE const* const GRID ;
      size_t GEO_STATE ;

      PEL_Vector* const EXCLUDE_COLORS ;
      
      PEL_Vector const* SIDES ;
      PEL_VectorIterator* IT ;
      PDE_FaceFE const* SIDE ;
      bool OK ;
      
      PDE_LocalFEcell* const CELL_FE_0 ;
      PDE_LocalFEcell* const CELL_FE_1 ;
      
      GE_Vector* const NORMAL0 ;
      GE_Vector* const NORMAL1 ;
      GE_Vector* const DUMMY_VEC ;

      GE_CustomizedQR_provider const* const QRP0 ;
      GE_CustomizedQR_provider const* const QRP1 ;

      PDE_DiscreteField const* DF_ROW ;
      PDE_DiscreteField const* DF_COL ;
      size_t_vector ROW_NODES ;
      size_t_vector COL_NODES ;
      size_t_vector ROW_LOC0 ;
      size_t_vector ROW_LOC1 ;
      size_t_vector COL_LOC0 ;
      size_t_vector COL_LOC1 ;

      bool SAME_ROW_AND_COL ;
      mutable doubleArray2D DIST ;
      
      mutable PDE_HelperFIC* HELPER_FIC ;
} ;

#endif
