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

#ifndef PDE_LOCAL_FE_BOUND_HH
#define PDE_LOCAL_FE_BOUND_HH

#include <PDE_LocalFEmulti.hh>

#include <PDE_DomainAndFields.hh>

class PEL_VectorIterator ;

class GE_Point ;
class GE_Matrix ;
class GE_SimplePolygon2D ;
class GE_Vector ;

class PDE_BoundFE ;
class PDE_GridFE ;

/*
Special PDE_LocalFE servers for which the geometrical entity is a 
boundary mesh of a volumic domain.

INSTANCE DELIVERY : `PDE_DomainAndFields::create_LocalFEbound'
*/

class PEL_EXPORT PDE_LocalFEbound : public PDE_LocalFEmulti
{
   public: //-----------------------------------------------------------
      
   //-- Instance delivery and initialization

      // Update the specific attributes of `self' in case of modification
      // of the underlying grid.
      virtual void update( void ) ;       
    
   //-- Mesh-iterator movement

      virtual size_t nb_meshes( void ) const ;

      virtual void go_i_th( size_t an_id_mesh ) ;
      
      virtual void start( void ) ;

      virtual bool is_valid( void ) const ;

      virtual void go_next( void ) ;

      virtual GE_Mpolyhedron const* polyhedron( void ) const ;

      virtual size_t refinement_level( void ) const ;

      virtual GE_Color const* color( void ) const  ;

      virtual size_t mesh_id( void ) const ;

      // the unit outward normal vector with respect to the bounded domain
      GE_Vector const* outward_normal( void ) const ;

      GE_Mpolyhedron const* adjacent_cell_polyhedron( void ) const ;

      GE_Color const* adjacent_cell_color( void ) const ;
      
      size_t adjacent_cell_id( void ) const ;

      // signed distance from the current bound to the finite volume center
      // of the adjacent cell (positive if this finite volume center is
      // inside the half-plane delimited by `self' which contains its
      // adjacent cell)
      double distance_to_adjacent_finite_volume_center( void ) const ;

   //-- Computations at a given point of the current mesh

      virtual void set_calculation_point( GE_Point const* pt ) ;

   //-- Input - Output

      virtual void print_current_mesh( std::ostream& os, 
                                       size_t indent_width ) const ;

   protected: //--------------------------------------------------------
      
   private: //----------------------------------------------------------

      PDE_LocalFEbound( void ) ;
     ~PDE_LocalFEbound( void ) ;
      PDE_LocalFEbound( PDE_LocalFEbound const& other ) ;
      PDE_LocalFEbound& operator=( PDE_LocalFEbound const& other ) ;

      friend PDE_LocalFEbound* PDE_DomainAndFields::create_LocalFEbound(
                                     PEL_Object* a_owner ) const ;

      PDE_LocalFEbound( PEL_Object* a_owner,
                        PDE_GridFE const* a_grid ) ;

   //-- Internals on the current mesh

      void compute_current_distance( void ) const ;

      void do_internals_on_current_mesh( void ) ;

      virtual bool is_in_mesh( PDE_DiscreteField const* ff,
                               size_t i,
                               PDE_BasisFunction const* bf ) const ;

      virtual GE_QuadratureRule const* quadrature_rule( 
                                        GE_QRprovider const* qrp ) const ;

      virtual void compute_itg_pts( GE_QuadratureRule const* qr ) ;

      bool rejected_bound( void ) const ;

      static PDE_BoundFE const* the_bound( PEL_VectorIterator const* it ) ;

   //-- Attributes

      // Mesh iterator datas :
      PEL_Vector const* BOUNDS ;
      PEL_VectorIterator* IT ;
      PDE_BoundFE const* BOUND ;
      bool OK ;

      // Grid :
      PDE_GridFE const* const GRID ;
      size_t GEO_STATE ;

      // internal variables
      GE_Point* const PT ;
      GE_Point* const BD_PT_REF ;
      GE_Point* const CELL_PT_REF ;
      GE_Matrix* CELL_tr_JAC ;
      doubleArray3D* HESSIAN ;

      GE_Vector* const DUMMY_VEC ;
      mutable doubleVector DIST ;
} ;

#endif
