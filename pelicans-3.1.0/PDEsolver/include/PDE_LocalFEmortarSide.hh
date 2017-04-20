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

#ifndef PDE_LOCAL_FE_MORTAR_SIDE_HH
#define PDE_LOCAL_FE_MORTAR_SIDE_HH

#include <PDE_LocalFEmulti.hh>

#include <PDE_InterfaceAndFields.hh> 

class PEL_VectorIterator ;

class GE_Color ;
class GE_Matrix ;
class GE_Vector ;

class PDE_DiscreteField ;
class PDE_GridFE ;
class PDE_MortarSideFE ;

/*
INSTANCE DELIVERY : `PDE_InterfaceAndFields::create_LocalFEmortarSide'
*/

class PEL_EXPORT PDE_LocalFEmortarSide : public PDE_LocalFEmulti
{
   public: //-----------------------------------------------------------
  
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
      
   //-- Computations at a given point of the current mesh

      virtual void set_calculation_point( GE_Point const* pt ) ;

   //-- Input - Output

      virtual void print_current_mesh( std::ostream& os, 
                                       size_t indent_width ) const ; 

   protected: //--------------------------------------------------------
      
   private: //----------------------------------------------------------

      PDE_LocalFEmortarSide( void ) ;
     ~PDE_LocalFEmortarSide( void ) ;
      PDE_LocalFEmortarSide( PDE_LocalFEmortarSide const& other ) ;
      PDE_LocalFEmortarSide& operator=( 
                             PDE_LocalFEmortarSide const& other ) ;

      PDE_LocalFEmortarSide( PEL_Object* a_owner,
			     size_t nb_sp_dims, 
                             PEL_Vector const* mortar_sides ) ;

      // to simulate selective export
      friend PDE_LocalFEmortarSide* 
      PDE_InterfaceAndFields::create_LocalFEmortarSide( 
                                                 PEL_Object* a_owner ) const ;

   //-- Internals on the current mesh

      void do_internals_on_current_mesh( void ) ;

      virtual bool is_in_mesh( PDE_DiscreteField const* ff,
                               size_t i,
                               PDE_BasisFunction const* bf ) const ;

      void append_one_domain_local_fields( PEL_Vector const* domain_bounds ) ;

      virtual void compute_itg_pts( GE_QuadratureRule const* qr ) ;

      void connect_with_one_domain_meshes( bool is_IP,
                                           GE_Point const* pt,
                                           PEL_Vector const* domain_bounds ) ;

      bool rejected_side( void ) const ;

      static PDE_MortarSideFE const* the_side( PEL_VectorIterator const* it ) ;

   //-- Attributes

      size_t NB_SIDES ;
      PEL_VectorIterator* IT ;
      PDE_MortarSideFE const* SIDE ;
      bool OK ;
      GE_Point* const PT ;
      GE_Point* const BD_PT_REF ;
      GE_Point* const CELL_PT_REF ;
      GE_Matrix* CELL_tr_JAC ;
      doubleArray3D* CELL_HESSIAN ;
} ;

#endif
