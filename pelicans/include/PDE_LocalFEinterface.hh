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

#ifndef PDE_LOCAL_FE_INTERFACE_HH
#define PDE_LOCAL_FE_INTERFACE_HH

#include <PDE_LocalFEmulti.hh>

#include <PDE_DomainAndFields.hh> 

class PEL_Vector ;
class PEL_VectorIterator ;

class GE_Color ;
class GE_Matrix ;
class GE_Vector ;

class PDE_CellFE ;
class PDE_FaceFE ;
class PDE_DiscreteField ;
class PDE_GridFE ;

/*
Special PDE_LocalFE servers for which the geometrical entity is the
common face of two volumic meshes, one being the inward adjacent mesh and 
the other one the outward adjacent mesh.

INSTANCE DELIVERY : `PDE_DomainAndFields::create_LocalFEinterface' 
*/

class PEL_EXPORT PDE_LocalFEinterface : public PDE_LocalFEmulti
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
      
      // the unit normal vector directed from the inward adjacent volumic
      // mesh towards the outward adjacent volumic mesh
      GE_Vector const* outward_normal( void ) const ;

   //-- Computations at a given point of the current mesh

      virtual void set_calculation_point( GE_Point const* pt ) ;

   //-- Input - Output

      virtual void print_current_mesh( std::ostream& os, 
                                       size_t indent_width ) const ; 

   protected: //--------------------------------------------------------
      
   private: //----------------------------------------------------------

      PDE_LocalFEinterface( void ) ;
     ~PDE_LocalFEinterface( void ) ;
      PDE_LocalFEinterface( PDE_LocalFEinterface const& other ) ;
      PDE_LocalFEinterface& operator=( PDE_LocalFEinterface const& other ) ;

      friend 
      PDE_LocalFEinterface*
      PDE_DomainAndFields::create_LocalFEinterface( PEL_Object* a_owner,
                                      GE_Color const* inward_domain,
                                      GE_Color const* outward_domain ) const ;

      PDE_LocalFEinterface( PEL_Object* a_owner,
                            PDE_GridFE const* grid,
                            GE_Color const* inward_domain,
                            GE_Color const* outward_domain ) ;

   //-- Internals on the current mesh

      void do_internals_on_current_mesh( void ) ;

      virtual void compute_itg_pts( GE_QuadratureRule const* qr ) ;

      bool rejected_side( void ) const ;

      static PDE_FaceFE const* the_side( PEL_VectorIterator const* it ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Attributes

      PEL_VectorIterator* IT ;
      PDE_FaceFE const* SIDE ;
      PDE_CellFE const* CELL ;
      bool OK ;
      size_t NB_SIDES ;

      GE_Color const* IN_COLOR ;
      GE_Color const* OUT_COLOR ;
      PEL_Vector* INTERFACE ;

      // internal variables
      GE_Matrix* CELL_tr_JAC ;
      doubleArray3D* HESSIAN ;
      GE_Point* const CELL_PT_REF ;
      GE_Point* const PT ; 
      GE_Vector* NORMAL ;
     
} ;

#endif
