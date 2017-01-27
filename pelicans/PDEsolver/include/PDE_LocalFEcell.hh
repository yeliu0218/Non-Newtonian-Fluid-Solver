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

#ifndef PDE_LOCAL_FE_CELL_HH
#define PDE_LOCAL_FE_CELL_HH

#include <PDE_LocalFEsingle.hh>

#include <PDE_DomainAndFields.hh>

class GE_Matrix ;
class GE_Point ;

class PEL_ModuleExplorer ;
class PEL_VectorIterator ;

class PDE_CellFE ;
class PDE_CFootFinder ;
class PDE_GridFE ;
class PDE_IDomainOnGrid ;
class PDE_ReferenceElement ;

class size_t_vector ;

/*
Special PDE_LocalFE servers for which the geometrical entity is a 
volumic mesh.

INSTANCE DELIVERY : `PDE_DomainAndFields::create_LocalFEcell'
                    `PDE_CursorFEside::adjacent_localFEcell'

*/

class PEL_EXPORT PDE_LocalFEcell : public PDE_LocalFEsingle
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
      
      // identifiers of the cells that share a common side with the 
      // current cell
      // (may be used as argument when calling `PDE_LocalFEcell::go_i_th()')
      size_t_vector const& adjacent_cell_ids( void ) const ;
      
      // identifiers of the sides that are adjacent to the current cell
      // (may be used as argument when calling `PDE_CursorFEside::go_i_th()')
      size_t_vector const& adjacent_side_ids( void ) const ;

      // identifiers of the bounds that are adjacent to the current cell,
      // if any
      // (may be used as argument when calling `PDE_LocalFEbound::go_i_th()')
      size_t_vector const& adjacent_bound_ids( void ) const ;

   //-- Computations at a given point of the current mesh

      virtual void set_calculation_point( GE_Point const* pt ) ;

   //-- Transport

      void set_foot_finder( PEL_ModuleExplorer const* exp ) ;

      bool has_foot_finder( void ) const ;

      PDE_CFootFinder* foot_finder( void ) const ;

      void synchronize_foot_finder( void ) ;
      
   //-- Input - Output

      virtual void print_current_mesh( std::ostream& os, 
                                       size_t indent_width ) const ;

   //-- Hidden

      PDE_ReferenceElement const* reference_element(
                                   PDE_DiscreteField const* ff ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_LocalFEcell( void ) ;
     ~PDE_LocalFEcell( void ) ;
      PDE_LocalFEcell( PDE_LocalFEcell const& other ) ;
      PDE_LocalFEcell& operator=( PDE_LocalFEcell const& other ) ;

      friend PDE_LocalFEcell* PDE_DomainAndFields::create_LocalFEcell(
                                     PEL_Object* a_owner ) const ;

      PDE_LocalFEcell( PEL_Object* a_owner,
                       PDE_GridFE const* a_grid ) ;

   //-- Internals on the current mesh

      void do_internals_on_current_mesh( void ) ;
      
      virtual bool is_in_mesh( PDE_BasisFunction const* bf ) const ;

      virtual GE_QuadratureRule const* quadrature_rule( 
                                          GE_QRprovider const* qrp ) const ;

      virtual void compute_itg_pts( GE_QuadratureRule const* qr ) ;

      bool rejected_cell( void ) const;

      static PDE_CellFE* the_cell( PEL_VectorIterator const* it ) ;

   //-- Attributes
      
      // Mesh iterator datas :
      PEL_Vector const* CELLS ;
      PEL_VectorIterator* IT ;
      PDE_CellFE* CELL ;
      bool OK ;

      // Grid :
      PDE_GridFE const* const GRID ;
      size_t GEO_STATE ;

      // Foot finder :
      PDE_CFootFinder* FFINDER ;

      // internal variables
      GE_Matrix* tr_JAC ;
      doubleArray3D* HESSIAN ;
      GE_Point* const PT ;
      GE_Point* const PT_REF ;
      mutable size_t_vector ADJ_CELL ;
      mutable size_t_vector ADJ_SIDE ;
      mutable size_t_vector ADJ_BOUND ;
} ;


#endif

