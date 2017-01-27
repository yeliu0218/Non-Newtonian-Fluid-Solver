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

#ifndef PDE_DOMAIN_BUILDER_HH
#define PDE_DOMAIN_BUILDER_HH

#include <PEL_Object.hh>

#include <size_t_array2D.hh>
#include <size_t_vector.hh>
#include <stringVector.hh>

#include <map>

class PEL_List ;
class PEL_ListIterator ;
class PEL_ModuleExplorer ;
class PEL_Vector ;
class PEL_VectorIterator ;

class GE_Color ;
class GE_Meshing ;
class GE_Mpolyhedron ;
class GE_Point ;
class GE_SetOfPoints ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_GridFE ;
class PDE_IDomainOnGrid ;
class PDE_ReferenceElement ;
class PDE_SetOfBasisFunctions ;
class PDE_SetOfBCs ;
class PDE_SetOfDiscreteFields ;
class PDE_SetOfFieldCompositions ;

#include <string>
#include <vector>

class PEL_EXPORT PDE_DomainBuilder : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_DomainBuilder* object( std::string const& a_name ) ;

      static PDE_DomainBuilder* create( PEL_Object* a_owner, 
                                        PEL_ModuleExplorer const* exp,
                                        std::string const& a_name  ) ;

   //-- Characteristics

      std::string const& name( void ) const ;

   //-- Modifiers

      // Create a new discrete field called `name_of_new_field', identical
      // to the already existing discrete field called `model_name',
      // and add it to `::set_of_discrete_fields()'.
      void duplicate_field( std::string const& model_name,
                            std::string const& name_of_new_field ) const ;
      
      // Create a new discrete field called `name_of_new_field', with the
      // same discretization than the already existing discrete field called
      // `model_name', and add it to `::set_of_discrete_fields()'. 
      // The values of all the DOFs of the newly created discrete field 
      // (of name `name_of_new_field') are set to zero.
      void duplicate_field( std::string const& model_name,
                            std::string const& name_of_new_field,
                            size_t nb_components_of_new_field,
                            size_t storage_depth_of_new_field ) const ;

      // Create new discrete field according to the Module Hierarchy
      // attached to `exp', and add them to `::set_of_discrete_fields()'.
      void append_fields( PEL_ModuleExplorer const* exp ) ;

   //-- Associated Module Hierarchy

      bool has_explorer( std::string const& path_and_name ) const ;

      PEL_ModuleExplorer* create_explorer( 
                                  PEL_Object* a_owner,
                                  std::string const& path_and_name ) const ;

   //-- Associated PDE_DomainAndField instance

      void set_facade( PDE_DomainAndFields* a_dom ) ;

      PDE_DomainAndFields* facade( void ) const ;

   //-- Conformal adjacencies

      size_t nb_conformal_adjacencies( void ) const ;

      PDE_DomainBuilder* conformal_adjacency( size_t i ) const ;

      void insert_conformal_adjacency( PDE_DomainBuilder* other ) ;
 
   //-- Elements also accessible through PDE_DomainAndFields instance

      size_t nb_space_dimensions( void ) const ;

      PDE_SetOfDiscreteFields* set_of_discrete_fields( void ) const ;

      PDE_SetOfFieldCompositions const* set_of_field_compositions( void ) const ;
      
      bool is_defined_on_boundary( std::string const& fname ) const ;

      PDE_SetOfBCs const* set_of_boundary_conditions( void ) const ;

      GE_SetOfPoints* set_of_vertices( void ) const ;
      
      PDE_SetOfBasisFunctions* set_of_basis_functions( void ) const ;
      
   //-- Inner PELICANS objects

      PDE_GridFE* finite_element_grid( void ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_DomainBuilder( void ) ;
     ~PDE_DomainBuilder( void ) ;
      PDE_DomainBuilder( PDE_DomainBuilder const& other ) ;
      PDE_DomainBuilder& operator=( PDE_DomainBuilder const& other ) ;
     
      PDE_DomainBuilder( PEL_Object* a_owner,
                         PEL_ModuleExplorer const* exp,
                         std::string const& a_name ) ;

      void build_all( void ) ;
      
      GE_Meshing* create_meshing( PEL_Object* a_owner ) const ;

      void build_set_of_discrete_fields( PEL_ModuleExplorer const* exp ) ;
      
      void build_grid( GE_Meshing* tria ) ;
      
      void build_special_colors( void ) ;

      void build_discretizations( PEL_ModuleExplorer const* exp ) ;
      
      void build_field_compositions( void ) ;
      
      void build_one_interior_field( PDE_DiscreteField* ff,
                                     PEL_Vector const* elms ) ;

      void build_one_boundary_field( PDE_DiscreteField* ff,
                                     PEL_Vector const* elms ) ;
      
      PEL_Vector* create_elements( PEL_Object* a_owner,
                                   PEL_ModuleExplorer const* ee ) const ;
      
      PDE_ReferenceElement const* ref_element(
                                      std::string const& field_name,
                                      PEL_Vector const* elms,
                                      GE_Mpolyhedron const* poly ) const ;

      void build_field_nodes( PDE_DiscreteField const* ff,
                              PEL_Vector const* elms,
                              size_t& nb_nodes,
                              size_t_vector& nb_meshnodes,
                              size_t_array2D& node_of_mesh ) const ;

   //-- Class Attributes & Methods

      static std::map< std::string, PDE_DomainBuilder* >& INSTANCES( void ) ;

   //-- Attributes

      PEL_ModuleExplorer const* EXP ;
      std::string NAME ;
      PDE_DomainAndFields* DOM ;
      int VERB ;
      size_t DIM ;
      PDE_GridFE* GRID ;
      PDE_SetOfDiscreteFields* FIELDS ;
      PDE_SetOfFieldCompositions* FIELD_COMPOS ;
      PDE_SetOfBCs const* BCS ;
      PDE_SetOfBasisFunctions* BFS ;

      size_t NVERT ;
      GE_SetOfPoints* VMAN ;
      stringVector FIELDS_ON_BOUNDARY ;
      
      GE_Point* W_PT ;

      std::vector< PDE_DomainBuilder* > ADJ ;
} ;

#endif
