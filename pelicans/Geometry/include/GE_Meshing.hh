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

#ifndef GE_MESHING_HH
#define GE_MESHING_HH

#include <PEL_Object.hh>
#include <doubleVector.hh>
#include <size_t_vector.hh>

class PEL_DataWithContext ;
class PEL_Double ;
class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;

class GE_Color ;
class GE_Colorist ;
class GE_ReferencePolyhedron ;

class stringVector ;
class size_t_array2D ;

/*
Meshings of geometric domains.

A Meshing is a collection of finitely many cells such that :
   - the union of all cells is the closure of the domain ;
   - each cell is a polyhedron that can represented by a
     `GE_Mpolyhedron::' instance ;
   - two distinct cells have disjoint interiors ;
   - any face of any cell is either a face of another cell or a subset
     of the domain boundary.

Three sets of geometrical entities are involved in meshings :
   - the set of cells ;
   - the set of cell vertices ;
   - the set of the cell faces.

Each geometrical entity in these sets is identified with an intance
of `GE_Color::' called its color.

Each of these sets might be traversed by an iterator that visits every
geometrical entity exactly once. The traversal order defines an implicit 
numbering of the geometrical entities within their sets, numbering
used in the representation of the connections betwenn cells and vertices
and between cells and faces.

FRAMEWORK INSTANTIATION :
   1. Derive a concrete subclass, say `MyMeshing' (derivation of abstract
      subclasses is improbable although posible).
   2. Choose a name for `MyMeshing', say "MyMeshing".
   3. Implement a private destructor.
   4. Declare all constructors private.
   5. Define the prototype to be registered :
      5.1 Implement a default constructor that initializes the
          `GE_Meshing::' subobject by calling
               `GE_Meshing( std::string const& )'
          with "MyMeshing" as argument.
      5.2 Define and initialize a static instance by calling the default
          constructor.
             declaration :
                `static MyMeshing const* prototype ;'
             definition :
                `MyMeshing const* MyMeshing::prototype = new MyMeshing() ;'
   6. Implement the `create_replica' methods, that call a constructor.
   7. Implemenent the constructor called by `create_replica'. The `MyMeshing'
      subobject is initialized by calling
      `GE_Meshing( PEL_Object*, PEL_ModuleExplorer const*, size_t )'
      with "MyMeshing" as a second argument.
   8. Implement all pure virtual methods.
*/

class PEL_EXPORT GE_Meshing : public PEL_Object
{
   public: //------------------------------------------------------------
      
   //-- Instance delivery and initialization

      // Considering the module hierarchy associated to `exp', create and 
      // return an instance of a derived class whose name is the data
      // associated to the keyword "concrete_name" in `exp'->module()
      static GE_Meshing* create( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp,
                                 size_t dim_space ) ;
      
   //-- Measurement

      // number of space dimension of the geometric domain
      size_t nb_space_dimensions( void ) const ;

      // size of the set of cell vertices
      virtual size_t nb_vertices( void ) const = 0 ;

      // number of polyhedron subdividing the geometric domain
      virtual size_t nb_cells( void ) const = 0 ;

      // size of the set of cell faces
      virtual size_t nb_faces( void ) const = 0 ;
      
   //-- Vertex-iterator movement

      // Move vertex-iterator to first position.
      virtual void start_vertex_iterator( void ) = 0 ;

      // Is vertex-iterator position valid ?
      virtual bool valid_vertex( void ) const = 0 ;

      // Move vertex-iterator one position.
      virtual void go_next_vertex( void ) = 0 ;

      // coordinates of the current vertex
      virtual doubleVector const& vertex_coordinates( void ) const = 0 ;

      // color of the current vertex, given by `GE_Colorist::vertex_color' if 
      // a `GE_Colorist::' object has been associated to self (in the
      // `::create' method), or by `::default_vertex_color' otherwise
      GE_Color const* vertex_color( void ) const ;

   //-- Cell-iterator movement

      // Move cell-iterator to first position.
      virtual void start_cell_iterator( void ) = 0 ;

      // Is cell-iterator position valid ?
      virtual bool valid_cell( void ) const = 0 ;

      // Move cell-iterator one position.
      virtual void go_next_cell( void ) = 0 ;

      // polyhedron name of the current cell
      virtual std::string const& cell_polyhedron_name( void ) const = 0 ;

      // reference polyhedron associated to the current cell
      virtual GE_ReferencePolyhedron const* cell_reference_polyhedron( 
                                                               void ) const ;

      // color of the current cell, given by `GE_Colorist::cell_color' if 
      // a `GE_Colorist::' object has been associated to self (in the
      // `::create' method), or by `::default_cell_color' otherwise
      GE_Color const* cell_color( void ) const ;

      // indices of the current cell vertices, using the implicit numbering
      // defined by the traversal order of the vertex-iterator
      virtual size_t_vector const& cell_vertices( void ) const = 0 ;

      // indices of the current cell faces, using the implicit numbering
      // defined by the traversal order of the face-iterator
      virtual size_t_vector const& cell_faces( void ) const = 0 ;

   //-- Face-iterator movement

      // Move face-iterator to first position.
      virtual void start_face_iterator( void ) = 0 ;

      // Is face-iterator position valid ?
      virtual bool valid_face( void ) const = 0 ;

      // Move face-iterator one position.
      virtual void go_next_face( void ) = 0 ;

      // polyhedron name of the current face
      virtual std::string const& face_polyhedron_name( void ) const = 0 ;

      // reference polyhedron associated to the current face
      virtual GE_ReferencePolyhedron const* face_reference_polyhedron( 
                                                               void ) const ;

      // color of the current face, given by `GE_Colorist::face_color' if 
      // a `GE_Colorist::' object has been associated to self (in the
      // `::create' method), or by `::default_face_color' otherwise
      GE_Color const* face_color( void ) const ;

      // indices of the current face vertices, using the implicit numbering
      // defined by the traversal order of the vertex-iterator
      virtual size_t_vector const& face_vertices( void ) const = 0 ;
            
   //-- Input - Output
     
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
      // Write a textual representation of `self' to `os' with
      // `indent_width' indentation.
      virtual void display( std::ostream& os, size_t indent_width ) ;
      
      // Write `self' so as to be an entry for `GE_ExplicitMeshing::'.
      // Bounds are not described.
      void print_as_an_explicit_meshing( std::ostream& os ) ;

   protected: //---------------------------------------------------------
      
      virtual ~GE_Meshing( void ) ;

      // for prototype registration only
      GE_Meshing( std::string const& name ) ;

      GE_Meshing( PEL_Object* a_owner,
                  PEL_ModuleExplorer const* exp,
                  size_t dim_space ) ;
      
      // Is `self' a prototype ?
      bool is_a_prototype( void ) const ;

      // Create and return a replica of `self' on the basis of `exp'.
      virtual GE_Meshing* create_replica( PEL_Object* a_owner,
                                          PEL_ModuleExplorer const* exp,
                                          size_t dim_space ) const = 0 ;

      // current vertex color if `self' has not been associated to a
      // `GE_Colorist::' object in the `::create' method
      virtual GE_Color const* default_vertex_color( void ) const = 0 ;

      // current cell color if `self' has not been associated to a
      // `GE_Colorist::' object in the `::create' method
      virtual GE_Color const* default_cell_color( void ) const = 0 ;

      // current face color if `self' has not been associated to a
      // `GE_Colorist::' object in the `::create' method
      virtual GE_Color const* default_face_color( void ) const = 0 ;

      // Read the argument of the entry whose keyword is mesh_polyhedron
      // in the module associated to `exp' and raise a fatal error if the 
      // dimensions are not consistent with `dim_space' space dimensions.
      void read_mesh_polyhedron( PEL_ModuleExplorer const* exp,
                                 size_t dim_space,
                                 std::string& face_poly_name,
                                 std::string& cell_poly_name ) const ;

      // Check if the meshes of names `cell_poly_name' and `face_poly_name' are
      // consistent with the meshing to be built.
      // IMPLEMENTATION: do nothing.
      virtual void check_mesh_polyhedron(
                         size_t dim_space,
                         std::string const& face_poly_name,
                         std::string const& cell_poly_name ) const ;

      // Raise an fatal error with a message saying that `dim' does not
      // belong to `allowed_dims'.
      void raise_invalid_nb_space_dimensions(
                                size_t dim, std::string allowed_dims ) const ;

      // Raise an fatal error with a message saying that `poly_name' does not
      // represent a valid cell polyhedron.
      void raise_invalid_cell_polyhedron(
                         std::string const& poly_name,
                         std::string const& error = "" ) const ;
      
      // Raise an fatal error with a message saying that `poly_name' does not
      // represent a valid face polyhedron.
      void raise_invalid_face_polyhedron(
                         std::string const& poly_name,
                         std::string const& error = "" ) const ;

   //-- Roundoff strategy

      void initialize_rounding_strategy( PEL_ModuleExplorer const* exp ) ;

      double roundoff( double x ) const ;

   //-- Preconditions, Postconditions, Invariant

      bool go_next_vertex_PRE( void ) const ;

      bool vertex_coordinates_PRE( void ) const ;
      bool vertex_coordinates_POST( doubleVector const& result ) const ;

      bool default_vertex_color_PRE( void ) const ;
      bool default_vertex_color_POST( GE_Color const* result ) const ;

      bool start_cell_iterator_PRE( void ) const ;
      bool start_cell_iterator_POST( void ) const ;

      bool go_next_cell_PRE( void ) const ;

      bool valid_cell_POST( bool result ) const ;

      bool cell_polyhedron_name_PRE( void ) const ;
      bool cell_polyhedron_name_POST( std::string const& result ) const ;

      bool cell_reference_polyhedron_PRE( void ) const ;
      bool cell_reference_polyhedron_POST( 
                               GE_ReferencePolyhedron const* result ) const ;

      bool default_cell_color_PRE( void ) const ;
      bool default_cell_color_POST( GE_Color const* result ) const ;

      bool cell_vertices_PRE( void ) const ;
      bool cell_vertices_POST( size_t_vector const& result ) const ;

      bool cell_faces_PRE( void ) const ;
      bool cell_faces_POST( size_t_vector const& result ) const ;

      bool start_face_iterator_PRE( void ) const ;
      bool start_face_iterator_POST( void ) const ;

      bool go_next_face_PRE( void ) const ;

      bool valid_face_POST( bool result ) const ;

      bool face_polyhedron_name_PRE( void ) const ;
      bool face_polyhedron_name_POST( std::string const& result ) const ;

      bool face_reference_polyhedron_PRE( void ) const ;
      bool face_reference_polyhedron_POST( 
                               GE_ReferencePolyhedron const* result ) const ;

      bool default_face_color_PRE( void ) const ;
      bool default_face_color_POST( GE_Color const* result ) const ;

      bool face_vertices_PRE( void ) const ;
      bool face_vertices_POST( size_t_vector const& result ) const ;

      bool check_mesh_polyhedron_PRE( 
                         size_t dim_space,
                         std::string const& face_poly_name,
                         std::string const& cell_poly_name ) const ;

      virtual bool create_replica_PRE( PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp,
                                       size_t dim_space ) const ;
      virtual bool create_replica_POST( GE_Meshing const* result,
                                        PEL_Object* a_owner,
                                        size_t dim_space ) const ;

   private: //-----------------------------------------------------------

      GE_Meshing( void ) ;
      GE_Meshing( GE_Meshing const& other ) ;
      GE_Meshing& operator=( GE_Meshing const& other ) ;

      void raise_invalid_mesh_polyhedron( std::string const& error ) const ;

      static void print_colors( std::ostream& os, 
                                stringVector const& colors, 
                                size_t_array2D const& item_color ) ;

      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Attributes

      size_t const DIM ;
      bool const IS_PROTO ;
      GE_Colorist* COLORIST ;
      PEL_DataWithContext* RDOFF ;
      PEL_Double* XX ;
} ;

#endif


