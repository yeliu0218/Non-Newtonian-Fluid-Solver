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

#ifndef PDE_LOCAL_FE_HH
#define PDE_LOCAL_FE_HH

#include <PEL_Object.hh>

#include <vector>

class PEL_Vector ;
class doubleArray2D ;
class doubleArray3D ;
class doubleVector ;
class size_t_vector ;

class GE_Color ;
class GE_Matrix ;
class GE_Point ;
class GE_Mpolyhedron ;
class GE_QuadratureRule ;
class GE_QRprovider ;

class PDE_DiscreteField ;

/*
Servers, which always have a current position on a geometrical entity
accessible through an iterator, for building Discrete Variational Problems. 
The geometrical entity is identified by a mesh, whose nature depends on the 
concrete subclasses.

The statement of the discrete Variational Problem is : 
   find u in U such that for all v in V, a given condition involving 
   u and v is fulfilled, where U (called space of trial solutions) 
   and V (called space of weighting functions) are finite dimensional 
   functional spaces spanned by basis functions of local support.

Instances of PDE_DiscreteField are handled, each of them being
associated to a Lagrange finite element. These finite elements are 
providing the basis functions that span the spaces of trial 
solutions and weighting functions. The basis functions are real
valued. For non scalar functional spaces and fields with multiple 
components, the tensor product of the scalar basis functions related 
to each component is implicitly performed. 
*/

class PEL_EXPORT PDE_LocalFE : public PEL_Object
{

   public: //-----------------------------------------------------------

   //-- Indicators(80)

      enum field_id { row = 0, col = 1 } ;

      // indicator for node discretization
      // (implementation : node = 0 ) 
      static int const node ;
      
      // indicator for zero-th order of spatial derivation
      // (implementation : N = 1 ) 
      static int const N ;

      // indicator for first order of spatial derivation
      // (implementation : dN = 2 ) 
      static int const dN ;

      // indicator for second order of spatial derivation
      // (implementation : d2N = 4 ) 
      static int const d2N ;

   //-- Geometry(90)

      // number of space dimensions
      size_t nb_space_dimensions( void ) const ;
      
   //-- Handled fields(95)

      // Notify that the `order'-th spatial derivatives of  
      // `ff' basis functions will be requested during mesh iteration.
      void require_field_calculation( PDE_DiscreteField const* ff,
                                      int order ) ;
      
      // Will the `order'-th spatial derivatives of `ff' basis function be
      // available during mesh iteration ?
      bool field_calculation_is_handled( PDE_DiscreteField const* ff,
                                         int order ) const ;

      // Is `ff' handled ? 
      bool field_is_handled( PDE_DiscreteField const* ff ) const ;

      // number_of_handled fields
      size_t nb_handled_fields( void ) const ;

      // iF-th handled field
      PDE_DiscreteField const* handled_field( size_t iF ) const ;
   
   //-- Mesh-iterator movement(100)

      // number of meshes that the mesh-iterator can traverse
      virtual size_t nb_meshes( void ) const = 0 ;

      // Exclude meshes for which color is `a_color'.
      void exclude_color( GE_Color const* a_color ) ;

      // Include excluded meshes for which color is `a_color'.
      void include_color( GE_Color const* a_color ) ;

      // Is `a_color' excluded ?
      bool is_excluded( GE_Color const* a_color ) const ;

      // Move mesh-iterator to the mesh of identifier `an_id_mesh'.
      virtual void go_i_th( size_t an_id_mesh ) = 0 ;
      
      // Move mesh-iterator to the first position.
      virtual void start( void ) = 0 ;

      // Is mesh-iterator position valid ? 
      virtual bool is_valid( void ) const = 0 ;

      // Move mesh-iterator one position.
      virtual void go_next( void ) = 0 ;

      // geometry of the current mesh
      virtual GE_Mpolyhedron const* polyhedron( void ) const = 0 ;
      
      // level of refinement of the current mesh (the coarsest level is 0)
      virtual size_t refinement_level( void ) const = 0 ;

      // color of the current mesh
      virtual GE_Color const* color( void ) const = 0 ;

      // identifier of the current mesh
      virtual size_t mesh_id( void ) const = 0 ;

   //-- Local discretization of fields(101)

      // number of nodes of `ff' involved in the computation, 
      // at any arbitrary point in the current mesh, of the reconstruction
      // defined by `ff' and its associated finite element
      virtual size_t nb_local_nodes( PDE_DiscreteField const* ff ) const = 0 ;

      // pointer to an object whose value is the location of the node 
      // associated to the `i'-th basis function involved 
      // in the computation, at any arbitrary point in the current mesh, 
      // of the reconstruction defined by `ff' and its associated 
      // finite element
      // WARNING
      //    the value of the object pointed by the result of this
      //    method might change each time this latter is called : the same 
      //    object is used whatever `ff' and `i', only its value changes
      virtual GE_Point const* local_node_location( PDE_DiscreteField const* ff,
                                                   size_t i ) const = 0 ;

      // Is the node of `ff' associated to the `i'-th basis function (involved 
      // in the computation, at any arbitrary point in the current mesh, 
      // of the reconstruction defined by `ff' and its associated 
      // finite element) located in the current mesh ?
      virtual bool local_node_is_in_mesh( PDE_DiscreteField const* ff,
                                          size_t i ) const = 0 ;

      // node of `ff' associated to the `i'-th basis function involved 
      // in the computation, at any arbitrary point in the current mesh, 
      // of the reconstruction defined by `ff' and its associated 
      // finite element
      virtual size_t global_node( PDE_DiscreteField const* ff,
                                  size_t i ) const = 0 ;
      
      // refinement level of the node of `ff' associated to the `i'-th basis 
      // function involved in the computation, at any arbitrary point in 
      // the current mesh, of the reconstruction defined by `ff' and its 
      // associated finite element (0 being the coarsest refinement level)
      virtual size_t node_refinement_level( PDE_DiscreteField const* ff,
                                            size_t i ) const = 0 ;

   //-- Computations at a given point of the current mesh(101.5)

      // Set `pt' as the calculation point for subsequent calls to
      // `::value_at_pt', `::gradient_at_pt', `::hessian_at_pt', 
      // `::N_at_pt', `::dN_at_pt', `::d2N_at_pt'.
      virtual void set_calculation_point( GE_Point const* pt ) = 0 ;

      virtual GE_Point const* calculation_point( void ) const = 0 ;

      // value at the point given by `::set_calculation_point' of the `ic'-th 
      // component of the reconstruction obtained from the `level'-th storage 
      // of `ff' and its associated finite element
      virtual double value_at_pt( PDE_DiscreteField const* ff,
                                  size_t level,
                                  size_t ic=0 ) const = 0 ;

      // value at the point given by `::set_calculation_point' of the 
      // derivative with respect to the `a'-th direction of the `ic'-th 
      // component of the reconstruction obtained from the `level'-th storage 
      // of `ff' and its associated finite element 
      virtual double gradient_at_pt( PDE_DiscreteField const* ff,
                                     size_t level,
                                     size_t a,
                                     size_t ic=0 ) const = 0 ;

      // value at the point given by `::set_calculation_point' of the double
      // derivative with respect to the `a'-th and `b'-th directions of 
      // the `ic'-th component of the reconstruction obtained from 
      // the `level'-th storage of `ff' and its associated finite element 
      virtual double hessian_at_pt( PDE_DiscreteField const* ff, 
                                    size_t level,
                                    size_t a,
                                    size_t b,
                                    size_t ic=0 ) const = 0 ;

      // value at the point given by `::set_calculation_point' of the 
      // `i'-th basis function (local numbering) in 
      // the finite element space of `ff'
      virtual double N_at_pt( PDE_DiscreteField const* ff, 
                              size_t i ) const = 0 ;

      // value at the point given by `::set_calculation_point' of the
      // derivative with respect to the `a'-th direction of the 
      // `i'-th basis function (local numbering) in the finite element 
      // space of `ff'
      virtual double dN_at_pt( PDE_DiscreteField const* ff, 
                               size_t i, 
                               size_t a ) const = 0 ;

      // value at the point given by `::set_calculation_point' of the 
      // double derivative with respect to the `a'-th and `b'-th directions
      //  of the `i'-th basis function (local numbering) in the finite element 
      // space of `ff'
      virtual double d2N_at_pt( PDE_DiscreteField const* ff,
                                size_t i, 
                                size_t a, 
                                size_t b ) const = 0 ;

   //-- Local discrete variational problem definition(102)

      // Identify the space of trial solutions with `row_field' and the
      // space of weighting functions with `col_field'.
      virtual void set_row_and_col_fields( 
                                   PDE_DiscreteField const* row_field, 
                                   PDE_DiscreteField const* col_field ) = 0 ;

      // field identifying the space of weighting functions if `sf' is 0 
      // or the space of trial solutions if `sf' is 1
      virtual PDE_DiscreteField const* field( field_id sf ) const = 0 ;

      // for the space of weighting functions if `sf' is 0 and the space 
      // of trial solutions if `sf' is 1, number of basis functions whose 
      // support intersects the geometrical entity defined by the current
      // mesh 
      virtual size_t nb_basis_functions( field_id sf ) const = 0 ;

      // for the space of weighting functions, global numbers of the basis 
      // functions whose support intersects the geometrical entity defined 
      // by the current mesh 
      virtual size_t_vector const& row_field_node_connectivity( void ) const = 0 ;

      // for the space of trial solutions, global numbers of the basis 
      // functions whose support intersects the geometrical entity defined 
      // by the current mesh 
      virtual size_t_vector const& col_field_node_connectivity( void ) const = 0 ;

   //-- IP-iterator movement (IP standing for Integration Points)(103)

      // Move IP-iterator to the first IP of the quadrature rule provided
      // by `qrp' on the current geometrical entity.  
      virtual void start_IP_iterator( GE_QRprovider const* qrp ) = 0 ;

      // Move IP-iterator one position (within the current geometrical
      // entity).  
      virtual void go_next_IP( void ) = 0 ;

      // Is IP-iterator position valid ?  
      virtual bool valid_IP( void ) const = 0 ;

      // weight of current IP  
      virtual double weight_of_IP( void ) const = 0 ;

      // absolute coordinates of current IP  
      virtual GE_Point const* coordinates_of_IP( void ) const = 0 ;

      // value at current IP of the `ic'-th component of the reconstruction
      // obtained from the `level'-th storage of `ff' and its
      // associated finite element
      virtual double value_at_IP( PDE_DiscreteField const* ff,
                                  size_t level,
                                  size_t ic=0 ) const = 0 ;

      // value at current IP of the derivative with respect to the `a'-th
      // direction of the `ic'-th component of the reconstruction
      // obtained from the `level'-th storage of `ff' and its
      // associated finite element 
      virtual double gradient_at_IP( PDE_DiscreteField const* ff, 
                                     size_t level,
                                     size_t a, 
                                     size_t ic=0 ) const = 0 ;

      // value at current IP of the double derivative with respect to 
      // the `a'-th and `b'-th directions of the `ic'-th component of 
      // the reconstruction obtained from the `level'-th storage of `ff' 
      // and its associated finite element 
      virtual double hessian_at_IP( PDE_DiscreteField const* ff, 
                                    size_t level,
                                    size_t a,
                                    size_t b,
                                    size_t ic=0 ) const = 0 ;

      // value at current IP of the `i'-th basis function (local numbering)
      // of the space of weighting functions (resp. trial solutions) 
      // if `sf' is 0 (resp. 1)  
      virtual double N_at_IP( field_id sf, size_t i ) const = 0 ;

      // value at current IP of the derivative with respect to the
      // `a'-th direction of the `i'-th basis function (local numbering)
      // of the space of weighting functions (resp. trial solutions) 
      // if `sf' is 0 (resp. 1)  
      virtual double dN_at_IP( field_id sf, size_t i, size_t a ) const = 0 ;

      // value at current IP of the double derivative with respect 
      // to the `a'-th and `b'-th directions of the `i'-th basis function 
      // (local numbering) of the space of weighting functions 
      // (resp. trial solutions) if `sf' is 0 (resp. 1)  
      virtual double d2N_at_IP( field_id sf, 
                                size_t i, size_t a, size_t b ) const = 0 ;

      // value at current IP of all the basis functions 
      // (organized according to their local numbering)
      // of the space of weighting functions (resp. trial solutions) 
      // if `sf' is 0 (resp. 1)  
      virtual doubleVector const& Ns_at_IP( field_id sf ) const = 0 ;
      
      // value at current IP of the derivative with respect to all
      // directions of all the basis functions of the space of 
      // weighting functions (resp. trial solutions) if `sf' is 0 (resp. 1)  
      virtual doubleArray2D const& dNs_at_IP( field_id sf ) const = 0 ;
      
      // value at current IP of the double derivative with respect 
      // to all directions of all basis functions of the space of 
      // weighting functions (resp. trial solutions) if `sf' is 0 (resp. 1)  
      virtual doubleArray3D const& d2Ns_at_IP( field_id sf ) const = 0 ;
      
      // Mask the value at current IP of the reconstruction obtained from 
      // the `level'-th storage of `ff' and its associated finite element.  
      virtual void mask_value_at_IP( PDE_DiscreteField const* ff,
                                     size_t level,
                                     double value,
                                     size_t ic=0 ) = 0 ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

      void print_handled_fields( std::ostream& os, 
                                 size_t indent_width ) const ;

      virtual void print_current_mesh( std::ostream& os, 
                                       size_t indent_width ) const = 0 ;

      virtual void print_local_discretization_of_fields( 
                                             std::ostream& os, 
                                             size_t indent_width ) const ;

      void print_local_discrete_variational_problem(
                                             std::ostream& os, 
                                             size_t indent_width ) const ;

      virtual void print_current_IP( std::ostream& os, 
                                     size_t indent_width ) const = 0 ;

      virtual void print_values_at_current_IP( 
                                             std::ostream& os, 
                                             size_t indent_width ) const = 0 ;

   protected: //--------------------------------------------------------

      virtual ~PDE_LocalFE( void ) ;

      PDE_LocalFE( PEL_Object* a_owner, size_t aNbSpDims ) ;

   //-- Handled fields

      size_t field_local_index( PDE_DiscreteField const* ff ) const ;
      
      int handled_field_deri( size_t iF ) const ;

      size_t handled_fields_max_nb_comps( void ) const ;

      bool jacobian_of_mapping_required( void ) const ;
      bool hessian_of_mapping_required( void ) const ;

   //-- Errors

      void raise_no_discretization_error( PDE_DiscreteField const* ff ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool go_i_th_PRE( size_t an_id_mesh ) const ;
      virtual bool go_i_th_POST( size_t an_id_mesh ) const ;

      virtual bool start_POST( void ) const ;

      virtual bool go_next_PRE( void ) const ;
      virtual bool go_next_POST( void ) const ;

      virtual bool polyhedron_PRE( void ) const ;
      virtual bool polyhedron_POST( GE_Mpolyhedron const* result ) const ;
      
      virtual bool refinement_level_PRE( void ) const ;

      virtual bool color_PRE( void ) const ;
      virtual bool color_POST( GE_Color const* result ) const ;

      virtual bool mesh_id_PRE( void ) const ;
      virtual bool mesh_id_POST( size_t result ) const ;

      virtual bool nb_local_nodes_PRE( PDE_DiscreteField const* ff ) const ;
      virtual bool nb_local_nodes_POST( size_t result,
                                        PDE_DiscreteField const* ff ) const ;

      virtual bool local_node_is_in_mesh_PRE( PDE_DiscreteField const* ff,
                                              size_t i ) const ;
      virtual bool local_node_is_in_mesh_POST( bool result,
                                               PDE_DiscreteField const* ff,
                                               size_t i ) const ;

      virtual bool local_node_location_PRE( PDE_DiscreteField const* ff,
                                            size_t i ) const ;
      virtual bool local_node_location_POST( GE_Point const* result,
                                             PDE_DiscreteField const* ff,
                                             size_t i ) const ;
      
      virtual bool global_node_PRE( PDE_DiscreteField const* ff,
                                    size_t i ) const ;
      virtual bool global_node_POST( size_t result,
                                     PDE_DiscreteField const* ff,
                                     size_t i ) const ;

      virtual bool node_refinement_level_PRE( 
                                     PDE_DiscreteField const* ff,
                                     size_t i ) const ;
      virtual bool node_refinement_level_POST( 
                                     size_t result,
                                     PDE_DiscreteField const* ff,
                                     size_t i ) const ;

      virtual bool set_calculation_point_PRE( GE_Point const* pt ) const ;
      virtual bool set_calculation_point_POST( GE_Point const* pt ) const ;

      virtual bool calculation_point_POST( GE_Point const* result ) const ;

      virtual bool value_at_pt_PRE( PDE_DiscreteField const* ff,
                                    size_t level,
                                    size_t ic ) const ;
      
      virtual bool gradient_at_pt_PRE( PDE_DiscreteField const* ff,
                                       size_t level,
                                       size_t a,
                                       size_t ic ) const ;

      virtual bool hessian_at_pt_PRE( PDE_DiscreteField const* ff, 
                                      size_t level,
                                      size_t a,
                                      size_t b,
                                      size_t ic ) const ;

      virtual bool N_at_pt_PRE( PDE_DiscreteField const* ff,
                                size_t i ) const ;

      virtual bool dN_at_pt_PRE( PDE_DiscreteField const* ff, 
                                 size_t i, 
                                 size_t a ) const ;

      virtual bool d2N_at_pt_PRE( PDE_DiscreteField const* ff,
                                  size_t i, 
                                  size_t a, 
                                  size_t b ) const ;

      virtual bool set_row_and_col_fields_PRE( 
                              PDE_DiscreteField const* row_field, 
                              PDE_DiscreteField const* col_field  ) const ;
      virtual bool set_row_and_col_fields_POST( 
                              PDE_DiscreteField const* row_field, 
                              PDE_DiscreteField const* col_field  ) const ;
      
      virtual bool row_field_node_connectivity_PRE( void ) const ;

      virtual bool row_field_node_connectivity_POST( 
                                      size_t_vector const& result ) const ;

      virtual bool col_field_node_connectivity_PRE( void ) const ;

      virtual bool col_field_node_connectivity_POST( 
                                      size_t_vector const& result ) const ;

      virtual bool nb_basis_functions_PRE( field_id sf ) const ;
      virtual bool nb_basis_functions_POST( size_t result, 
                                            field_id sf ) const ;

      virtual bool start_IP_iterator_PRE( GE_QRprovider const* qrp ) const ;
      virtual bool start_IP_iterator_POST( GE_QRprovider const* qrp ) const ;

      virtual bool go_next_IP_PRE( void ) const ;

      virtual bool valid_IP_PRE( void ) const ;

      virtual bool weight_of_IP_PRE( void ) const ;

      virtual bool coordinates_of_IP_PRE( void ) const ;
      virtual bool coordinates_of_IP_POST( GE_Point const* result ) const ;

      virtual bool value_at_IP_PRE( PDE_DiscreteField const* ff,
                                    size_t level,
                                    size_t ic=0 ) const ;

      virtual bool gradient_at_IP_PRE( PDE_DiscreteField const* ff, 
                                       size_t level,
                                       size_t a, 
                                       size_t ic=0 ) const ;

      virtual bool hessian_at_IP_PRE( PDE_DiscreteField const* ff, 
                                      size_t level,
                                      size_t a, 
                                      size_t b,
                                      size_t ic=0 ) const ;

      virtual bool N_at_IP_PRE( field_id sf, size_t i ) const ;

      virtual bool dN_at_IP_PRE( field_id sf, size_t i, size_t a ) const ;

      virtual bool d2N_at_IP_PRE( field_id sf,
                                  size_t i, size_t a, size_t b ) const ;
      
      virtual bool Ns_at_IP_PRE( field_id sf ) const ;
      virtual bool Ns_at_IP_POST( doubleVector const& result, 
                                  field_id sf ) const ;

      virtual bool dNs_at_IP_PRE( field_id sf ) const ;
      virtual bool dNs_at_IP_POST( doubleArray2D const& result, 
                                   field_id sf ) const ;

      virtual bool d2Ns_at_IP_PRE( field_id sf ) const ;
      virtual bool d2Ns_at_IP_POST( doubleArray3D const& result, 
                                    field_id sf ) const ;
      
      virtual bool mask_value_at_IP_PRE( PDE_DiscreteField const* ff,
                                         size_t level,
                                         double value,
                                         size_t ic=0 ) const ;
      virtual bool mask_value_at_IP_POST( PDE_DiscreteField const* ff,
                                          size_t level,
                                          double value,
                                          size_t ic=0 ) const ;

      virtual bool print_current_mesh_PRE( std::ostream& os, 
                                           size_t indent_width ) const ;

      virtual bool print_local_discretization_of_fields_PRE( 
                                             std::ostream& os, 
                                             size_t indent_width ) const ;

      virtual bool print_current_IP_PRE( std::ostream& os, 
                                         size_t indent_width ) const ;

      virtual bool print_values_at_current_IP_PRE( 
                                             std::ostream& os, 
                                             size_t indent_width ) const ;

      virtual bool invariant( void ) const ;

   private: //----------------------------------------------------------

      PDE_LocalFE( void ) ;
      PDE_LocalFE( PDE_LocalFE const& other ) ;
      PDE_LocalFE& operator=( PDE_LocalFE const& other ) ;

   //-- Attributes

      // NB_FIELDs    : nombre de champs (objets PDE_DiscreteField) gérées
      // FIELDS       : vecteur des champs
      // GLOB_2_iF(i) : indice (position dans FIELDS) du champ de 
      //                numero global i (si il y est, sinon PEL::bad_index())
      // F_DERIS(iF)  : ordre de derivation requis pour le champ d'indice iF
      // MAX_NB_COMPs : nombre de composantes maximum des champs de FIELDS
      size_t NB_FIELDs ;
      std::vector< size_t > GLOB_2_iF ;
      std::vector< int > F_DERIS ;
      std::vector< PDE_DiscreteField const* > FIELDS ;
      size_t MAX_NB_COMPs ;
      size_t const NB_SP_DIMS ;
      PEL_Vector* EXCLUDE_COLORS ;
      bool JACOBIAN ;
      bool HESSIAN ;

} ;

#endif

