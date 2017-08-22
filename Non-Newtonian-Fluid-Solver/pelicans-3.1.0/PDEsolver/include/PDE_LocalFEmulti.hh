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

#ifndef PDE_LOCAL_FE_MULTI_HH
#define PDE_LOCAL_FE_MULTI_HH

#include <PDE_LocalFE.hh>

#include <boolVector.hh>
#include <boolArray2D.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>
#include <doubleVector.hh>
#include <intArray2D.hh>
#include <intVector.hh>
#include <size_t_array2D.hh>
#include <size_t_array3D.hh>
#include <size_t_vector.hh>

class PDE_BFvalues ;
class PDE_BasisFunction ;
class PDE_MeshFE ;

class PEL_EXPORT PDE_LocalFEmulti : public PDE_LocalFE
{
   public: //-----------------------------------------------------------

   //-- Local discretization of fields

      virtual size_t nb_local_nodes( PDE_DiscreteField const* ff ) const ;

      virtual GE_Point const* local_node_location( PDE_DiscreteField const* ff,
                                                   size_t i ) const ;

      virtual bool local_node_is_in_mesh( PDE_DiscreteField const* ff,
                                          size_t i ) const ;

      virtual size_t global_node( PDE_DiscreteField const* ff,
                                  size_t i ) const ;

      virtual size_t node_refinement_level( PDE_DiscreteField const* ff,
                                            size_t i ) const ;

   //-- Computations at a given point of the current mesh(101.5)

      virtual GE_Point const* calculation_point( void ) const ;

      virtual double value_at_pt( PDE_DiscreteField const* ff,
                                  size_t level,
                                  size_t ic=0 ) const ;

      virtual double gradient_at_pt( PDE_DiscreteField const* ff,
                                     size_t level,
                                     size_t a,
                                     size_t ic=0 ) const ;

      virtual double hessian_at_pt( PDE_DiscreteField const* ff, 
                                    size_t level,
                                    size_t a,
                                    size_t b,
                                    size_t ic=0 ) const ;

      virtual double N_at_pt( PDE_DiscreteField const* ff, 
                              size_t i ) const ;

      virtual double dN_at_pt( PDE_DiscreteField const* ff, 
                               size_t i, 
                               size_t a ) const ;

      virtual double d2N_at_pt( PDE_DiscreteField const* ff,
                                size_t i, 
                                size_t a, 
                                size_t b ) const ;

   //-- Local discrete variational problem definition(102)

      virtual void set_row_and_col_fields( PDE_DiscreteField const* row_field, 
                                   PDE_DiscreteField const* col_field ) ;

      virtual PDE_DiscreteField const* field( field_id sf ) const ;

      virtual size_t nb_basis_functions( field_id sf ) const ;

      virtual size_t_vector const& row_field_node_connectivity( void ) const ;

      virtual size_t_vector const& col_field_node_connectivity( void ) const ;

   //-- IP-iterator movement (IP standing for Integration Points)(103)

      virtual void start_IP_iterator( GE_QRprovider const* qrp ) ;

      virtual void go_next_IP( void ) ;

      virtual bool valid_IP( void ) const ;

      virtual double weight_of_IP( void ) const ;

      virtual GE_Point const* coordinates_of_IP( void ) const ;

      virtual double value_at_IP( PDE_DiscreteField const* ff,
                                  size_t level,
                                  size_t ic=0 ) const ;

      virtual double gradient_at_IP( PDE_DiscreteField const* ff, 
                                     size_t level,
                                     size_t a, 
                                     size_t ic=0 ) const ;

      virtual double hessian_at_IP( PDE_DiscreteField const* ff, 
                                    size_t level,
                                    size_t a,
                                    size_t b,
                                    size_t ic=0 ) const ;

      virtual double N_at_IP( field_id sf, size_t i ) const ;

      virtual double dN_at_IP( field_id sf, size_t i, size_t a ) const ;

      virtual double d2N_at_IP( field_id sf, 
                                size_t i, size_t a, size_t b ) const ;

      virtual doubleVector const& Ns_at_IP( field_id sf ) const ;
      
      virtual doubleArray2D const& dNs_at_IP( field_id sf ) const ;
      
      virtual doubleArray3D const& d2Ns_at_IP( field_id sf ) const ;
      
      virtual void mask_value_at_IP( PDE_DiscreteField const* ff,
                                     size_t level,
                                     double value,
                                     size_t ic=0 ) ;

   //-- Input - Output

      virtual void print_local_discretization_of_fields( 
                                             std::ostream& os, 
                                             size_t indent_width ) const ;

      virtual void print_current_IP( std::ostream& os, 
                                     size_t indent_width ) const ;

      virtual void print_values_at_current_IP( std::ostream& os, 
                                               size_t indent_width ) const ;

   protected: //--------------------------------------------------------

      virtual ~PDE_LocalFEmulti( void ) ;

      PDE_LocalFEmulti( PEL_Object* a_owner, size_t aNbSpDims ) ;

   //-- Internals on the current mesh
      
      void clear_local_fields( void ) ;

      void append_local_field( size_t iF, PDE_MeshFE const* lf ) ;

      void terminate_local_fields( void ) ;

      virtual bool is_in_mesh( PDE_DiscreteField const* ff,
                               size_t i,
                               PDE_BasisFunction const* bf ) const ;

      size_t nb_nested_meshes( PDE_MeshFE const* mesh ) const ;

      void reset_row_and_col_fields( void ) ;

      void reset_calculation_point( void ) ;

      void set_CP( GE_Point const* pt ) ;

      void connect_CP( PDE_MeshFE const* mesh,
                       size_t coarse_level,
                       GE_Point const* pt_ref,
                       GE_Matrix const* tr_jac,
                       doubleArray3D const* hessian ) ;

      void reset_IP_iterator( void ) ;

      void set_single_fictive_IP( GE_Point const* pt ) ;
      
      void set_nb_IPs( size_t aNbItgPts ) ;

      void append_IP( GE_Point const* pt, double weight ) ;

      void connect_IP( PDE_MeshFE const* mesh,
                       size_t coarse_level,
                       GE_Point const* pt_ref,
                       GE_Matrix const* tr_jac,
                       doubleArray3D const* hessian ) ;

      virtual GE_QuadratureRule const* quadrature_rule( 
                                          GE_QRprovider const* qrp ) const ;

      virtual void compute_itg_pts( GE_QuadratureRule const* qr ) = 0 ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool quadrature_rule_PRE( GE_QRprovider const* qrp ) const ;
      virtual bool quadrature_rule_POST( GE_QuadratureRule const* result,
                                         GE_QRprovider const* qrp ) const ;

   private: //----------------------------------------------------------

      PDE_LocalFEmulti( void ) ;
      PDE_LocalFEmulti( PDE_LocalFEmulti const& other ) ;
      PDE_LocalFEmulti& operator=( PDE_LocalFEmulti const& other ) ;

      void locate_discretization_at_pt( size_t iF,
                                        size_t& im, size_t& ee ) const ;

      PDE_BFvalues* C_BFvalues( size_t im, size_t cl, size_t ee ) const ;

      void locate_discretization_at_IP( size_t iF, size_t ip, 
                                        size_t& im, size_t& ee ) const ;

      PDE_BFvalues* I_BFvalues( size_t im, 
                                size_t cl, size_t ee, size_t ip ) const ;

      size_t local_mesh_index( PDE_MeshFE const* mesh ) const ;

      void build_connectivities( void ) ;

      void prepare_for_Ns_requests( void ) ;

      void build_Ns( size_t iF, size_t im, size_t ee,
                     doubleVector& val,
                     doubleArray2D& grad,
                     doubleArray3D& grad_grad ) ;
      
      static void initialize_IP_vector( size_t new_length,
                                        size_t nb_dims,
                                        PEL_Vector* vec,
                                        PEL_Vector* pool ) ;

      static void initialize_ELM_BF_vectors_1( size_t new_length,
                                               PEL_Vector* vec2,
                                               PEL_Vector* vec3 ) ;

      void initialize_BF_vector_2( size_t nb_pts,
                                   PEL_Vector* vec,
                                   PEL_Vector* pool ) ;

   //-- Attributes

      /* indices :
            iF<nb_handled_fields()  im<NB_MESHES e<NB_ELMS(im)  
            ip<NB_IPs  il<NB_LOCAL_NODES(iF)

      NB_MESHES   
         number of PDE_MeshFE instances involved in the local discrete 
         representation of the handled field (all the meshes intersecting 
         the current mesh)
      MESHES(im)  
         im-th observed PDE_MeshFE instance
      LOCAL_FIELDS 
         vector of all PDE_LocalField instances
      LF_index(im,iF) 
         for the iF-th handled field, index in LOCAL_FIELDS of the 
         PDE_LocalField instance associated to the im-th mesh

      NB_ELMS(im)
         number of PDE_ReferenceElement objects of the im-th mesh
      ELMS(im) : PEL_Vector
         ELMS(im)(e) e-th PDE_ReferenceElement object of the im-th mesh
      ELMS_DERIS(im,e)
         calculation requirement for the e-th PDE_ReferenceElement object
         of the im-th mesh
      ELM_index(im,iF)
         index in ELMS(im) of the PDE_ReferenceElement object associated to 
         the iF-th handled field on the im-th mesh (PEL::bad_index() if none)
         ELMS(im)(ELM_index(im,iF)) : PDE_ReferenceElement object of the
         iF-th handled field in the im-th mesh
      NB_LOCAL_NODES(iF)
         number of nodes involved in the discrete local representation of 
         the iF-th handled field (also called local nodes)
      GLOBAL_NODES(iF,il)
         global number of the il-th local node of the iF-th handled field
      NODE_LOCATIONS(iF,il)
         location (GE_Point) of the il-th local node of the iF-th handled 
         field
      ELM_BF_index(im,iF,il)
         index in PDE_ReferenceElement of the il-th node involved in the 
         discrete representation on the im-th mesh of the iF-th handled field 
         (PEL::bad_index() if that il-th node is not in the im-th mesh)

      CP_LOCATION
         location of the calculation point
      CP_NB_MESHES
         number of meshes that contain the calculation point
      CP_MESHES(ii)
         index in MESHES of the ii-th mesh that contains the calculation point
      C_BF_VALS(im) : PEL_Vector*
         C_BF_VALS(im)(e) : PDE_BFvalues object of reference element e
         at the calculation point in the im-th mesh

      NB_IPs
         number of itg points
      IP_LOCATIONS(ip)
         location (GE_Point*) of the ip-th itg points
      IP_WEIGHTS(ip)
         weight of the ip-th itg point
      IP_NB_MESHES(ip)
         number of meshes that contain the ip-th itg points
      IP_MESHES(ip,ii)
         index in MESHES of the ii-th mesh that contains the ip-th itg point
         ii<IP_MESHES(ip)
      BF_VALS(im) : PEL_Vector*
         BF_VALS(im)(e+NB_ELMs(im)*ip) : PDE_BFvalues object of reference
            element e at the ip-th itg point in the im-th mesh
      */

      size_t NB_MESHES ;
      PEL_Vector* MESHES ;
      size_t_vector NB_NESTED_LEVELS ;
      size_t MAX_NB_NODES ;
      boolArray2D HAS_DISC ;

      intArray2D ELM_DERIS ;
      size_t_array2D ELM_index ;
      size_t_vector NB_LOCAL_NODES ;
      size_t_array2D GLOBAL_NODE ;
      PEL_Vector* BASIS_FUNC ;
      PEL_Vector* NODE_LOCATIONS ;
      mutable boolVector NODE_LOC_OK ;
      size_t_array3D ELM_BF_index ;
      size_t_array3D REFINEMENT_LEVEL ;

      bool HAS_CP ;
      GE_Point* CP_LOCATION ;
      size_t CP_NB_MESHES ;
      size_t_vector CP_MESHES ;
      PEL_Vector* C_BF_VALS ;
      PEL_Vector* C_BF_VALS_POOL ;

      GE_QRprovider const* QR_PROVIDER ;
      size_t QR_PROVIDER_STAT ;

      size_t NB_IPs ;
      PEL_Vector* IP_LOCATIONS ;
      PEL_Vector* IP_LOCATIONS_POOL ;
      doubleVector IP_WEIGHTS ;
      size_t_vector IP_NB_MESHES ;
      size_t_array2D IP_MESHES ;
      PEL_Vector* I_BF_VALS ;
      PEL_Vector* I_BF_VALS_POOL ;
      size_t i_IP ;

      doubleArray3D MASKED_VALUES ;
      size_t MASKED_LEVEL ;

      size_t ROW_iF ;
      size_t COL_iF ;
      size_t ROW_im ;
      size_t COL_im ;
      size_t ROW_ee ;
      size_t COL_ee ;
      size_t_vector ROW_NODES ;
      size_t_vector COL_NODES ;
      doubleVector ROW_Ns ;
      doubleVector COL_Ns ;
      doubleArray2D ROW_dNs ;
      doubleArray2D COL_dNs ;
      doubleArray3D ROW_d2Ns ;
      doubleArray3D COL_d2Ns ;
} ;

#endif

