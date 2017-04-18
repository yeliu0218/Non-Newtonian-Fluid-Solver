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

#ifndef PDE_LOCAL_FE_SINGLE_HH
#define PDE_LOCAL_FE_SINGLE_HH

#include <PDE_LocalFE.hh>

#include <boolVector.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <size_t_array2D.hh>
#include <size_t_vector.hh>

class PEL_Vector ;

class PDE_BasisFunction ;
class PDE_BFvalues ;
class PDE_DiscreteField ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalEquation ;
class PDE_MeshFE ;
class PDE_ReferenceElement ;

class GE_Color ;
class GE_Point ;
class GE_Matrix ;
class GE_Mpolyhedron ;
class GE_QuadratureRule ;
class GE_QRprovider ;

class PEL_EXPORT PDE_LocalFEsingle : public PDE_LocalFE
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

   //-- Computations at a given point of the current mesh

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

   //-- Local discrete variational problem definition

      virtual void set_row_and_col_fields( 
                                   PDE_DiscreteField const* row_field, 
                                   PDE_DiscreteField const* col_field  ) ;

      virtual PDE_DiscreteField const* field( field_id sf ) const ;

      virtual size_t nb_basis_functions( field_id sf ) const ;

      virtual size_t_vector const& row_field_node_connectivity( void ) const ;

      virtual size_t_vector const& col_field_node_connectivity( void ) const ;

   //-- IP-iterator movement (IP standing for Integration Points)

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

      virtual void print_values_at_current_IP( 
                                             std::ostream& os, 
                                             size_t indent_width ) const ;

   protected: //--------------------------------------------------------

      virtual ~PDE_LocalFEsingle( void ) ;

      PDE_LocalFEsingle( PEL_Object* a_owner, size_t aNbSpDims ) ;
      
   //-- Internals on the current mesh
      
      void set_single_mesh( PDE_MeshFE* a_mesh ) ;

      virtual bool is_in_mesh( PDE_BasisFunction const* bf ) const = 0 ;

      size_t nb_nested_meshes( void ) const ;

      void reset_row_and_col_fields( void ) ;

      void reset_calculation_point( void ) ;

      void set_CP( GE_Point const* pt ) ;
             
      void connect_CP( size_t coarse_level, 
                       GE_Point const* pt_ref,
                       GE_Matrix const* tr_jac,
                       doubleArray3D const* hessian ) ;

      void reset_IP_iterator( void ) ;

      void set_nb_IPs( size_t a_nb_ips ) ;

      void append_IP( double weight ) ;

      void connect_IP( size_t coarse_level,
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

      PDE_LocalFEsingle( void ) ;
      PDE_LocalFEsingle( PDE_LocalFEsingle const& other ) ;
      PDE_LocalFEsingle const& operator=( PDE_LocalFEsingle const& other ) ;

      void init_C_BFvalues( size_t nb_cl, size_t nb_ee ) ;

      PDE_BFvalues* C_BFvalues( size_t cl, size_t ee ) const ;
      
      void init_I_BFvalues( size_t nb_cl, size_t nb_ee, size_t nb_ip ) ;

      PDE_BFvalues* I_BFvalues( size_t cl, size_t ee, size_t ip ) const ;

      void build_connectivities( void ) ;

      void prepare_for_Ns_requests( void ) ;

      void build_Ns( size_t iF, doubleVector& val, doubleArray2D& grad, 
                     doubleArray3D& grad_grad ) ;

   //-- Attributes

      /* indices
            iF<nb_handled_fields()  e<NB_ELMS(im)  
     
      LOCAL_FIELDS 
         vector of all PDE_LocalField instances
      NB_ELMS
         number of PDE_ReferenceElement objects
      ELMS(e)
         e-th PDE_ReferenceElement object
      ELM_DERIS(e)
         calculation requirement for the e-th PDE_ReferenceElement object
      ELM_index(iF)
         index in ELMS of the PDE_ReferenceElement object of the iF-th field
         ELMS( ELM_index(iF) ) : PDE_ReferenceElement* of field(iF)

      C_BF_VALS(e) 
         PDE_BFvalues object of reference element ELMS(e) at the calcultion 
         point
      
      I_BF_VALS(e+NB_ELMS*ip)
         PDE_BFvalues object of reference element ELMS(e) at the ip-th 
         itg point
      */

      size_t NB_NESTED_LEVELS ;
      PDE_MeshFE* THE_MESH ;

      size_t MAX_NB_FIELDS ;
      size_t MAX_NB_NODES ;

      intVector ELM_DERIS ;
      size_t_vector ELM_index ;
      size_t_vector NB_LOCAL_NODES ;
      size_t_array2D GLOBAL_NODE ;
      PEL_Vector* BASIS_FUNC ;
      size_t_array2D ELM_BF_index ;
      PEL_Vector* NODE_LOCATIONS ;
      mutable boolVector NODE_LOC_OK ;
      size_t_array2D REFINEMENT_LEVEL ;

      bool HAS_CP ;
      GE_Point* CP_LOCATION ;
      PEL_Vector* C_BF_VALS ;

      GE_QRprovider const* QR_PROVIDER ;
      size_t QR_PROVIDER_STAT ;
      GE_QuadratureRule const* QR ;

      size_t NB_IPs ;
      GE_Point* IP_LOCATION ;
      doubleVector IP_WEIGHTS ;
      PEL_Vector* I_BF_VALS ;
      size_t i_IP ;

      doubleArray3D MASKED_VALUES ;
      size_t MASKED_LEVEL ;

      size_t_vector SF_2_iF ;
      size_t_vector ROW_NODES ;
      size_t_vector COL_NODES ;
      PDE_BFvalues const* ROW_bfv ;
      PDE_BFvalues const* COL_bfv ;
      doubleVector ROW_Ns ;
      doubleVector COL_Ns ;
      doubleArray2D ROW_dNs ;
      doubleArray2D COL_dNs ;
      doubleArray3D ROW_d2Ns ;
      doubleArray3D COL_d2Ns ;

      GE_Point* PT ;
} ;

#endif

