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

#ifndef PDE_BF_VALUES_HH
#define PDE_BF_VALUES_HH

#include <PEL_Object.hh>

#include <doubleVector.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>

class GE_Matrix ;
class GE_Point ;

class PDE_DiscreteField ;
class PDE_MeshFE ;
class PDE_ReferenceElement ;

// Basis function evaluation for a given reference element at a
// given point in reference element.
// First and second derivative are given with respect to global coordinates.

class PEL_EXPORT PDE_BFvalues : public PEL_Object
{
   public: //-----------------------------------------------------------

      // indicator for zero-th order of spatial derivation
      static int const N ;

      // indicator for first order of spatial derivation
      static int const dN ;

      // indicator for second order of spatial derivation
      static int const d2N ;

  //-- Instance delivery and initialization

      static PDE_BFvalues* create( PEL_Object* a_owner ) ;
      
      // New PDE_BFvalues object evaluating values of basis
      // functions of reference element `r' with a given `tr_jac' transformation
      // matrix from global coordinate to reference one at reference point
      // `pt_ref'.
      // `order' value is used to indicate if value (N), first derivative (dN)
      // and second one (d2N) has to be computed.
      // When only basis function value are requested, matrix can be omitted.
      static PDE_BFvalues* create( PEL_Object* a_owner,
                                   PDE_ReferenceElement const* r,
                                   int order,
                                   GE_Point const* pt_ref,
                                   GE_Matrix const* tr_jac,
                                   doubleArray3D const* hess ) ;
      
      // Re-initialize self with current geometrical data.
      void re_initialize( PDE_ReferenceElement const* r,
                          int order,
                          GE_Point const* pt_ref,
                          GE_Matrix const* tr_jac,
                          doubleArray3D const* hess ) ;

   //-- Values retrieving

      // Basis function `i' value.
      double N_at_pt( size_t i ) const ;
      
      // Basis function values.
      doubleVector const& Ns_at_pt( void ) const ;
      
      // Basis function `i' first derivative with respect
      //  to `a' direction.
      double  dN_at_pt( size_t i, size_t a ) const ;
      
      // Basis function gradient.
      doubleArray2D const& dNs_at_pt( void ) const ;
      
      // Basis function second derivative for `i' basis function
      // with respect to `a' and `b' directions.
      double d2N_at_pt( size_t i, size_t a, size_t b ) const ;

      // Basis function second derivative.
      doubleArray3D const& d2Ns_at_pt( void ) const ;

      //-- Structure
      PDE_ReferenceElement const* reference_element( void ) const ;
      size_t nb_basis_function( void ) const ;
      size_t space_dimension( void ) const ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_BFvalues( void ) ;
     ~PDE_BFvalues( void ) ;
      PDE_BFvalues( PDE_BFvalues const& other ) ;
      PDE_BFvalues& operator=( PDE_BFvalues const& other ) ;

      PDE_BFvalues( PEL_Object* a_owner ) ;

      void   computeN( GE_Point const* pt_ref ) ;
      void  computedN( GE_Matrix const* tr_jac ) ;      
      void computed2N( GE_Matrix const* tr_jac,
                       doubleArray3D const* hess ) ;
      void  computedNref( GE_Point const* pt_ref ) ;      
      void computed2Nref( GE_Point const* pt_ref ) ;
      bool same_point( GE_Point const* pt_ref ) const ;
      
   //-- Attributes

      PDE_ReferenceElement const* ELM ;
      int ELM_DERI ;
      size_t DIM ;
      GE_Point* PT_REF ;
      
      // number of basis functions
      size_t NB_BFs ;

      // BFs( i ) : i-th basis function
      doubleVector BFs ;

      // d_BFs( i, a ) : a-th derivative of the i-th ref BF
      doubleArray2D d_BFs_REF ;

      // d_BFs( i, a ) : a-th derivative of the i-th real BF
      doubleArray2D d_BFs ;

      // d2_BFs_REF_vec( i, a, b ) : (a,b)-th 2nd derivative of the ref BF
      doubleArray3D d2_BFs_REF ;

      // d2_BFs( i, a, b ) : (a,b)-th 2nd derivative of the real BF
      doubleArray3D d2_BFs ;
} ;

#ifndef OUTLINE
   #include <PDE_BFvalues.icc>
#endif

#endif
