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

#include <PDE_BFvalues.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_MeshFE.hh>
#include <PDE_ReferenceElement.hh>

#include <GE_Point.hh>
#include <GE_Matrix.hh>
#include <GE_Mpolyhedron.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Root.hh>

#include <LA_DenseMatrix.hh>
#include <LA_GaussLU_DS.hh>
#include <LA_SeqVector.hh>

#ifdef OUTLINE
   #define inline
   #include <PDE_BFvalues.icc>
   #undef inline
#endif

int const PDE_BFvalues::N   = 1 ;
int const PDE_BFvalues::dN  = 2 ;
int const PDE_BFvalues::d2N = 4 ;

//----------------------------------------------------------------------
PDE_BFvalues* 
PDE_BFvalues:: create( PEL_Object* a_owner )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BFvalues:: create( PEL_Object* )" ) ;

   PDE_BFvalues* result = new PDE_BFvalues( a_owner ) ;

   PEL_CHECK_POST( result!=0 && result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_BFvalues* 
PDE_BFvalues:: create( PEL_Object* a_owner,
                       PDE_ReferenceElement const* r,
                       int order,
                       GE_Point const* pt_ref,
                       GE_Matrix const* tr_jac,
                       doubleArray3D const* hess )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BFvalues:: create" ) ;
   PEL_CHECK_PRE( r!=0 ) ;

   PDE_BFvalues* result = new PDE_BFvalues( a_owner ) ;
   result->re_initialize( r, order, pt_ref, tr_jac, hess ) ;  

   PEL_CHECK_POST( result!=0 && result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->reference_element()==r ) ;
   PEL_CHECK_POST( result->nb_basis_function()==r->nb_nodes() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_BFvalues:: PDE_BFvalues( PEL_Object* a_owner )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , ELM( 0 )
   , ELM_DERI( 0 )
   , DIM( PEL::bad_index() )
   , PT_REF( 0 )
   , NB_BFs( PEL::bad_index() )
   , BFs( 0 )
   , d_BFs_REF( 0, 0 )
   , d_BFs( 0, 0 )
   , d2_BFs_REF( 0, 0, 0 )
   , d2_BFs( 0, 0, 0 )
{
}

//----------------------------------------------------------------------
PDE_BFvalues:: ~PDE_BFvalues( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PDE_BFvalues:: re_initialize( PDE_ReferenceElement const* r,
                              int order,
                              GE_Point const* pt_ref,
                              GE_Matrix const* tr_jac,
                              doubleArray3D const* hess )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BFvalues:: re_initialize(r,order,...)" ) ;
   PEL_CHECK_PRE( r!=0 ) ;
   PEL_CHECK_PRE( pt_ref!=0 ) ;
   PEL_CHECK_PRE( IMPLIES( order != N, 
                           tr_jac!=0 && tr_jac->determinant()!=0 ) ) ;
   
   bool build_ref_vals = ( r!=ELM || ELM_DERI!=order || !same_point(pt_ref) ) ;
   ELM = r ;
   NB_BFs = r->nb_nodes() ;
   DIM = r->dimension() ; //??????????????? FAUX
   ELM_DERI = order ;
   if( build_ref_vals )
   {
      if( order & N )
      {
         computeN( pt_ref ) ;
      }
      if( order & dN )
      {
         computedNref( pt_ref ) ;
      }
      if( order & d2N )
      {
         computed2Nref( pt_ref ) ;
      }
   }

   if( order!=N)
   {
      if( ELM_DERI & dN )
      {
         computedN( tr_jac ) ;
      }
      if( ELM_DERI & d2N )
      {
         computed2N( tr_jac, hess ) ;
      }
   }
   if( PT_REF==0 )
   {
      PT_REF = pt_ref->create_clone( this ) ;
   }
   else
   {
      PT_REF->copy( pt_ref ) ;
   }
}

//----------------------------------------------------------------------
void
PDE_BFvalues:: computeN( GE_Point const* pt_ref )
//----------------------------------------------------------------------
{
   if( BFs.size()!=NB_BFs )
   {
      BFs.re_initialize( NB_BFs ) ;
   }
   for( size_t i=0 ; i<NB_BFs ; i++ )
   {
      BFs( i ) = ELM->N_local( i, pt_ref ) ;
   }
}

//----------------------------------------------------------------------
void
PDE_BFvalues:: computedN( GE_Matrix const* tr_jac )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BFvalues:: computedN" ) ;
   PEL_CHECK( ELM_DERI & dN ) ;
   
   tr_jac->invert( d_BFs_REF, d_BFs ) ;
}

//----------------------------------------------------------------------
void
PDE_BFvalues:: computed2N( GE_Matrix const* tr_jac,
                           doubleArray3D const* hess )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BFvalues:: computed2N" ) ;
   PEL_CHECK( ELM_DERI & dN ) ;
   PEL_CHECK( ELM_DERI & d2N ) ;

   static LA_GaussLU_DS* solver =
                    LA_GaussLU_DS::create( PEL_Root::object(), 1.E-30 ) ;
   static LA_DenseMatrix* JAC_JAC =
                   LA_DenseMatrix::create( PEL_Root::object(), DIM, DIM ) ;
   static LA_SeqVector* vec_x = LA_SeqVector::create( PEL_Root::object(), DIM ) ;
   static LA_SeqVector* vec_b = LA_SeqVector::create( PEL_Root::object(), DIM ) ;
   if( vec_b->nb_rows()!=DIM*DIM )
   {
      vec_b->re_initialize( DIM*DIM ) ;
      vec_x->re_initialize( DIM*DIM ) ;
      JAC_JAC->re_initialize( DIM*DIM, DIM*DIM ) ;
   }

   // Set JAC_JAC :
   for( size_t a = 0 ; a<DIM ; ++a )
   {
      for( size_t b = 0 ; b<DIM ; ++b )
      {
         size_t ii = a*DIM + b ;
         for( size_t l = 0 ; l<DIM ; ++l )
         {
            for( size_t m = 0 ; m<DIM ; ++m )
            {
               size_t jj = l*DIM + m ;
               double xx = tr_jac->item( a, l ) * tr_jac->item( b, m ) ;
               JAC_JAC->set_item( ii, jj, xx ) ;
            }
         }
      }
   }

   // LU factorization of JAC_JAC :
   JAC_JAC->synchronize() ;
   solver->set_matrix( JAC_JAC ) ;

   for( size_t i=0 ; i<NB_BFs ; i++ )
   {
      // Compute RHS :
      for( size_t a=0 ; a<DIM ; a++ )
      {
         for( size_t b=0 ; b<DIM ; b++ )
         {
            double xx = d2_BFs_REF( i, a, b ) ;
            if( hess!=0 )
            {
               for( size_t d=0 ; d<DIM ; ++d )
               {
                  xx -= (*hess)( d, a, b )*d_BFs( i, d ) ;
               }
            }
            vec_b->set_item( a*DIM+b, xx ) ;
         }
      }

      // Solve :
      vec_b->synchronize() ;
      solver->solve( vec_b, vec_x ) ;

      // Set d2_BFs :
      for( size_t a=0 ; a<DIM ; a++ )
      {
         for( size_t b=0 ; b<DIM ; b++ )
         {
            d2_BFs( i, a , b ) = vec_x->item( a*DIM+b ) ;
         }
      }
   }
   solver->unset_matrix() ;
}

//----------------------------------------------------------------------
void
PDE_BFvalues:: computedNref( GE_Point const* pt_ref )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BFvalues:: computedNref" ) ;
   if( d_BFs.index_bound(0)!=NB_BFs || d_BFs.index_bound(1)!=DIM )
   {
      d_BFs.re_initialize( NB_BFs, DIM ) ;
      d_BFs_REF.re_initialize( NB_BFs, DIM ) ;
   }
   
   for( size_t a=0 ; a<DIM ; a++ )
   {
      for( size_t i=0 ; i<NB_BFs ; i++ )
      {
         d_BFs_REF( i, a ) = ELM->dN_local( i, a, pt_ref ) ;
      }
   }
}

//----------------------------------------------------------------------
void
PDE_BFvalues:: computed2Nref( GE_Point const* pt_ref )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BFvalues:: computed2Nref" ) ;

   if( d2_BFs.index_bound(0)!=NB_BFs || d2_BFs.index_bound(1)!=DIM )
   {
      d2_BFs.re_initialize( NB_BFs, DIM, DIM ) ;
      d2_BFs_REF.re_initialize( NB_BFs, DIM, DIM ) ;
   }
   
   for( size_t i=0 ; i<NB_BFs ; i++ )
   {
      for( size_t a=0 ; a<DIM ; a++ )
      {
         for( size_t b=0 ; b<DIM ; b++ )
         {
            d2_BFs_REF( i, a, b ) = ELM->d2N_local( i, a, b, pt_ref ) ;
	 }
      }
   }   
}

//----------------------------------------------------------------------
PDE_ReferenceElement const* 
PDE_BFvalues:: reference_element( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BFvalues:: reference_element" ) ;

   PDE_ReferenceElement const*  result = ELM ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->nb_nodes()==nb_basis_function() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_BFvalues:: same_point( GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BFvalues:: same_point" ) ;

   bool result = PT_REF!=0 ;
   for( size_t i=0 ; result && i<DIM ; i++ )
   {
      result = result && PT_REF->coordinate(i)==pt_ref->coordinate(i) ;
   }
   return result ;
}


