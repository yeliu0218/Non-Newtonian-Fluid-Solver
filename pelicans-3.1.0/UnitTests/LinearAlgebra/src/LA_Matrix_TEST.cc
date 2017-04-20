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

#include <LA_Matrix_TEST.hh>

#include <LA_SeqVector.hh>
#include <LA_DenseMatrix.hh>
#include <LA_Matrix.hh>
#include <LA_PelMatrix.hh>
#include <LA_Scatter.hh>
#include <LA_SeqMatrix.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_System.hh>


#include <string>
#include <iostream>
#include <iomanip>

#include <PEL_Exceptions.hh>

#define TEST_FAILED( X, Y ) {}

/* #define TEST_FAILED( X, Y ) \
   try \
   { \
      PEL::out() << "Je vais faire ca " << #X << std::endl ; \
      X ; \
      PEL::out() << "J'ai fais ca " << #X << std::endl ; \
      notify_one_test_result( Y , false ) ; \
   } \
   catch( PEL_Exceptions::Error ) \
   { \
      notify_one_test_result( Y , true ) ; \
   }
*/


LA_Matrix_TEST const*
LA_Matrix_TEST::PROTOTYPE = new LA_Matrix_TEST()  ;

//-------------------------------------------------------------------------
LA_Matrix_TEST:: LA_Matrix_TEST(void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "LA_Matrix", "LA_Matrix_TEST" )
   , MY_EPS( PEL::bad_double() )
   , MY_MIN( PEL::bad_double() )
   , TYPE( "" )
   , MAT_PROTO( 0 )
{
}

//-------------------------------------------------------------------------
LA_Matrix_TEST:: LA_Matrix_TEST( std::string const& matrixClass,
                                 std::string const& registration )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( matrixClass, registration )
   , MY_EPS( 1.E-8 )
   , MY_MIN( 1.E-10 )
   , TYPE( "" )
   , MAT_PROTO( 0 )
{
}

//-------------------------------------------------------------------------
LA_Matrix_TEST:: ~LA_Matrix_TEST( void )
//-------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_Matrix_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix_TEST:: process_one_test" ) ;

   PEL_ModuleExplorer* proto_exp = exp->create_subexplorer( 0, "PROTOTYPE" ) ;
   MAT_PROTO = LA_Matrix::make( this, proto_exp ) ;
   TYPE = proto_exp->string_data("concrete_name") ;
   proto_exp->destroy() ;
   notify_one_test_result(  TYPE+" make", MAT_PROTO!=0 ) ;

   if( exp->has_entry( "dbl_epsilon" ) )
   {
      MY_EPS = exp->double_data( "dbl_epsilon" ) ;
      exp->test_data( "dbl_epsilon", "dbl_epsilon>0." ) ;
      exp->set_default( "dbl_epsilon", "1.E-8" ) ;
   }
   if( exp->has_entry( "dbl_minimum" ) )
   {
      MY_MIN = exp->double_data( "dbl_minimum" ) ;
      exp->test_data( "dbl_minimum", "dbl_minimum>0." ) ;
      exp->set_default( "dbl_minimum", "1.E-10" ) ;
   }


   std::string const& test_matrix = exp->string_data( "test_matrix" ) ;
   MAT_PROTO->readMM( test_matrix ) ;
   notify_one_test_result( TYPE+ " readMM synchronized",
                           MAT_PROTO->is_synchronized() ) ;
   test( exp, MAT_PROTO, test_matrix ) ;

   destroy_possession( MAT_PROTO ) ; MAT_PROTO = 0 ;
   TYPE = "" ;
   MY_EPS = 1.E-8 ;
   MY_MIN = 1.E-10 ;
}

//-------------------------------------------------------------------------
bool
LA_Matrix_TEST:: matrix_equality( LA_Matrix const* mat1,
                                  LA_SeqMatrix const* mat2 ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix_TEST:: matrix_equality" ) ;
   PEL_CHECK_PRE( mat1 != 0 ) ;
   PEL_CHECK_PRE( mat1->is_synchronized() ) ;
   PEL_CHECK_PRE( mat2 != 0 ) ;

   bool ok = mat1->nb_rows() == mat2->nb_rows() &&
             mat1->nb_cols() == mat2->nb_cols() ;
   if( ok )
   {
      LA_SeqMatrix const* loc1 = mat1->create_local_matrix(0) ;
      LA_MatrixIterator* it_mat = mat1->create_stored_item_iterator( 0 ) ;
      ok = loc1->nb_rows() == mat2->nb_rows() &&
           loc1->nb_cols() == mat2->nb_cols() ;
      if( ok )
      {
         for( it_mat->start_all_items() ; it_mat->is_valid() ; it_mat->go_next() )
         {
            size_t i = it_mat->row() ;
            size_t j = it_mat->col() ;
            if( !mat1->is_symmetric() || i<=j )
            {
               bool okt = double_equality( it_mat->item(), mat2->item(i,j) ) ;
               if( !okt )
               {
                  out() << "    matrix_equality i " << i << " j " << j << std::endl ;
               }
               ok = ok && okt ;
            }
         }
      }
      loc1->destroy() ; it_mat->destroy() ;
   }
   return ok ;
}

//-------------------------------------------------------------------------
bool
LA_Matrix_TEST:: vector_equality( LA_SeqVector const* vec1,
                                  LA_SeqVector const* vec2 ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix_TEST:: vector_equality" ) ;
   PEL_CHECK_PRE( vec1 != 0 ) ;
   PEL_CHECK_PRE( vec2 != 0 ) ;

   bool ok = vec1->nb_rows() == vec2->nb_rows() ;
   if( ok )
   {
      for( size_t i=0 ; i<vec1->nb_rows() ; i++ )
      {
         bool okt = double_equality( vec1->item(i), vec2->item(i) ) ;
         if( !ok )
         {
            out() << " vector_equality diff " << i << std::endl ;
         }
         ok = ok && okt ;
      }
   }
   return ok ;
}

//-------------------------------------------------------------------------
bool
LA_Matrix_TEST:: double_equality( double const x, double const y ) const
//-------------------------------------------------------------------------
{
   bool result = PEL::double_equality( x, y, MY_EPS, MY_MIN ) ;
   if( !result )
   {
      out() << std::setprecision( 10 )
            << std::setiosflags( std::ios::scientific )
            << "failed : " << x << " <-> " << y << std::endl ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
void
LA_Matrix_TEST:: test( PEL_ModuleExplorer const* exp,
                       LA_Matrix const* matrix,
                       std::string const& test_matrix )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix_TEST:: test" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   PEL_CHECK_PRE( matrix != 0 ) ;

   size_t const N = 100 ;

   std::string label = TYPE+" "+ PEL_System::basename( test_matrix ) ;

   LA_Matrix* clone = matrix->create_matrix( 0 ) ;
   clone->set( matrix ) ;
   clone->synchronize() ;

   LA_Matrix* other_mat = matrix->create_matrix( clone ) ;
   other_mat->set( matrix ) ;
   other_mat->synchronize() ;
   bool square = matrix->nb_rows()==matrix->nb_cols() ;

   LA_SeqMatrix* ref = LA_DenseMatrix::create( clone, 0, 0 ) ;
   if( dynamic_cast<LA_DenseMatrix const*>( matrix ) != 0 )
   {
      ref = LA_PelMatrix::create( clone, 0, 0 ) ;
   }

   // IO tests
   {
      std::string tmpmat = "__TemporaryMatrix__" ;
      if( matrix->is_desynchronizable() || 
          PEL_Exec::communicator()->rank() == 0 )
      {
         matrix->writeMM( tmpmat ) ;
      }

      // Barrier
      PEL_Exec::communicator()->boolean_and( true ) ;

      ref->readMM( test_matrix ) ;
      bool ok1 = matrix_equality( matrix, ref ) ;
      notify_one_test_result(  label+" readMM", ok1 ) ;

      ref->readMM( tmpmat ) ;
      bool ok2 = matrix_equality( matrix, ref ) ;
      notify_one_test_result(  label+" writeMM", ok2 ) ;

      bool ok3 = matrix_equality( clone, ref ) ;
      notify_one_test_result(  label+" clone", ok3 ) ;

      if( ok1 && ok2 && ok3 ) PEL_System::erase( tmpmat ) ;
   }

   LA_SeqMatrix* other_ref = ref->create_matrix( this ) ;
   other_ref->set( ref ) ;
   other_ref->synchronize() ;

   // Vector and scatter tests
   LA_Vector* vec = matrix->create_vector( clone ) ;
   {
      notify_one_test_result(  label+" create_vector",
                               vec->nb_rows()==matrix->nb_rows() &&
                               vec->implementation()==matrix->implementation()  &&
                               vec->is_synchronized() ) ;
      test_vector( vec ) ;
      LA_Vector* vec2 = vec->create_vector( clone ) ;
      vec2->re_initialize( N ) ;
      for( size_t i=0 ; i<N ; i++ ) vec2->set_item( i, 1.0*i ) ;
      TEST_FAILED( vec2->scale(1.0), "vec unsynchonized after set" ) ;
      vec2->synchronize() ;
      test_vector( vec2 ) ;
   }

   LA_Scatter* S = create_identity_scatter( clone, vec ) ;
   size_t n = vec->nb_rows() ;
   size_t m = matrix->nb_cols() ;
   LA_SeqVector * BV = LA_SeqVector::create( clone, n ) ;
   LA_SeqVector * BV2 = LA_SeqVector::create( clone, n ) ;
   LA_SeqVector * BV3 = LA_SeqVector::create( clone, m ) ;

   // Matrix basic tests
   {
      LA_Vector* y=vec->create_vector(clone) ;
      LA_SeqVector* by=LA_SeqVector::create( clone, n ) ;
      if( !square )
      {
         y->re_initialize( n ) ;
         by->re_initialize( n ) ;
         vec->re_initialize( m ) ;
         BV2->re_initialize( m ) ;
      }

      if( square )
      {
         clone->extract_diag(vec) ;
         S->get(vec,BV) ;
         ref->extract_diag(BV2) ;
         notify_one_test_result(  label+" extract_diag",
                                  vector_equality( BV, BV2) ) ;
      }
      else
      {
         vec->set(1.0) ;
         BV2->set(1.0) ;
      }

      double alpha = 0.1 ;
      double beta = 1.0/3.0 ;

      matrix->multiply_vec_then_add( vec, y, alpha, beta ) ;
      S->get(y,BV) ;
      ref->multiply_vec_then_add( BV2, by, alpha, beta ) ;
      notify_one_test_result(  label+" multiply_vec_then_add",
                               vector_equality( BV, by) ) ;

      LA_Scatter* S2 = create_identity_scatter( clone, vec ) ;
      matrix->tr_multiply_vec_then_add( y, vec, alpha, beta ) ;
      S2->get(vec,BV3) ;
      BV2->synchronize() ;
      ref->tr_multiply_vec_then_add( by, BV2, alpha, beta ) ;
      notify_one_test_result(  label+" tr_multiply_vec_then_add",
                               vector_equality( BV3, BV2 ) ) ;

      double zeta = 1.0/5.0 ;
      clone->multiply_vec_then_add( vec, y ) ;
      ref->multiply_vec_then_add( BV2, by ) ;
      clone->scale( zeta ) ;
      ref->scale( zeta ) ;
      notify_one_test_result(  label+" scale",
                               matrix_equality( clone, ref ) ) ;

      if( !clone->is_symmetric() )
      {
         clone->scale_as_diag_mat_mat(y) ;
         ref->scale_as_diag_mat_mat( by ) ;
         notify_one_test_result(  label+" scale_as_diag_mat_mat",
                                  matrix_equality( clone, ref ) ) ;

         clone->scale_as_mat_diag_mat( vec ) ;
         ref->scale_as_mat_diag_mat( BV2 ) ;
         notify_one_test_result(  label+" scale_as_mat_diag_mat",
                                  matrix_equality( clone, ref ) ) ;
      }
      if( square )
      {
         clone->add_to_diag( vec ) ;
         //clone->synchronize() ;
         ref->add_to_diag( BV2 ) ;
         //ref->synchronize() ;
         notify_one_test_result(  label+" add_to_diag",
                                  matrix_equality( clone, ref ) ) ;
      }

      LA_Matrix*  tclone = create_matrix( clone, m, n ) ;
      LA_SeqMatrix* tref = ref->create_matrix( clone ) ;
      tref->re_initialize( m, n ) ;
      tclone->add_tMat( clone, alpha ) ;
      tref->add_tMat( ref, alpha ) ;
      tclone->synchronize() ;
      tref->synchronize() ;
      notify_one_test_result(  label+" add_tMat",
                               matrix_equality( tclone, tref ) ) ;

      clone->add_Mat( other_mat, alpha ) ;
      //clone->synchronize() ;
      ref->add_Mat( other_ref, alpha ) ;
      //ref->synchronize() ;
      notify_one_test_result(  label+" add_Mat",
                               matrix_equality( clone, ref ) ) ;

      if( !clone->is_symmetric() &&
          !( exp->has_entry( "add_Mat_Mat" ) && !exp->bool_data( "add_Mat_Mat" ) ) )
      {
         LA_Matrix* sclone = create_matrix( clone, n, n ) ;
         LA_SeqMatrix* sref = ref->create_matrix( clone ) ;
         sref->re_initialize( n, n ) ;
         sclone->add_Mat_Mat( other_mat, tclone, alpha ) ;
         sclone->synchronize() ;
         sref->add_Mat_Mat( other_ref, tref, alpha ) ;
         sref->synchronize() ;
         //ref->synchronize() ;
         notify_one_test_result(  label+" add_Mat_Mat",
                                  matrix_equality( sclone, sref ) ) ;
      }

      clone->nullify() ;
      ref->nullify() ;
      clone->synchronize() ;
      ref->synchronize() ;
      notify_one_test_result(  label+" nullify",
                               matrix_equality( clone, ref ) ) ;

      for( size_t i=0 ; i<clone->nb_rows()*clone->nb_cols()/2 ; i++ )
      {
         size_t r = PEL::rand()%clone->nb_rows() ;
         size_t c = PEL::rand()%clone->nb_cols() ;

         if( !clone->is_symmetric() || r<=c )
         {
            clone->set_item(r,c,1.0 ) ;
            ref->set_item(r,c,1.0 ) ;
            if( clone->is_symmetric() && r!=c ) ref->set_item(c,r,1.0 ) ;
         }

      }
      BV2->set(1.0) ;
      vec->set(1.0) ;
      TEST_FAILED( clone->scale(1.0), "mat unsynchonized after set" ) ;
      clone->synchronize() ;
      ref->synchronize() ;
      notify_one_test_result(  label+" set_item",
                               matrix_equality( clone, ref ) ) ;

      for( size_t i=0 ; i<clone->nb_rows()*clone->nb_cols()/2 ; i++ )
      {
         size_t r = PEL::rand()%clone->nb_rows() ;
         size_t c = PEL::rand()%clone->nb_cols() ;
         //out() << "["<<PEL_Exec::communicator()->rank()<<"] r "<<r<<" c "<<c<<std::endl;

         if( !clone->is_symmetric() || r<=c )
         {
            if( !clone->is_desynchronizable() || PEL_Exec::communicator()->rank()==0 )
            {
               clone->add_to_item(r,c,1.0 ) ;
            }

            ref->add_to_item(r,c,1.0 ) ;
            if( clone->is_symmetric() && r!=c ) ref->add_to_item(c,r,1.0 ) ;
         }

      }
      TEST_FAILED( clone->scale(1.0), "mat unsynchonized after add" ) ;
      clone->synchronize() ;
      notify_one_test_result(  label+" add_item",
                               matrix_equality( clone, ref ) ) ;

      // other_mat->re_initialize( clone->nb_rows(), clone->nb_cols() ) ;
      other_mat->set( clone ) ;
      other_mat->nullify() ;
      other_mat->set( clone, true ) ;
      other_mat->synchronize() ;
      notify_one_test_result(  label+" set same_pattern",
                               matrix_equality( other_mat, ref ) ) ;

      other_mat->nullify() ;
      other_mat->add_Mat( clone, 1.0, true ) ;
      other_mat->synchronize() ;
      notify_one_test_result(  label+" add_Mat same_pattern",
                               matrix_equality( other_mat, ref ) ) ;
   }

   if( clone->is_resizable() )
   {
      size_t const nn = 30 ;
      size_t const pp = ( clone->is_symmetric() ? nn : 20 ) ;
      clone->re_initialize( nn, pp ) ;
      notify_one_test_result( label+ " re_initialize",
                              clone->nb_rows()==nn && clone->nb_cols()==pp ) ;
   }

   clone->destroy() ;
}

//-------------------------------------------------------------------------
void
LA_Matrix_TEST:: test_vector( LA_Vector const* vec )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix_TEST:: test_vector" ) ;
   PEL_CHECK( vec != 0 && vec->is_synchronized() ) ;

   LA_Vector* V = vec->create_vector( 0 ) ;
   LA_Vector* V2 = vec->create_vector( V ) ;
   size_t n = V->nb_rows() ;

   LA_SeqVector * BV = LA_SeqVector::create( V, n ) ;
   LA_SeqVector * BV2 = LA_SeqVector::create( V, n ) ;

   for( size_t i=0 ; i<n ; i++)
   {
      BV->set_item(i,0.1*i) ;
   }
   BV->synchronize() ;

   LA_Scatter* S = create_identity_scatter( V, vec ) ;

   std::string label = S->type_name() ;
   notify_one_test_result(  label+" create_scatter",
                            S->size()==n ) ;
   S->set( BV, V ) ;

   notify_one_test_result(  label+" scatter set ",
                            double_equality( V->two_norm(), BV->two_norm() ) ) ;

   S->get( V, BV2 ) ;
   notify_one_test_result(  label+" scatter get",
                            vector_equality( BV, BV2) ) ;

   for( size_t i=0 ; i<n ; i++)
   {
      double d =  1.0/( 1.0+i ) ;
      V->set_item( i, d ) ;
      BV->set_item( i, d ) ;
   }

   label = V->type_name() ;
   if( V->is_desynchronizable() )
   {
      notify_one_test_result(  label+" locally ! synchronized after set",
                               V->state() != LA::Sync ) ;
      notify_one_test_result(  label+" globally ! synchronized after set",
                               V->state() != LA::Sync ) ;
      V->synchronize() ;
      BV->synchronize() ;
      notify_one_test_result(  label+" locally synchronized after synchronize",
                               V->state() == LA::Sync ) ;
      notify_one_test_result(  label+" globally synchronized after synchronize",
                               V->state() == LA::Sync ) ;
   }
   for( size_t i=0 ; i<n ; i++)
   {
      double d =  3.0/( 1.0+i ) ;
      if( !V->is_desynchronizable() || PEL_Exec::communicator()->rank()==0 )
         V->add_to_item( i, d ) ;
      BV->add_to_item( i, d ) ;
   }
   if( V->is_desynchronizable() )
   {
      notify_one_test_result(  label+" globally ! synchronized after set",
                               !V->is_synchronized() ) ;
      V->synchronize() ;
      BV->synchronize() ;
      notify_one_test_result(  label+" locally synchronized after synchronize",
                               V->state() == LA::Sync ) ;
      notify_one_test_result(  label+" globally synchronized after synchronize",
                               V->state() == LA::Sync ) ;
   }
   S->get( V, BV2 ) ;
   notify_one_test_result(  label+" set_item add_to_item",
                            vector_equality( BV, BV2) ) ;

   notify_one_test_result(  label+" dot",
                            double_equality( V->dot(V), BV->dot(BV) ) ) ;

   notify_one_test_result(  label+" two_norm",
                            double_equality( V->two_norm(), BV->two_norm() ) ) ;

   notify_one_test_result(  label+" max_norm",
                            double_equality( V->max_norm(), BV->max_norm() ) ) ;

   V->scale( 1.1 ) ;
   BV->scale( 1.1 ) ;

   notify_one_test_result(  label+" scale",
                            double_equality( V->two_norm(), BV->two_norm() ) ) ;
   V2->nullify() ;
   notify_one_test_result(  label+" nullify",
                            V2->two_norm() == 0.0 ) ;

   V2->re_initialize( V->nb_rows() ) ;
   V2->set( V ) ;
   notify_one_test_result(  label+" set vec",
                            double_equality( V2->two_norm(), V->two_norm() ) ) ;
   V2->set( 1.0 ) ;
   notify_one_test_result(  label+" set double",
                            double_equality( V2->two_norm(), PEL::sqrt(V2->nb_rows()) ) ) ;

   LA_Vector* V3 = V2->create_vector(V) ;
   V3->set_as_v_product( V, V2 ) ;
   notify_one_test_result(  label+" set_as_v_product",
                            double_equality( V3->two_norm(), V->two_norm() ) ) ;

   V3->set_as_reciprocal( V ) ;
   V2->set_as_v_product( V, V3 ) ;
   notify_one_test_result(  label+" set_as_reciprocal",
                            double_equality( V2->two_norm(), PEL::sqrt(V2->nb_rows()) ) ) ;
   V->sum( V, -1.0 ) ;
   notify_one_test_result(  label+" sum",
                            double_equality( V->two_norm(), 0.0 ) ) ;

   if( PEL_Exec::communicator()->nb_ranks()==0 || vec->is_desynchronizable() )
   {
      std::string tmp( "vector.txt" ) ;

      vec->write( tmp ) ;
      LA_Vector* V4 = V->create_vector(V) ;
      V4->read( tmp ) ;
      bool ok = double_equality( V4->two_norm(), vec->two_norm() ) ;
      notify_one_test_result(  label+" read/write", ok ) ;

      LA_Vector* V5 = LA_SeqVector::create(V,0) ;
      V5->read( tmp ) ;
      ok = double_equality( V5->two_norm(), vec->two_norm() ) ;
      notify_one_test_result(  label+" write same format", ok ) ;
      if( ok )
      {
         if( PEL_Exec::communicator()->rank()==0 )
         {
            PEL_System::erase( tmp.c_str() ) ;
         }
      }
      PEL_Exec::communicator()->barrier() ;
   }
   bool res ;
   if( V->is_desynchronizable() )
   {
      if( PEL_Exec::communicator()->rank()==0 )
      {
         V->add_to_item(0,0.0) ;
         res = ( V->state() != LA::Sync ) ;
      }
      else
      {
         res = ( V->state() == LA::Sync ) ;
      }
      notify_one_test_result(  label+" is locally synchronized after set on 0", res ) ;
      res = !V->is_synchronized() ;
      notify_one_test_result(  label+" globally ! synchronized after set on 0",
                               res ) ;
   }

   V->destroy() ;
}

//-------------------------------------------------------------------------
LA_Scatter*
LA_Matrix_TEST:: create_identity_scatter(
                                PEL_Object* a_owner, LA_Vector const* vec )
//-------------------------------------------------------------------------
{
   size_t n = vec->nb_rows() ;

   size_t_vector from(n) ;
   size_t_vector to(n) ;

   for( size_t i=0 ; i<n ; i++)
   {
      // Identity
      from(i) = i ;
      to(i) = i ;
   }
   LA_Scatter* result = vec->create_scatter( a_owner, from, to ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
LA_Matrix*
LA_Matrix_TEST:: create_matrix(
                            PEL_Object* a_owner, size_t n, size_t m ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix_TEST:: create_matrix" ) ;

   if( !MAT_PROTO->is_resizable() )
   {
      PEL_Error::object()->raise_internal( "not implemented" ) ;
   }

   LA_Matrix* result = MAT_PROTO->create_matrix( a_owner ) ;
   result->re_initialize( n, m ) ;

   return( result ) ;
}
