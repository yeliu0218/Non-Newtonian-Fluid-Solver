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

#include <LA_Matrix.hh>

#include <fstream>
#include <iomanip>

#include <LA_SeqMatrix.hh>
#include <LA_MatrixIterator.hh>
#include <LA_Vector.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Communicator.hh>
#include <PEL_DistributedPartition.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>

#ifdef OUTLINE
#define inline
   #include <LA_Matrix.icc>
#undef inline
#endif

//----------------------------------------------------------------------
LA_Matrix*
LA_Matrix:: make( PEL_Object* a_owner,
                  PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: make" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   std::string name = exp->string_data( "concrete_name" ) ;
   LA_Matrix const* proto =
      static_cast<LA_Matrix const*>( plugins_map()->item( name ) ) ;

   PEL_ASSERT( proto!=0 ) ;
   PEL_ASSERT( proto->is_a_prototype() ) ;

   LA_Matrix* result = proto->create_replica( a_owner, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == exp->string_data( "concrete_name" ) ) ;
   PEL_CHECK_POST( IMPLIES( result->is_desynchronizable(),
                            result->state() == LA::NotSync_undef ) ) ;
   PEL_CHECK_POST( result->distribution_strategy() !=
                   LA::InvalidDistribution ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Matrix:: LA_Matrix( PEL_Object* a_owner, 
                       std::string const& a_name,
                       LA::DistributionStrategy a_dist_strat )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , IS_PROTO( false )
   , NAME( a_name )
   , SYNC_STATE( LA::Sync )
   , RESIZABLE( true )
   , ONLY_LOCAL_MODIFS( false )
   , DIST_STRAT( a_dist_strat )
{
   PEL_LABEL( "LA_Matrix:: LA_Matrix" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;

   PEL_CHECK_POST( !is_a_prototype() ) ;
   PEL_CHECK_POST( owner() == a_owner ) ;
   PEL_CHECK_POST( name() == a_name ) ;
   PEL_CHECK_POST( is_resizable() ) ;
   PEL_CHECK_POST( distribution_strategy() == a_dist_strat ) ;
}

//----------------------------------------------------------------------
LA_Matrix:: LA_Matrix( std::string const& a_name )
//----------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , IS_PROTO( true )
   , NAME( a_name )
   , SYNC_STATE( LA::Sync )
   , RESIZABLE( false )
   , ONLY_LOCAL_MODIFS( false )
   , DIST_STRAT( LA::InvalidDistribution )
{
   PEL_LABEL( "LA_Matrix:: LA_Matrix( prototype )" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;

   plugins_map()->register_item( a_name, this ) ;

   PEL_CHECK_POST( is_a_prototype() ) ;
   PEL_CHECK_POST( name() == a_name ) ;
   PEL_CHECK_POST( ! is_resizable() ) ;
   PEL_CHECK_POST( distribution_strategy() == LA::InvalidDistribution ) ;
}

//----------------------------------------------------------------------
LA_Matrix:: ~LA_Matrix( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
std::string const&
LA_Matrix:: name( void ) const
//----------------------------------------------------------------------
{
   return( NAME ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: is_resizable( void ) const
//----------------------------------------------------------------------
{
   return( RESIZABLE ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: is_symmetric( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: is_symmetric" ) ;

   bool result = false ;

   PEL_CHECK_POST( is_symmetric_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_Matrix:: make_non_resizable( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: make_non_resizable" ) ;

   RESIZABLE = false ;

   PEL_CHECK_POST( !is_resizable() ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: only_local_modifs( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: only_local_modifs" ) ;

   bool result = ONLY_LOCAL_MODIFS ;

   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_Matrix:: set_only_local_modifs_state( bool only_local )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: set_only_local_modifs_state" ) ;

   ONLY_LOCAL_MODIFS = only_local ;

   PEL_CHECK_POST( only_local_modifs() == only_local ) ;
}

//----------------------------------------------------------------------
void
LA_Matrix:: start_local_modifs( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: start_local_modifs" ) ;
   PEL_CHECK_PRE( start_local_modifs_PRE() ) ;

   set_only_local_modifs_state( true ) ;

   PEL_CHECK_POST( start_local_modifs_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_Matrix:: stop_local_modifs( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: stop_local_modifs" ) ;
   PEL_CHECK_PRE( stop_local_modifs_PRE() ) ;

   set_only_local_modifs_state( false ) ;

   PEL_CHECK_POST( stop_local_modifs_POST() ) ;
}

//----------------------------------------------------------------------
LA::DistributionStrategy
LA_Matrix:: distribution_strategy( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: distribution_strategy" ) ;

   return( DIST_STRAT ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: is_synchronized( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: is_synchronized" ) ;
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;

   bool result = ( is_desynchronizable() ?
                      PEL_Exec::communicator()->boolean_and(
                                                SYNC_STATE == LA::Sync ) :
                      true ) ;

   PEL_CHECK_POST(
      IMPLIES( !is_desynchronizable(), result == true ) ) ;
   PEL_CHECK_POST(
      IMPLIES( is_desynchronizable(),
               result == PEL_Exec::communicator()->boolean_and(
                                               state() == LA::Sync ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
LA_Matrix:: nb_local_rows( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: nb_local_rows" ) ;

   size_t result = row_distribution()->local_number() ;

   PEL_CHECK_POST( result == row_distribution()->local_number() ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
size_t
LA_Matrix:: nb_local_cols( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: nb_local_cols" ) ;

   size_t result = col_distribution()->local_number() ;

   PEL_CHECK_POST( result == col_distribution()->local_number() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_Matrix:: synchronize( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: synchronize" ) ;
   PEL_CHECK_PRE( synchronize_PRE() ) ;

   SYNC_STATE = LA::Sync ;

   PEL_CHECK_POST( synchronize_POST() ) ;
}

//----------------------------------------------------------------------
LA_MatrixIterator*
LA_Matrix:: create_stored_item_iterator( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: create_stored_item_iterator" ) ;
   PEL_CHECK_PRE( create_stored_item_iterator_PRE( a_owner ) ) ;

   LA_SeqMatrix* matrix = create_local_matrix( 0 ) ;
   LA_MatrixIterator* result = matrix->create_stored_item_iterator( a_owner ) ;
   matrix->set_owner( result ) ;

   PEL_CHECK_POST( create_stored_item_iterator_POST( a_owner, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
LA_Matrix:: extract_diag( LA_Vector* diag ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: extract_diag" ) ;
   PEL_CHECK_PRE( extract_diag_PRE( diag ) ) ;

   PEL_Error::object()->raise_not_implemented( this, "extract_diag" ) ;

   PEL_CHECK_POST( extract_diag_POST( diag ) ) ;
}

//----------------------------------------------------------------------
void
LA_Matrix:: tr_multiply_vec_then_add(
                        LA_Vector const* x, LA_Vector* y,
                        double alpha, double beta )  const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: tr_multiply_vec_then_add" ) ;
   PEL_CHECK_PRE( tr_multiply_vec_then_add_PRE( x, y, alpha, beta ) ) ;

   PEL_Error::object()->raise_not_implemented( this,
                                               "tr_multiply_vec_then_add" ) ;

   PEL_CHECK_POST( tr_multiply_vec_then_add_POST( y ) ) ;
}


//----------------------------------------------------------------------
void
LA_Matrix:: scale_as_diag_mat_mat( LA_Vector const* lvec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: scale_as_diag_mat_mat" ) ;
   PEL_CHECK_PRE( scale_as_diag_mat_mat_PRE( lvec ) ) ;

   PEL_Error::object()->raise_not_implemented( this,
                                               "scale_as_diag_mat_mat" ) ;

   PEL_CHECK_POST( scale_as_diag_mat_mat_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_Matrix:: scale_as_mat_diag_mat( LA_Vector const* rvec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: scale_as_mat_diag_mat" ) ;
   PEL_CHECK_PRE( scale_as_mat_diag_mat_PRE( rvec ) ) ;

   PEL_Error::object()->raise_not_implemented( this,
                                               "scale_as_mat_diag_mat" ) ;

   PEL_CHECK_POST( scale_as_mat_diag_mat_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_Matrix:: add_to_diag( LA_Vector const* vec )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: add_to_diag" ) ;
   PEL_CHECK_PRE( add_to_diag_PRE( vec ) ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;

   PEL_Error::object()->raise_not_implemented( this, "add_to_diag" ) ;

   PEL_CHECK_POST( add_to_diag_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_Matrix:: set( LA_Matrix const* A, bool same_pattern )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: set" ) ;
   PEL_CHECK_PRE( set_PRE( A, same_pattern ) ) ;

   nullify() ;
   add_Mat( A, 1.0, same_pattern ) ;
   //synchronize() ;

   PEL_CHECK_POST( set_POST( A ) ) ;
}

//----------------------------------------------------------------------
void
LA_Matrix:: add_tMat( LA_Matrix const* A, double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: add_tMat" ) ;
   PEL_CHECK_PRE( add_tMat_PRE( A, alpha ) ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;

   PEL_Error::object()->raise_not_implemented( this, "add_tMat" ) ;

   PEL_CHECK_POST( add_tMat_POST( OLD( state ) ) ) ;
}

//----------------------------------------------------------------------
void
LA_Matrix:: add_Mat_Mat( LA_Matrix const* A, LA_Matrix const* B,
                         double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: add_Mat_Mat" ) ;
   PEL_CHECK_PRE( add_Mat_Mat_PRE( A, B, alpha ) ) ;
   PEL_SAVEOLD( LA::SyncState, state, state() ) ;

   PEL_Error::object()->raise_not_implemented( this, "add_Mat_Mat" ) ;

   PEL_CHECK_POST( add_Mat_Mat_POST( OLD( state ) ) ) ;
}
//----------------------------------------------------------------------
void
LA_Matrix:: readMM( std::string const& file )
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: readMM" ) ;
   PEL_CHECK_PRE( readMM_PRE( file ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( !is_resizable() )
   {
      PEL_Error::object()->raise_internal(
         "Reading of non resizable matrix is not implemented" ) ;
   }

   std::ifstream in( file.c_str(), std::ios_base::binary ) ; //???
   if( !in )
   {
      PEL_Error::object()->raise_file_handling( file, "open" ) ;
   }
   bool array, sym ;
   size_t n, m, p ;
   PEL_ASSERT( readMMHeader( in, array, sym, n, m, p ) ) ;
   if( n*m == 0 && ( n>0 || m>0 ) )
   {
      PEL_Error::object()->raise_plain(
         "*** LA_Matrix:: readMM error:\n"
         "    Bad number of row and column in file "+file ) ;
   }
   if( is_symmetric() && !sym )
   {
      PEL_Error::object()->raise_plain(
         "*** LA_Matrix:: readMM error:\n"
         "    No symmetric matrix in file "+file+"\n"
         "    can't be read by a symmetric one" ) ;
   }
   if( sym && n!=m )
   {
      PEL_Error::object()->raise_plain(
         "*** LA_Matrix:: readMM error:\n"
         "    Symmetric matrix not square in file "+file ) ;
   }

   re_initialize( n, m ) ;

   double v ;
   size_t nb = ( is_desynchronizable() ? 
                    PEL_Exec::communicator()->nb_ranks() : 1 ) ; //???
   size_t r =  ( is_desynchronizable() ? 
                    PEL_Exec::communicator()->rank() : 0 ) ; //???
   size_t s = p / nb ;
   size_t first = r*s ;
   size_t last = (r+1)*s ;
   if( r==nb-1 ) last = p ;

   if( !array )
   {
      size_t i, j ;
      for( size_t k = 0 ; k < last ; k++ )
      {
         in >> i >> j >> v ;
         if( k >= first )
         {
            i-- ;
            j-- ;
            if( sym && i<j )
            {
               PEL_Error::object()->raise_plain(
                  "*** LA_Matrix:: readMM error:\n"
                  "    Symmetric matrix in file "+file+"\n"
                  "    should contain only terms below the main diagonal" ) ;
            }

            if( is_symmetric() )
            {
               set_item( j, i, v ) ;
            }
            else
            {
               set_item( i, j, v ) ;
               if( sym && i!=j )
               {
                  set_item( j, i, v ) ;
               }
            }
         }
      }
   }
   else
   {
      for( size_t j=0 ; j<nb_cols() ; j++ )
      {
         size_t lower = 0 ;
         if( sym )
         {
            lower = j ;
         }
         for( size_t i=lower; i<nb_rows() ; i++ )
         {
            in >> v ;
            if( is_symmetric() )
            {
               set_item( j, i, v ) ;
            }
            else
            {
               set_item( i, j, v ) ;
               if( sym && i!=j )
               {
                  set_item( j, i, v ) ;
               }
            }
         }
      }
   }
   synchronize() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( readMM_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_Matrix:: writeMM( std::string const& file ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: writeMM" ) ;
   PEL_CHECK_PRE( writeMM_PRE( file ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t rank = PEL_Exec::communicator()->rank() ;
   size_t size = PEL_Exec::communicator()->nb_ranks() ;
   int dummy ;
   LA_SeqMatrix* mat = create_local_matrix(0) ;
   LA_MatrixIterator* it = mat->create_stored_item_iterator( mat ) ;
   size_t total_nb = mat->nb_stored_items() ;
   if( is_symmetric() )
   {
      total_nb= 0 ;
      for( it->start_all_items() ; it->is_valid() ; it->go_next() )
      {
         if(it->row()<=it->col()) total_nb++ ;
      }
   }
   else
   {
      total_nb = mat->nb_stored_items() ;
   }

   total_nb = (size_t)( PEL_Exec::communicator()->sum( total_nb ) ) ;

   if( rank>0 ) PEL_Exec::communicator()->receive( rank-1, dummy ) ;

   std::ofstream out(
      file.c_str(),
      ( rank==0 ? std::ios::out : std::ios::out|std::ios::app ) ) ;
   if( !out )
   {
      PEL_Error::object()->raise_file_handling( file, "open" ) ;
   }
   if( rank==0 )
   {
      out << "%%MatrixMarket matrix coordinate real "
          << ( is_symmetric() ? "symmetric" : "general" ) ;
      out << std::endl << nb_rows() << " " << nb_cols()
          << " " << total_nb << std::endl ;
   }

   for( it->start_all_items() ; it->is_valid() ; it->go_next() )
   {
      if( !is_symmetric() )
      {
         out << std::setprecision(16) << it->row()+1
             << " " << it->col()+1 << " " << it->item() << std::endl ;
      }
      else if(it->row()<=it->col())
      {
         out << std::setprecision(16)  << it->col()+1
             << " " << it->row()+1 << " " << it->item() << std::endl ;
      }
   }
   mat->destroy() ;

   out.close() ;

   if( rank!=size-1 )
   {
      PEL_Exec::communicator()->send( rank+1, dummy ) ;
      PEL_Exec::communicator()->receive( size-1, dummy ) ;
   }
   else
   {
      for( size_t i=0 ; i<size-1 ; ++i )
         PEL_Exec::communicator()->send( i, dummy ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( writeMM_POST() ) ;
}

//----------------------------------------------------------------------
void
LA_Matrix:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string s( indent_width, ' ' ) ;
   os << s << "matrix: \"" << name() << "\"" << std::endl ;
}

//----------------------------------------------------------------------
void
LA_Matrix:: print_items( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: print_items" ) ;
   PEL_CHECK_PRE( print_items_PRE( os, indent_width ) ) ;

   std::string space( indent_width, ' ' ) ;
   os << space << "nb_rows:" << nb_rows() << " nb_cols:" << nb_cols()
      << std::endl ;
   if( nb_rows()>0 && nb_cols()>0 )
   {
      std::ios_base::fmtflags original_flags = os.flags() ;
      os.setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
      std::streamsize p = os.precision() ;
      os << std::setprecision( 7 ) ;
      LA_MatrixIterator* it = create_stored_item_iterator( 0 ) ;
      for( size_t i=0 ; i<nb_rows() ; i++ )
      {
         bool first = true ;
         for( it->start_row_items(i) ; it->is_valid() ; it->go_next() )
         {
            if( first )
            {
               os << space << "Rows n°" << i << std::endl ;
               first = false ;
            }
            os << space
               << "   ( i=" << i << "  j=" << it->col()
               << "  val=" << std::setw(15) << it->item() << " )"
               << std::endl ;
         }
      }
      it->destroy() ; it = 0 ;
      os << std::setprecision(p) ;
      os.flags( original_flags ) ;
   }
}

//----------------------------------------------------------------------
bool
LA_Matrix:: is_a_prototype( void ) const
//----------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: readMMHeader( std::ifstream& in,
                          bool& array,
                          bool& sym,
                          size_t& n,
                          size_t& m,
                          size_t& nb ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "LA_Matrix:: readMMHeader" ) ;
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_CHECK_PRE( in ) ;

   std::string error_mess ;

   if( !in )
   {
      error_mess = "Reading matrix header : file not opened " ;
   }
   else
   {
      std::string label, storage, content, type, kind, str ;
      in >> label >> content >> storage >> type >> kind ;
      getline( in, str ) ;

      if( ( label != "%%MatrixMarket" )
          || ( content != "matrix" )
          || ( type != "real" )
          || ( !( storage == "coordinate" || storage == "array" ) )
          || ( !( kind == "general" || kind == "symmetric" ) ) )
      {
         error_mess = "Reading matrix header : " ;

         if( label != "%%MatrixMarket" )
         {
            error_mess += " Bad label : " + label ;
         }
         if( content != "matrix" )
         {
            error_mess += " Bad content : " + content ;
         }
         if( type != "real" )
         {
            error_mess += " Bad type : " + type ;
         }
         if( !( storage == "coordinate" || storage == "array" ) )
         {
            error_mess += " Bad storage : " + storage ;
         }
         if( !( kind == "general" || kind == "symmetric" ) )
         {
            error_mess += " Bad matrix kind : " + kind ;
         }
      }
      else
      {
         sym = ( kind == "symmetric" ) ;
         array = ( storage == "array" ) ;
         std::string ss ;
         std::streampos pos = in.tellg() ;

         while( !in.eof() )
         {
            pos = in.tellg() ;
            getline( in, ss ) ;

            if( ss[0] != '%' )
            {
               break ;
            }
         }
         in.seekg( pos ) ;
         if( storage == "coordinate" )
         {
            in >> n >> m >> nb ;
         }
         else
         {
            in >> n >> m ;
            nb = n*m ;
         }
      }
   }
   if( !error_mess.empty() )
   {
      PEL_Error::object()->raise_plain( error_mess ) ;
   }

   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: create_matrix_PRE( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: create_matrix_POST(
                          LA_Matrix* result, PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( result->name()==name() ) ;
   PEL_ASSERT( result->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( result->nb_cols() == nb_cols() ) ;
   PEL_ASSERT( result->implementation() == implementation() ) ;
   PEL_ASSERT( result->name() == name() ) ;
   PEL_ASSERT( result->is_resizable() == is_resizable() ) ;
   PEL_ASSERT( IMPLIES( result->is_desynchronizable(),
                        result->state() == LA::NotSync_undef ) ) ;
   PEL_ASSERT( result->row_distribution()->is_compatible(
                                             row_distribution() ) ) ;
   PEL_ASSERT( result->col_distribution()->is_compatible(
                                             col_distribution() ) ) ;
   PEL_ASSERT( IMPLIES( result->is_desynchronizable(),
                        !result->is_synchronized() ) ) ;
   PEL_ASSERT( result->distribution_strategy() == distribution_strategy() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: create_vector_PRE( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: create_vector_POST( LA_Vector const* result, PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( result->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( result->implementation() == implementation() ) ;
   PEL_ASSERT( result->is_resizable() == is_resizable() ) ;
   PEL_ASSERT( result->is_desynchronizable() == is_desynchronizable() ) ;
   PEL_ASSERT( result->state() == LA::Sync ) ;
   PEL_ASSERT(
      result->row_distribution()->is_compatible( row_distribution() ) ) ;
   PEL_ASSERT( result->is_synchronized() ) ;
   PEL_ASSERT( result->distribution_strategy() == distribution_strategy() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: re_initialize_PRE(
                  size_t a_nb_rows, size_t a_nb_cols,
                  size_t a_nb_local_rows, size_t a_nb_local_cols ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( is_resizable() ) ;

   PEL_ASSERT( IMPLIES( distribution_strategy() != LA::FromLocalSize,
                        a_nb_rows != PEL::bad_index() ) ) ;
   PEL_ASSERT( IMPLIES( distribution_strategy() != LA::FromLocalSize,
                        a_nb_cols != PEL::bad_index() ) ) ;
   
   PEL_ASSERT( IMPLIES( distribution_strategy() == LA::FromLocalSize,
                        a_nb_local_rows != PEL::bad_index() ) ) ;
   PEL_ASSERT( IMPLIES( distribution_strategy() == LA::FromLocalSize,
                        a_nb_local_cols != PEL::bad_index() ) ) ;
   
   PEL_ASSERT(
     IMPLIES( distribution_strategy() == LA::FromGlobalSize,
        PEL_Exec::communicator()->same_value_everywhere( (int) a_nb_rows ) ) ) ;
   PEL_ASSERT(
     IMPLIES( distribution_strategy() == LA::FromGlobalSize,
        PEL_Exec::communicator()->same_value_everywhere( (int) a_nb_cols ) ) ) ;

   PEL_ASSERT(
     IMPLIES( is_symmetric() && distribution_strategy() != LA::FromLocalSize ,
              a_nb_rows == a_nb_cols ) ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: re_initialize_POST( size_t a_nb_rows, size_t a_nb_cols,
                                size_t a_nb_local_rows, size_t a_nb_local_cols ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( distribution_strategy() != LA::FromLocalSize,
                        nb_rows() == a_nb_rows ) ) ;
   PEL_ASSERT( IMPLIES( distribution_strategy() != LA::FromLocalSize,
                        nb_cols() == a_nb_cols ) ) ;
   PEL_ASSERT( IMPLIES( distribution_strategy() == LA::FromLocalSize,
                        nb_local_rows() ==  a_nb_local_rows ) ) ;
   PEL_ASSERT( IMPLIES( distribution_strategy() == LA::FromLocalSize,
                        nb_local_cols() ==  a_nb_local_cols ) ) ;
   PEL_ASSERT( IMPLIES( is_desynchronizable(),
                        state() == LA::NotSync_undef ) ) ;
   PEL_ASSERT( IMPLIES( is_desynchronizable(), !is_synchronized() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: implementation_POST( LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: is_symmetric_POST( bool result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( result, nb_rows() == nb_cols() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: start_local_modifs_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( !only_local_modifs() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: start_local_modifs_POST( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( only_local_modifs() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: stop_local_modifs_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( only_local_modifs() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: stop_local_modifs_POST( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( !only_local_modifs() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: synchronize_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( IMPLIES( is_desynchronizable(), !is_synchronized() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: synchronize_POST( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT( state() == LA::Sync ) ;
   return( true ) ;
}


//----------------------------------------------------------------------
bool
LA_Matrix:: row_distribution_PRE( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: row_distribution_POST(
                          PEL_DistributedPartition const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == this ) ;
   PEL_ASSERT( result->partitioning().size()==
                               PEL_Exec::communicator()->nb_ranks() ) ;
   PEL_ASSERT( result->global_number()==nb_rows() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: col_distribution_PRE( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: col_distribution_POST(
                          PEL_DistributedPartition const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == this ) ;
   PEL_ASSERT( result->partitioning().size()==
                               PEL_Exec::communicator()->nb_ranks() ) ;
   PEL_ASSERT( result->global_number()==nb_cols() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: create_stored_item_iterator_PRE( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: create_stored_item_iterator_POST(
                  PEL_Object* a_owner, LA_MatrixIterator* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( result->is_valid() == false ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: create_local_matrix_PRE( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: create_local_matrix_POST( LA_SeqMatrix const* result,
                                      PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( result->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( result->nb_cols() == nb_cols() ) ;
   PEL_ASSERT( result->is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: nb_stored_items_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: extract_diag_PRE( LA_Vector const* diag ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT( nb_rows() == nb_cols() ) ;
   PEL_ASSERT( diag != 0 ) ;
   PEL_ASSERT( diag->nb_rows() == nb_cols() ) ;
   PEL_ASSERT( diag->implementation() == implementation() ) ;
   PEL_ASSERT(
      diag->row_distribution()->is_compatible( row_distribution() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: extract_diag_POST( LA_Vector const* diag ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( diag->is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: set_item_PRE( size_t i, size_t j ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( i < nb_rows() ) ;
   PEL_ASSERT( j < nb_cols() ) ;
   PEL_ASSERT( IMPLIES( only_local_modifs(),
                        i >= row_distribution()->first_local_index() &&
                        i < row_distribution()->local_index_limit() ) ) ;
   PEL_ASSERT( state() != LA::NotSync_add ) ;
   PEL_ASSERT( IMPLIES( is_symmetric(), i<=j ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: set_item_POST( size_t i, size_t j,
                           LA::SyncState old_state ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT(
   IMPLIES( is_desynchronizable() && !only_local_modifs(),
            state() == LA::NotSync_set ) ) ;
   PEL_ASSERT(
   IMPLIES( !is_desynchronizable() || only_local_modifs(),
            state() == old_state ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: add_to_item_PRE( size_t i, size_t j ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( i < nb_rows() ) ;
   PEL_ASSERT( j < nb_cols() ) ;
   PEL_ASSERT( IMPLIES( only_local_modifs(),
                        i >= row_distribution()->first_local_index() &&
                        i < row_distribution()->local_index_limit() ) ) ;
   PEL_ASSERT( state() != LA::NotSync_set ) ;
   PEL_ASSERT( IMPLIES( is_symmetric(), i<=j ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: add_to_item_POST( size_t i, size_t j,
                              LA::SyncState old_state ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT(
      IMPLIES( is_desynchronizable() && !only_local_modifs(),
               state() == LA::NotSync_add ) ) ;
   PEL_ASSERT(
   IMPLIES( !is_desynchronizable() || only_local_modifs(),
            state() == old_state ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: nullify_PRE( void ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: nullify_POST( void ) const
//----------------------------------------------------------------------
{
   //??? Matrix is not synchronized because of PETSc:
   //???   synchronization freezes matrix pattern and avoids assembling
   PEL_ASSERT( IMPLIES( is_desynchronizable(),
                        state() == LA::NotSync_undef ) ) ;
   PEL_ASSERT( IMPLIES( is_desynchronizable(),
                        !is_synchronized() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: scale_PRE( double alpha ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT(
      IMPLIES( is_desynchronizable(),
               PEL_Exec::communicator()->same_value_everywhere( alpha ) ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: scale_POST( double alpha ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: multiply_vec_then_add_PRE( LA_Vector const* x,
                                       LA_Vector* y,
                                       double alpha,
                                       double beta ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT( x != 0 ) ;
   PEL_ASSERT( x->nb_rows() == nb_cols() ) ;
   PEL_ASSERT( x->implementation() == implementation() ) ;
   PEL_ASSERT( x->is_synchronized() ) ;
   PEL_ASSERT( x->row_distribution()->is_compatible( col_distribution() ) ) ;
   PEL_ASSERT( y != 0 ) ;
   PEL_ASSERT( y->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( y->implementation() == implementation() ) ;
   PEL_ASSERT( y->is_synchronized() ) ;
   PEL_ASSERT( y->row_distribution()->is_compatible( row_distribution() ) ) ;
   PEL_ASSERT( x != y ) ;
   PEL_ASSERT(
      IMPLIES( is_desynchronizable(),
               PEL_Exec::communicator()->same_value_everywhere( alpha ) ) ) ;
   PEL_ASSERT( alpha != 0. ) ;
   PEL_ASSERT(
      IMPLIES( is_desynchronizable(),
               PEL_Exec::communicator()->same_value_everywhere( beta ) ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: multiply_vec_then_add_POST( LA_Vector* y ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( y->is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: tr_multiply_vec_then_add_PRE( LA_Vector const* x,
                                          LA_Vector* y,
                                          double alpha,
                                          double beta ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT( x != 0 ) ;
   PEL_ASSERT( x->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( x->implementation() == implementation() ) ;
   PEL_ASSERT( x->is_synchronized() ) ;
   PEL_ASSERT( x->row_distribution()->is_compatible( row_distribution() ) ) ;
   PEL_ASSERT( y != 0 ) ;
   PEL_ASSERT( y->nb_rows() == nb_cols() ) ;;
   PEL_ASSERT( y->implementation() == implementation() ) ;
   PEL_ASSERT( y->is_synchronized() ) ;
   PEL_ASSERT( y->row_distribution()->is_compatible( col_distribution() ) ) ;
   PEL_ASSERT( x != y ) ;
   PEL_ASSERT(
      IMPLIES( is_desynchronizable(),
               PEL_Exec::communicator()->same_value_everywhere( alpha ) ) ) ;
   PEL_ASSERT( alpha != 0. ) ;
   PEL_ASSERT(
      IMPLIES( is_desynchronizable(),
               PEL_Exec::communicator()->same_value_everywhere( beta ) ) ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: tr_multiply_vec_then_add_POST( LA_Vector* y ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( y->is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: scale_as_diag_mat_mat_PRE( LA_Vector const* lvec ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT( !is_symmetric() ) ;
   PEL_ASSERT( lvec != 0 ) ;
   PEL_ASSERT( lvec->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( lvec->implementation() == implementation() ) ;
   PEL_ASSERT( lvec->is_synchronized() ) ;
   PEL_ASSERT( lvec->row_distribution()->is_compatible( row_distribution() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: scale_as_diag_mat_mat_POST( void) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: scale_as_mat_diag_mat_PRE( LA_Vector const* rvec ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( is_synchronized() ) ;
   PEL_ASSERT( !is_symmetric() ) ;
   PEL_ASSERT( rvec != 0 ) ;
   PEL_ASSERT( rvec->nb_rows() == nb_cols() ) ;
   PEL_ASSERT( rvec->implementation() == implementation() ) ;
   PEL_ASSERT( rvec->is_synchronized() ) ;
   PEL_ASSERT( rvec->row_distribution()->is_compatible( col_distribution() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: scale_as_mat_diag_mat_POST( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: add_to_diag_PRE( LA_Vector const* vec ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( state() != LA::NotSync_set ) ;
   PEL_ASSERT( nb_cols() == nb_rows() ) ;
   PEL_ASSERT( vec != 0 ) ;
   PEL_ASSERT( vec->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( vec->implementation() == implementation() ) ;
   PEL_ASSERT( vec->is_synchronized() ) ;
   PEL_ASSERT( vec->row_distribution()->is_compatible( row_distribution() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: add_to_diag_POST( LA::SyncState old_state ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( state() == old_state ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: set_PRE( LA_Matrix const* A, bool same_pattern ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( A != 0 ) ;
   PEL_ASSERT( A->nb_rows()==nb_rows() ) ;
   PEL_ASSERT( A->nb_cols()==nb_cols() ) ;
   PEL_ASSERT( A->row_distribution()->is_compatible( row_distribution() ) ) ;
   PEL_ASSERT( A->col_distribution()->is_compatible( col_distribution() ) ) ;
   PEL_ASSERT( A->is_synchronized() ) ;
   PEL_ASSERT( A->implementation()==implementation() ) ;
   PEL_ASSERT( IMPLIES( is_symmetric() , A->is_symmetric() ) ) ;
   PEL_ASSERT( !is_a_prototype() ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: set_POST( LA_Matrix const* A ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( is_desynchronizable(), state() == LA::NotSync_undef ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: add_Mat_PRE( LA_Matrix const* A, double alpha,
                         bool same_pattern ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( state() != LA::NotSync_set ) ;
   PEL_ASSERT( A != 0 ) ;
   PEL_ASSERT( A != this ) ;
   PEL_ASSERT( A->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( A->nb_cols() == nb_cols() ) ;
   PEL_ASSERT( A->implementation() == implementation() ) ;
   PEL_ASSERT( A->is_synchronized() ) ;
   PEL_ASSERT( A->row_distribution()->is_compatible( row_distribution() ) ) ;
   PEL_ASSERT( A->col_distribution()->is_compatible( col_distribution() ) ) ;
   PEL_ASSERT(
       IMPLIES( is_desynchronizable(),
                PEL_Exec::communicator()->same_value_everywhere( alpha ) ) ) ;
   PEL_ASSERT( alpha != 0. ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: add_Mat_POST( LA::SyncState old_state ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( state() == old_state ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: add_tMat_PRE( LA_Matrix const* A,
                          double alpha ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( state() != LA::NotSync_set ) ;
   PEL_ASSERT( A != 0 ) ;
   PEL_ASSERT( A != this ) ;
   PEL_ASSERT( A->nb_rows() == nb_cols() ) ;
   PEL_ASSERT( A->nb_cols() == nb_rows() ) ;
   PEL_ASSERT( A->implementation() == implementation() ) ;
   PEL_ASSERT( A->is_synchronized() ) ;
   PEL_ASSERT( A->row_distribution()->is_compatible( col_distribution() ) ) ;
   PEL_ASSERT( A->col_distribution()->is_compatible( row_distribution() ) ) ;
   PEL_ASSERT(
       IMPLIES( is_desynchronizable(),
                PEL_Exec::communicator()->same_value_everywhere( alpha ) ) ) ;
   PEL_ASSERT( alpha != 0. ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: add_tMat_POST( LA::SyncState old_state ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( state() == old_state ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: add_Mat_Mat_PRE( LA_Matrix const* A,
                             LA_Matrix const* B,
                             double alpha ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( state() != LA::NotSync_set ) ;
   PEL_ASSERT( A != 0 ) ;
   PEL_ASSERT( A != this ) ;
   PEL_ASSERT( A->nb_rows() == nb_rows() ) ;
   PEL_ASSERT( A->implementation() == implementation() ) ;
   PEL_ASSERT( A->is_synchronized() ) ;
   PEL_ASSERT( A->row_distribution()->is_compatible( row_distribution() ) ) ;
   PEL_ASSERT( B != 0 ) ;
   PEL_ASSERT( B != this ) ;
   PEL_ASSERT( B->nb_rows() == A->nb_cols() ) ;
   PEL_ASSERT( B->nb_cols() == nb_cols() ) ;
   PEL_ASSERT( B->implementation() == implementation() ) ;
   PEL_ASSERT( B->is_synchronized() ) ;
   PEL_ASSERT( B->row_distribution()->is_compatible( A->col_distribution() ) ) ;
   PEL_ASSERT( B->col_distribution()->is_compatible( col_distribution() ) ) ;
   PEL_ASSERT(
      IMPLIES( is_desynchronizable(),
               PEL_Exec::communicator()->same_value_everywhere( alpha ) ) ) ;
   PEL_ASSERT( alpha != 0. ) ;
   PEL_ASSERT( !is_symmetric() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: add_Mat_Mat_POST( LA::SyncState old_state ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( state() == old_state ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: readMM_PRE( std::string const& file ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( !file.empty() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: readMM_POST( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: writeMM_PRE( std::string const& file ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_COLLECTIVE( is_desynchronizable() ) ;
   PEL_ASSERT( !file.empty() ) ;
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: writeMM_POST( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: print_items_PRE( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( os ) ;
   PEL_ASSERT( is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: create_replica_PRE( PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
//   PEL_CHECK_COLLECTIVE( is_distributed() ) ; les prototypes ne le sont pas
   PEL_ASSERT( is_a_prototype() ) ;
   PEL_ASSERT( exp != 0 ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Matrix:: create_replica_POST( LA_Matrix const* result,
                                 PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( result->name() == name() ) ;
   PEL_ASSERT( !result->is_a_prototype() ) ;
   PEL_ASSERT( IMPLIES( result->is_desynchronizable(),
                        !result->is_synchronized() ) ) ;
   PEL_ASSERT( result->distribution_strategy() != LA::InvalidDistribution ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
LA_Matrix:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
      PEL_ObjectRegister::create( PEL_Root::object(),
                                  "LA_Matrix descendant" ) ;
   return( result ) ;
}
