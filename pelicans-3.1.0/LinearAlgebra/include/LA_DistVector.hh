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

#ifndef LA_DIST_VECTOR_HH
#define LA_DIST_VECTOR_HH

#include <LA_Vector.hh>

#include <iosfwd>

#include <PEL_assertions.hh>
#include <PEL_Communicator.hh>
#include <PEL_DistributedPartition.hh>
#include <doubleVector.hh>

#include <LA_DistScatter.hh>
#include <LA_SeqVector.hh>
#include <vector>

class LA_DistScatter ;
class PEL_ModuleExplorer ;
class PEL_Vector ;
class size_t_vector ;

/*
   Sequential vectors of double, with a possible block structure in subvectors,
   compatibles with `LA_DistMatrix::' matrices.

   Local items are stored in contiguous vector and non-local
   ones are stored in a sparse structure.
*/

class PEL_EXPORT LA_DistVector : public LA_Vector
{
   public: //-----------------------------------------------------------

  //-- Instance delivery and initialization

      // Create and return an instance.
      static LA_DistVector* create( PEL_Object* a_owner,
                                    size_t a_nb_rows, size_t a_nb_local_rows,
                                    LA::DistributionStrategy dist_strat) ;

      virtual LA_DistVector* create_vector( PEL_Object* a_owner ) const ;

      virtual void re_initialize( size_t a_nb_rows,
                                  size_t a_nb_local_rows = PEL::bad_index() ) ;

   //-- Characteristics

      virtual LA_Implementation const* implementation( void ) const ;

   //-- Distributed processing

      virtual bool is_desynchronizable( void ) const ;

      virtual void synchronize( void ) ;

      virtual PEL_DistributedPartition const* row_distribution( void ) const ;

      virtual LA_DistScatter* create_scatter( PEL_Object* a_owner,
                                              size_t_vector const& from,
                                              size_t_vector const& to ) const ;

   //-- Access

      // `i'-th coefficient
      virtual double item( size_t i ) const ;

      virtual LA_SeqVector* create_local_vector( PEL_Object* a_owner ) const ;

   //-- Element change

      virtual void set_item( size_t i, double x ) ;

      virtual void add_to_item( size_t i, double x ) ;

   //-- BLAS level 1 : vector-vector operators

      virtual double dot( LA_Vector const* a ) const ;

      virtual double two_norm( void ) const ;

      virtual double max_norm( void ) const ;

      virtual void scale( double alpha ) ;

      virtual void set( LA_Vector const* a ) ;

      virtual void set( double value ) ;

      // Reinitialize by copying all coefficients of `a'.
      void set( LA_SeqVector const* a, size_t_vector const& local_to_global ) ;

      // Add all coefficients of `a'.
      void add( LA_SeqVector const* a, size_t_vector const& local_to_global ) ;

      virtual void set_as_v_product( LA_Vector const* a,
                                     LA_Vector const* b ) ;

      virtual void set_as_reciprocal( LA_Vector const* a,
                                      double smallest_inverted_item = 0.,
                                      double default_value = 1. ) ;

      virtual void sum( LA_Vector const* a, double alpha = 1. ) ;

      void recover_global_vector( LA_SeqVector* vec ) const ;

   //-- Input - Output

      virtual void print_items( std::ostream& os, size_t indent_width ) const ;

      virtual void write( std::string const& filename ) const ;

   //-- Hidden

      // pointer to internal data
      // (for time optimization, should be used with care)
      LA_SeqVector const* local_vector( void ) const ;
      LA_SeqVector* local_vector( void ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      LA_DistVector( void ) ;
     ~LA_DistVector( void ) ;
      LA_DistVector( LA_DistVector const& other ) ;
      LA_DistVector& operator=( LA_DistVector const& other ) ;

      LA_DistVector( PEL_Object* a_owner,
                     size_t a_nb_rows, size_t a_nb_local_rows,
                     LA::DistributionStrategy dist_strat ) ;

   //-- Distributed processing

      virtual void set_synchronized( void ) ;

      void send_subvector( size_t N, size_t i ) const ;

      void receive_subvector( size_t N, size_t i,
                              LA::SyncState effective_mode ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

      virtual bool implementation_POST(
                                LA_Implementation const* result ) const ;

   //-- Attributes

      PEL_DistributedPartition* DIST ;
      std::vector< intVector > GLOBAL_INDEXES ;
      std::vector< doubleVector > GLOBAL_VALUES ;
      LA_SeqVector* LOCAL_VECTOR ;
      bool VERBOSE ;
      bool INITIALIZED ;
      size_t FIRST ;
      size_t LAST ;
      size_t NB_RANKS ;
      bool HAS_BEEN_NULLIFIED ;
      size_t NB_EXTRA ;
} ;

#ifndef OUTLINE
   #include <LA_DistVector.icc>
#endif

#endif
