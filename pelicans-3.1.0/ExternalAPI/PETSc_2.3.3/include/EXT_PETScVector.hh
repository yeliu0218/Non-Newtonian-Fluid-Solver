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

#ifndef EXT_PETSC_VECTOR_HH
#define EXT_PETSC_VECTOR_HH

#include <PEL_Object.hh>

#include <iosfwd>

#include <PEL_assertions.hh>
#include <PEL_Communicator.hh>
#include <PEL_DistributedPartition.hh>
#include <doubleVector.hh>

#include <LA_SeqVector.hh>

#include <EXT_PETScScatter.hh>
#include <EXT_PETScAPI.hh>

class size_t_vector ;

/* PETSc vectors.
*/

class EXT_PETScVector : public LA_Vector
{
   public: //-----------------------------------------------------------

  //-- Instance delivery and initialization

      // Create and return an instance.
      static EXT_PETScVector* create( PEL_Object* a_owner,
                                      bool sequential,
                                      size_t a_nb_rows ) ;

      virtual EXT_PETScVector* create_vector( PEL_Object* a_owner ) const ;

      virtual void re_initialize( size_t a_nb_rows,
                                  size_t a_nb_local_rows = PEL::bad_index() ) ;

   //-- Characteristics

      virtual LA_Implementation const* implementation( void ) const ;

   //-- Distributed processing

      virtual bool is_desynchronizable( void ) const ;

      virtual void start_local_modifs( void ) ;

      virtual void stop_local_modifs( void ) ;

      virtual void synchronize( void ) ;

      virtual PEL_DistributedPartition const* row_distribution( void ) const ;

      virtual EXT_PETScScatter* create_scatter( PEL_Object* a_owner,
                                                size_t_vector const& from,
                                                size_t_vector const& to ) const ;

   //-- Access

      virtual LA_SeqVector* create_local_vector( PEL_Object* a_owner ) const ;

      // `i'-th coefficient
      virtual double item( size_t i ) const ;

   //-- Element change

      virtual void set_item( size_t i, double x ) ;

      virtual void add_to_item( size_t i, double x ) ;

      void set_verbosity( bool verbosity ) ;

   //-- BLAS level 1 : vector-vector operators

      virtual double dot( LA_Vector const* a ) const ;

      virtual double two_norm( void ) const ;

      virtual double max_norm( void ) const ;

      virtual void scale( double alpha ) ;

      virtual void set( LA_Vector const* a ) ;

      virtual void set( double value ) ;

      virtual void set_as_v_product( LA_Vector const* a,
                                     LA_Vector const* b ) ;

      virtual void set_as_reciprocal( LA_Vector const* a,
                                      double smallest_inverted_item=0.0,
                                      double default_value=1.0 ) ;

      virtual void sum( LA_Vector const* a, double alpha=1.0 ) ;

      LA_SeqVector const* local_vector( void ) const ;

   //-- Input - Output

      virtual void print_items( std::ostream& os, size_t indent_width ) const ;

      virtual void write( std::string const& filename ) const ;

   //-- Hidden

      // pointer to internal data
      // (for time optimization, should be used with care)
      Vec& vector( void ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      EXT_PETScVector( void ) ;
     ~EXT_PETScVector( void ) ;
      EXT_PETScVector( EXT_PETScVector const& other ) ;
      EXT_PETScVector& operator=(
                            EXT_PETScVector const& other ) ;

      EXT_PETScVector( PEL_Object* a_owner, bool sequential,
                       size_t a_nb_rows ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

      virtual bool implementation_POST(
                                LA_Implementation const* result ) const ;

   //-- Attributes

      Vec VECTOR ;
      PEL_DistributedPartition* DIST ;
      bool SEQ ;
      size_t FIRST ;
      size_t LAST ;
} ;


#endif
