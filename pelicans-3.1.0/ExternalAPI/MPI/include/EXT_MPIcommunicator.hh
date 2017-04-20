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

#ifndef EXT_MPI_COMMUNICATOR_HH
#define EXT_MPI_COMMUNICATOR_HH

#include <PEL_Communicator.hh>

/*
MPI communicator.

PUBLISHED
*/

class EXT_MPIcommunicator : public PEL_Communicator
{
   public: //-----------------------------------------------------------

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      EXT_MPIcommunicator( void ) ;
     ~EXT_MPIcommunicator( void ) ;
      EXT_MPIcommunicator( EXT_MPIcommunicator const& other ) ;
      EXT_MPIcommunicator& operator=( EXT_MPIcommunicator const& other ) ;
      
      enum { TAG_INT, TAG_DOUBLE, TAG_CHAR } ;

   //-- Characteristics
         
      virtual size_t nb_ranks( void ) const ;
         
      virtual size_t rank( void ) const ;
         
   //-- Point-to-point blocking communication

      virtual void send( size_t dest, int const* value, int nb ) const ;
      virtual void receive( size_t src, int* value, int nb ) const ;
      
      virtual void send( size_t dest, double const* value, int nb ) const ;
      virtual void receive( size_t src, double* value, int nb ) const ;

      virtual void send( size_t dest, char const* value, int nb ) const ;
      virtual void receive( size_t src, char* value, int nb ) const ;
      
   //-- Point-to-point non-blocking communication
      
      virtual void* Isend( size_t dest, int const* value, int nb ) const ;
      virtual void* Ireceive( size_t src, int* value, int nb ) const ;
      
      virtual void* Isend( size_t dest, double const* value, int nb ) const ;
      virtual void* Ireceive( size_t src, double* value, int nb ) const ;
      
   //-- Collective communication: broadcast
      
      virtual void broadcast( int* value, size_t nb, size_t root ) const ;

      virtual void broadcast( double* value, size_t nb, size_t root ) const ;
      
      virtual void broadcast( char* value, size_t nb, size_t root ) const ;
      
   //-- Collective communication: gather
      
      virtual void gather( double const* value, size_t nb, 
                           double* result, size_t root ) const ;
      
      virtual void gather( int const* value, size_t nb, 
                           int* result, size_t root ) const ;
      
   //-- Collective communication: gather to all
      
      virtual void all_gather( int const* value, size_t nb, 
                               int* result ) const ;
      
      virtual void all_gather_v( double const* values,
                                 size_t nb,
                                 double* result,
                                 intVector const& partition,
                                 intVector const& start ) const ;

   //-- Collective communication: all to all scatter/gather
      
      virtual void all_to_all( int const* value, size_t nb, 
                               int* result ) const ;
      
   //-- Collective communication: reductions
         
      virtual bool boolean_and( bool value ) const ;

      virtual bool boolean_or( bool value ) const ;
      
      virtual double sum( double value ) const ;

      virtual void sum_vector( double* values, int nb ) const ;
      
      virtual double min( double value ) const ;

      virtual double max( double value ) const ;
      
      virtual bool same_value_everywhere( double val ) const ;

      virtual bool same_value_everywhere( int val ) const ;
      
   //-- Collective communication: synchronization
      
      virtual void barrier( void ) const ;
      
      virtual void wait( void* request ) const ;
      
   //-- Class attributes
      
      static EXT_MPIcommunicator* SINGLETON ;

   //-- Attributes
      
      mutable int count ;
} ;

#endif
