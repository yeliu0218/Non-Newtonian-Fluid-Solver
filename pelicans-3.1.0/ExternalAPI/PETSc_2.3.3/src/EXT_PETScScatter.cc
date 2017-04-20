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

#include <EXT_PETScScatter.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_DistributedPartition.hh>

#include <LA_SeqVector.hh>

#include <EXT_PETScImplementation.hh>
#include <EXT_PETScVector.hh>

//----------------------------------------------------------------------
EXT_PETScScatter*
EXT_PETScScatter:: create( PEL_Object* a_owner,
                           EXT_PETScVector const* global_vector,
                           size_t_vector const& a_repatriated_items_table,
                           size_t_vector const& a_local_indices_table ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScScatter:: create" ) ;
   PEL_CHECK_PRE(
      FORALL(
         ( size_t i=0 ; i<a_repatriated_items_table.size() ; ++i ),
         a_repatriated_items_table(i)<global_vector->nb_rows() ) ) ;
   PEL_CHECK_PRE(
      a_repatriated_items_table.size()==a_local_indices_table.size() ) ;
   
   EXT_PETScScatter* result = new EXT_PETScScatter( a_owner,
                                                    global_vector,
                                                    a_repatriated_items_table,
                                                    a_local_indices_table ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST(
      IMPLIES( global_vector->is_desynchronizable(),
               result->distribution() != 0 ) ) ;
   PEL_CHECK_POST(
      IMPLIES( result->distribution() != 0,
               result->distribution()->is_compatible(
                               global_vector->row_distribution() ) ) ) ;
   PEL_CHECK_POST( result->size() == a_repatriated_items_table.size() ) ;
   PEL_CHECK_POST( result->repatriated_items() ==
                                            a_repatriated_items_table ) ;
   PEL_CHECK_POST( result->local_indices() == a_local_indices_table ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
EXT_PETScScatter:: EXT_PETScScatter(
                           PEL_Object* a_owner,
                           EXT_PETScVector const* global_vector,
                           size_t_vector const& a_repatriated_items_table,
                           size_t_vector const& a_local_indices_table )
//----------------------------------------------------------------------
   : LA_Scatter( a_owner )
   , SIZE( a_repatriated_items_table.size() )
   , DIST( PEL_DistributedPartition::create( this ) )
   , NEEDED( a_repatriated_items_table )
   , LOCAL( a_local_indices_table )
{
   PEL_LABEL( "EXT_PETScScatter:: EXT_PETScScatter" ) ;

   DIST->set( global_vector->row_distribution() ) ;
   
   size_t s = 0 ;
   for( size_t i=0 ; i<SIZE ; ++i ) if( s<=LOCAL(i) ) s = LOCAL(i) ;
   
   int* idx_from = new int [ SIZE ] ;
   int* idx_to = new int [ SIZE ] ;
   for( size_t i=0 ; i<SIZE ; ++i )
   {
      idx_from[i] = NEEDED(i) ;
      idx_to[i] = LOCAL(i) ;
   }
   IS from, to ;
   PETSc_do( ISCreateGeneral( PETSC_COMM_SELF, SIZE, idx_from, &from ) ) ;
   PETSc_do( ISCreateGeneral( PETSC_COMM_SELF, SIZE, idx_to, &to ) ) ;
   
   PETSc_do( VecCreateSeq( PETSC_COMM_SELF, s+1, &SEQ ) ) ;
   PETSc_do( VecScatterCreate(
                global_vector->vector(), from, SEQ, to, &SCATTER ) ) ;
   
   delete [] idx_from ;
   delete [] idx_to ;
   PETSc_do( ISDestroy( from ) ) ;
   PETSc_do( ISDestroy( to ) ) ;
}

//----------------------------------------------------------------------
EXT_PETScScatter:: ~EXT_PETScScatter( void )
//----------------------------------------------------------------------
{
   PETSc_do( VecScatterDestroy( SCATTER ) ) ;
   PETSc_do( VecDestroy( SEQ ) ) ;
}

//----------------------------------------------------------------------
LA_Implementation const*
EXT_PETScScatter:: implementation( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScScatter:: implementation" ) ;
   
   static LA_Implementation const* result = EXT_PETScImplementation::object() ;

   PEL_CHECK_POST( implementation_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
EXT_PETScScatter:: size( void ) const
//----------------------------------------------------------------------
{
   return( SIZE ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
EXT_PETScScatter:: repatriated_items( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScScatter:: repatriated_items" ) ;
   
   size_t_vector const& result = NEEDED ;

   PEL_CHECK_POST( repatriated_items_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
EXT_PETScScatter:: local_indices( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScScatter:: local_indices" ) ;
   
   size_t_vector const& result = LOCAL ;

   PEL_CHECK_POST( local_indices_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_DistributedPartition const*
EXT_PETScScatter:: distribution( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScScatter:: distribution" ) ;

   PEL_DistributedPartition const* result = DIST ;

   PEL_CHECK_POST( distribution_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScScatter:: get( LA_Vector const* source,
                        LA_SeqVector* dest ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScScatter:: get" ) ;
   PEL_CHECK_PRE( get_PRE( source, dest) ) ;

   EXT_PETScVector const* psource =
                        static_cast<EXT_PETScVector const*>( source ) ;
   PEL_CHECK( dynamic_cast<EXT_PETScVector const*>( source ) != 0 ) ;

   PETSc_do( VecScatterBegin( SCATTER, psource->vector(),
                              SEQ, INSERT_VALUES, SCATTER_FORWARD ) ) ;
   PETSc_do( VecScatterEnd( SCATTER, psource->vector(),
                            SEQ, INSERT_VALUES, SCATTER_FORWARD ) ) ;
   PetscScalar *ptr_x ;

   PETSc_do( VecGetArray( SEQ, &ptr_x ) ) ;
   for( size_t i=0 ; i<SIZE ; i++ )
   {
      dest->set_item( LOCAL(i), ptr_x[LOCAL(i)] ) ;
   }
   PETSc_do( VecRestoreArray(SEQ,&ptr_x) ) ;
   dest->synchronize() ;
   
   PEL_CHECK_POST( get_POST( source, dest) ) ;
}

//----------------------------------------------------------------------
void
EXT_PETScScatter:: set( LA_SeqVector const* source,
                        LA_Vector* dest ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "EXT_PETScScatter:: set" ) ;
   PEL_CHECK_PRE( set_PRE( source, dest) ) ;

   PEL_CHECK( dynamic_cast<EXT_PETScVector const*>( dest ) != 0 ) ;
   EXT_PETScVector const* pdest = static_cast<EXT_PETScVector const*>( dest ) ;

   for( size_t i=0 ; i<SIZE ; i++ )
   {
      VecSetValue( SEQ, LOCAL(i), source->item(LOCAL(i)), INSERT_VALUES ) ;
   }
   PETSc_do( VecScatterBegin( SCATTER, SEQ, pdest->vector(),
                              INSERT_VALUES, SCATTER_REVERSE ) ) ;
   PETSc_do( VecScatterEnd( SCATTER, SEQ, pdest->vector(),
                            INSERT_VALUES, SCATTER_REVERSE ) ) ;
   static_cast<LA_Vector*>(dest)->synchronize() ;
   
   PEL_CHECK_POST( set_POST( source, dest) ) ;
}

//----------------------------------------------------------------------
bool
EXT_PETScScatter:: implementation_POST(
                                 LA_Implementation const* result ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( LA_Scatter::implementation_POST( result ) ) ;
   PEL_ASSERT( result == EXT_PETScImplementation::object() ) ;
   return( true ) ;
}
