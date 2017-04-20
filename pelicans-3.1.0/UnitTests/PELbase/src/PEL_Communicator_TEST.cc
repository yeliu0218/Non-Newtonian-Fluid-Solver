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

#include <PEL_Communicator_TEST.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Communicator.hh>
#include <PEL_DoubleComparatorExact.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <boolArray2D.hh>
#include <boolVector.hh>
#include <doubleArray2D.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <intArray2D.hh>
#include <size_t_array2D.hh>
#include <size_t_vector.hh>
#include <stringVector.hh>

#include <string>
#include <iostream>

using std::endl ;

PEL_Communicator_TEST*
PEL_Communicator_TEST::REGISTRATOR = new PEL_Communicator_TEST() ;

//-------------------------------------------------------------------------
PEL_Communicator_TEST:: PEL_Communicator_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "PEL_Communicator", "PEL_Communicator_TEST" )
{
}

//-------------------------------------------------------------------------
PEL_Communicator_TEST:: ~PEL_Communicator_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
PEL_Communicator_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: process_one_test" ) ;

   PEL_Communicator const* com = PEL_Exec::communicator() ;

   out() << "| Communicator is " << com->name() << std::endl ;

   test_send_recv_size_t( com ) ;
   test_send_recv_size_t_vector( com ) ;
   test_send_recv_size_t_array2D( com ) ;
   test_send_recv_int( com ) ;
   test_send_recv_intVector( com ) ;
   test_send_recv_intArray2D( com ) ;
   test_send_recv_double( com ) ;
   test_send_recv_doubleVector( com ) ;
   test_send_recv_doubleArray2D( com ) ;
   test_send_recv_bool( com ) ;
   test_send_recv_boolVector( com ) ;
   test_send_recv_boolArray2D( com ) ;
   test_send_recv_string( com ) ;
   test_send_recv_stringVector( com ) ;
   test_send_recv_intPtr( com ) ;
   test_send_recv_doublePtr( com ) ;
   test_send_recv_charPtr( com ) ;

   test_NB_send_recv_intPtr( com ) ;
   test_NB_send_recv_doublePtr( com ) ;

   test_broadcast( com ) ;

   test_gather( com ) ;

   test_all_gather( com ) ;

   test_all_to_all( com ) ;

   test_boolean_and_or( com ) ;

   test_sum( com ) ;

   test_max_min( com ) ;

   test_same_value( com ) ;

   test_merge( com, exp ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_send_recv_size_t( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_send_recv_size_t" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   size_t val0 = 10 ;

   size_t val = PEL::bad_index() ;
   if( rank == 0 )
   {
      val = val0 ;
      com->send( next_proc, val ) ;
      com->receive( prev_proc, val ) ;
   }
   else
   {
      com->receive( prev_proc, val ) ;
      com->send( next_proc, val+1 ) ;
   }
   bool ok = ( rank==0 ? val==val0+nbrk-1 : val==val0+rank-1 ) ;

   notify_one_test_result( "send_receive(size_t)", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_send_recv_size_t_vector( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_send_recv_size_t_vector" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   size_t_vector val0( 4 ) ;
   for( size_t i=0 ; i<val0.size() ; ++i )
   {
      val0( i ) = 2*i ;
   }

   size_t_vector val( 0 ) ;
   if( rank == 0 )
   {
      val = val0 ;
      com->send( next_proc, val ) ;
      com->receive( prev_proc, val ) ;
   }
   else
   {
      com->receive( prev_proc, val ) ;
      size_t_vector valp1 = val ;
      for( size_t i=0 ; i<valp1.size() ; ++i )
      {
         ++valp1( i ) ;
      }
      com->send( next_proc, valp1 ) ;
   }
   bool ok = ( val.size() == val0.size() ) ;
   for( size_t i=0 ; i<val.size() ; ++i )
   {
      ok &= ( rank==0 ? val(i)==val0(i)+nbrk-1 : val(i)==val0(i)+rank-1 ) ;
   }

   notify_one_test_result( "send_receive(size_t_vector)", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_send_recv_size_t_array2D( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_send_recv_size_t_array2D" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   size_t_array2D val0( 4, 3 ) ;
   for( size_t i=0 ; i<val0.index_bound( 0 ) ; ++i )
   {
      for( size_t j=0 ; j<val0.index_bound( 1 ) ; ++j )
      {
         val0( i, j ) = 2*i + 8*j ;
      }
   }

   size_t_array2D val( 0, 0 ) ;
   if( rank == 0 )
   {
      val = val0 ;
      com->send( next_proc, val ) ;
      com->receive( prev_proc, val ) ;
   }
   else
   {
      com->receive( prev_proc, val ) ;
      size_t_array2D valp1 = val ;
      for( size_t i=0 ; i<valp1.index_bound( 0 ) ; ++i )
      {
         for( size_t j=0 ; j<valp1.index_bound( 1 ) ; ++j )
         {
            ++valp1( i, j ) ;
         }
      }
      com->send( next_proc, valp1 ) ;
   }
   bool ok = ( val.index_bound( 0 ) == val0.index_bound( 0 ) &&
               val.index_bound( 1 ) == val0.index_bound( 1 ) ) ;
   for( size_t i=0 ; i<val.index_bound( 0 ) ; ++i )
   {
      for( size_t j=0 ; j<val.index_bound( 1 ) ; ++j )
      {
         ok &= ( rank == 0 ?
                    val( i, j ) == val0( i, j ) + nbrk-1
               : val( i, j ) == val0( i, j ) + rank-1 ) ;
      }
   }

   notify_one_test_result( "send_receive(size_t_array2D)", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_send_recv_int( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_send_recv_int" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   int val0 = -10 ;

   int val = PEL::bad_int() ;
   if( rank == 0 )
   {
      val = val0 ;
      com->send( next_proc, val ) ;
      com->receive( prev_proc, val ) ;
   }
   else
   {
      com->receive( prev_proc, val ) ;
      com->send( next_proc, val+1 ) ;
   }
   bool ok = ( rank == 0 ?
                  val == val0 + (int)nbrk-1
               : val == val0 + (int)rank-1 ) ;

   notify_one_test_result( "send_receive(int)", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_send_recv_intVector( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_send_recv_intVector" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   intVector val0( 4 ) ;
   for( size_t i=0 ; i<val0.size() ; ++i )
   {
      val0( i ) = -2*i ;
   }

   intVector val( 0 ) ;
   if( rank == 0 )
   {
      val = val0 ;
      com->send( next_proc, val ) ;
      com->receive( prev_proc, val ) ;
   }
   else
   {
      com->receive( prev_proc, val ) ;
      intVector valp1 = val ;
      for( size_t i=0 ; i<valp1.size() ; ++i )
      {
         ++valp1( i ) ;
      }
      com->send( next_proc, valp1 ) ;
   }
   bool ok = ( val.size() == val0.size() ) ;
   for( size_t i=0 ; i<val.size() ; ++i )
   {
      ok &= ( rank == 0 ?
                 val( i ) == val0( i ) + (int)nbrk-1
              : val( i ) == val0( i ) + (int)rank-1 ) ;
   }

   notify_one_test_result( "send_receive(intVector)", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_send_recv_intArray2D(
                                               PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_send_recv_intArray2D" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   intArray2D val0( 4, 3 ) ;
   for( size_t i=0 ; i<val0.index_bound( 0 ) ; ++i )
   {
      for( size_t j=0 ; j<val0.index_bound( 1 ) ; ++j )
      {
         val0( i, j ) = 2*i - 8*j ;
      }
   }

   intArray2D val( 0, 0 ) ;
   if( rank == 0 )
   {
      val = val0 ;
      com->send( next_proc, val ) ;
      com->receive( prev_proc, val ) ;
   }
   else
   {
      com->receive( prev_proc, val ) ;
      intArray2D valp1 = val ;
      for( size_t i=0 ; i<valp1.index_bound( 0 ) ; ++i )
      {
         for( size_t j=0 ; j<valp1.index_bound( 1 ) ; ++j )
         {
            ++valp1( i, j ) ;
         }
      }
      com->send( next_proc, valp1 ) ;
   }
   bool ok = ( val.index_bound( 0 ) == val0.index_bound( 0 ) &&
               val.index_bound( 1 ) == val0.index_bound( 1 ) ) ;
   for( size_t i=0 ; i<val.index_bound( 0 ) ; ++i )
   {
      for( size_t j=0 ; j<val.index_bound( 1 ) ; ++j )
      {
         ok &= ( rank == 0 ?
                    val( i, j ) == val0( i, j ) + (int)nbrk-1
               : val( i, j ) == val0( i, j ) + (int)rank-1 ) ;
      }
   }

   notify_one_test_result( "send_receive(intArray2D)", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_send_recv_double( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_send_recv_double" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   double val0 = 10.0 ;

   double val = PEL::bad_double() ;
   if( rank == 0 )
   {
      val = val0 ;
      com->send( next_proc, val ) ;
      com->receive( prev_proc, val ) ;
   }
   else
   {
      com->receive( prev_proc, val ) ;
      com->send( next_proc, val+1.0 ) ;
   }
   bool ok = ( rank == 0 ?
                  PEL::abs( 1.0 - val/(val0+nbrk-1.0) ) < 1.e-10
             : PEL::abs( 1.0 - val/(val0+rank-1.0) ) < 1.e-10  ) ;

   notify_one_test_result( "send_receive(double)", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_send_recv_doubleVector(
                                             PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_send_recv_doubleVector" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   doubleVector val0( 4 ) ;
   for( size_t i=0 ; i<val0.size() ; ++i )
   {
      val0( i ) = 2.0*i + 1.0 ;
   }

   doubleVector val( 0 ) ;
   if( rank == 0 )
   {
      val = val0 ;
      com->send( next_proc, val ) ;
      com->receive( prev_proc, val ) ;
   }
   else
   {
      com->receive( prev_proc, val ) ;
      doubleVector valp1 = val ;
      for( size_t i=0 ; i<valp1.size() ; ++i )
      {
         valp1( i ) += 1.0 ;
      }
      com->send( next_proc, valp1 ) ;
   }
   bool ok = ( val.size() == val0.size() ) ;
   for( size_t i=0 ; i<val.size() ; ++i )
   {
      ok &= ( rank == 0 ?
                 PEL::abs( 1.0 - val(i)/(val0(i)+nbrk-1.0) ) < 1.e-10
            : PEL::abs( 1.0 - val(i)/(val0(i)+rank-1.0) ) < 1.e-10 ) ;
   }

   notify_one_test_result( "send_receive(doubleVector)", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_send_recv_doubleArray2D(
                                             PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_send_recv_doubleArray2D" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   doubleArray2D val0( 4, 3 ) ;
   for( size_t i=0 ; i<val0.index_bound( 0 ) ; ++i )
   {
      for( size_t j=0 ; j<val0.index_bound( 1 ) ; ++j )
      {
         val0( i, j ) = 2.0*i + 8.0*j + 1.0 ;
      }
   }

   doubleArray2D val( 0, 0 ) ;
   if( rank == 0 )
   {
      val = val0 ;
      com->send( next_proc, val ) ;
      com->receive( prev_proc, val ) ;
   }
   else
   {
      com->receive( prev_proc, val ) ;
      doubleArray2D valp1 = val ;
      for( size_t i=0 ; i<valp1.index_bound( 0 ) ; ++i )
      {
         for( size_t j=0 ; j<valp1.index_bound( 1 ) ; ++j )
         {
            valp1( i, j ) += 1.0 ;
         }
      }
      com->send( next_proc, valp1 ) ;
   }
   bool ok = ( val.index_bound( 0 ) == val0.index_bound( 0 ) &&
               val.index_bound( 1 ) == val0.index_bound( 1 ) ) ;
   for( size_t i=0 ; i<val.index_bound( 0 ) ; ++i )
   {
      for( size_t j=0 ; j<val.index_bound( 1 ) ; ++j )
      {
         ok &= ( rank == 0 ?
                    PEL::abs( 1.0 - val(i,j)/(val0(i,j)+nbrk-1.0) ) < 1.e-10
               : PEL::abs( 1.0 - val(i,j)/(val0(i,j)+rank-1.0) ) < 1.e-10 ) ;
      }
   }

   notify_one_test_result( "send_receive(doubleArray2D)", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_send_recv_bool( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_send_recv_bool" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   bool val0 = true ;

   bool val ;
   if( rank == 0 )
   {
      val = val0 ;
      com->send( next_proc, val ) ;
      com->receive( prev_proc, val ) ;
   }
   else
   {
      com->receive( prev_proc, val ) ;
      com->send( next_proc, !val ) ;
   }
   bool ok = ( rank==0 ? val==((int)nbrk%2==1) : val==((int)rank%2==1) ) ;

   notify_one_test_result( "send_receive(bool)", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_send_recv_boolVector( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_send_recv_boolVector" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   boolVector val0( 4 ) ;
   val0.set( true ) ;

   boolVector val( 0 ) ;
   if( rank == 0 )
   {
      val = val0 ;
      com->send( next_proc, val ) ;
      com->receive( prev_proc, val ) ;
   }
   else
   {
      com->receive( prev_proc, val ) ;
      boolVector nval = val ;
      for( size_t i=0 ; i<nval.size() ; ++i )
      {
         nval( i ) = !nval( i ) ;
      }
      com->send( next_proc, nval ) ;
   }
   bool ok = ( val.size() == val0.size() ) ;
   for( size_t i=0 ; i<val.size() ; ++i )
   {
      ok &= ( rank == 0 ?
                 val( i ) == ( (int)nbrk%2 == 1 )
            : val( i ) == ( (int)rank%2 == 1 ) ) ;
   }

   notify_one_test_result( "send_receive(boolVector)", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_send_recv_boolArray2D(
                                           PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_send_recv_boolArray2D" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   boolArray2D val0( 4, 3 ) ;
   val0.set( true ) ;

   boolArray2D val( 0, 0 ) ;
   if( rank == 0 )
   {
      val = val0 ;
      com->send( next_proc, val ) ;
      com->receive( prev_proc, val ) ;
   }
   else
   {
      com->receive( prev_proc, val ) ;
      boolArray2D nval = val ;
      for( size_t i=0 ; i<nval.index_bound( 0 ) ; ++i )
      {
         for( size_t j=0 ; j<nval.index_bound( 1 ) ; ++j )
         {
            nval( i, j ) = !nval( i, j ) ;
         }
      }
      com->send( next_proc, nval ) ;
   }
   bool ok = ( val.index_bound( 0 ) == val0.index_bound( 0 ) &&
               val.index_bound( 1 ) == val0.index_bound( 1 ) ) ;
   for( size_t i=0 ; i<val.index_bound( 0 ) ; ++i )
   {
      for( size_t j=0 ; j<val.index_bound( 1 ) ; ++j )
      {
         ok &= ( rank == 0 ?
                    val( i, j ) == ( (int)nbrk%2 == 1 )
               : val( i, j ) == ( (int)rank%2 == 1 ) ) ;
      }
   }

   notify_one_test_result( "send_receive(boolArray2D)", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_send_recv_string( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_send_recv_string" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   std::string val ;
   if( rank == 0 )
   {
      val = "" ;
      com->send( next_proc, val ) ;
      com->receive( prev_proc, val ) ;
   }
   else
   {
      com->receive( prev_proc, val ) ;
      val += "*" ;
      com->send( next_proc, val ) ;
   }
   std::string const s = std::string( (rank == 0 ? nbrk-1 : rank), '*' ) ;
   bool ok = ( val == s ) ;

   notify_one_test_result( "send_receive(string)", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_send_recv_stringVector( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_send_recv_stringVector" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   stringVector const val0 = "Hello, ,World" ;

   stringVector val(0) ;
   if( rank == 0 )
   {
      val = val0 ;
      com->send( next_proc, val ) ;
      com->receive( prev_proc, val ) ;
   }
   else
   {
      com->receive( prev_proc, val ) ;
      com->send( next_proc, val ) ;
   }
   bool ok = ( val == val0 ) ;

   notify_one_test_result( "send_receive(stringVector)", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_send_recv_intPtr( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_send_recv_intPtr" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   int const nb_vals = 4 ;
   int* val0 = new int[ nb_vals ] ;
   for( int i=0 ; i<nb_vals ; ++i )
   {
      val0[i] = -2*i ;
   }

   int* val = new int[ nb_vals ] ;
   if( rank == 0 )
   {
      for( int i=0 ; i<nb_vals ; ++i )
      {
         val[i] = val0[i] ;
      }
      com->send( next_proc, val, nb_vals ) ;
      com->receive( prev_proc, val, nb_vals ) ;
   }
   else
   {
      com->receive( prev_proc, val, nb_vals ) ;
      int* valp1 = new int[ nb_vals ] ;
      for( int i=0 ; i<nb_vals ; ++i )
      {
         valp1[ i ] = val[ i ] + 1 ;
      }
      com->send( next_proc, valp1, nb_vals ) ;
      delete [] valp1 ;
   }

   bool ok = true ;
   for( int i=0 ; i<nb_vals ; ++i )
   {
      ok &= ( rank == 0 ?
                 val[i] == val0[i] + (int)nbrk-1
            : val[i] == val0[i] + (int)rank-1 ) ;
   }
   notify_one_test_result( "send_receive(int*)", ok ) ;

   delete [] val  ;
   delete [] val0 ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_send_recv_doublePtr( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_send_recv_doublePtr" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   int const nb_vals = 4 ;
   double* val0 = new double[ nb_vals ] ;
   for( int i=0 ; i<nb_vals ; ++i )
   {
      val0[i] = 2.0*i + 1.0 ;
   }

   double* val = new double[ nb_vals ] ;
   if( rank == 0 )
   {
      for( int i=0 ; i<nb_vals ; ++i )
      {
         val[i] = val0[i] ;
      }
      com->send( next_proc, val, nb_vals ) ;
      com->receive( prev_proc, val, nb_vals ) ;
   }
   else
   {
      com->receive( prev_proc, val, nb_vals ) ;
      double* valp1 = new double[ nb_vals ] ;
      for( int i=0 ; i<nb_vals ; ++i )
      {
         valp1[ i ] = val[ i ] + 1.0 ;
      }
      com->send( next_proc, valp1, nb_vals ) ;
      delete [] valp1 ;
   }

   bool ok = true ;
   for( int i=0 ; i<nb_vals ; ++i )
   {
      ok &= ( rank == 0 ?
                 PEL::abs( 1.0 - val[i]/(val0[i]+nbrk-1.0) ) < 1.e-10
            : PEL::abs( 1.0 - val[i]/(val0[i]+rank-1.0) ) < 1.e-10 ) ;
   }
   notify_one_test_result( "send_receive(double*)", ok ) ;

   delete [] val  ;
   delete [] val0 ;
}


//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_send_recv_charPtr( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_send_recv_charPtr" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   int const nb_vals = 4 ;
   char* val0 = new char[ nb_vals ] ;
   for( int i=0 ; i<nb_vals ; ++i )
   {
      val0[i] = (char) 2*i ;
   }

   char* val = new char[ nb_vals ] ;
   if( rank == 0 )
   {
      for( int i=0 ; i<nb_vals ; ++i )
      {
         val[i] = val0[i] ;
      }
      com->send( next_proc, val, nb_vals ) ;
      com->receive( prev_proc, val, nb_vals ) ;
   }
   else
   {
      com->receive( prev_proc, val, nb_vals ) ;
      char* valp1 = new char[ nb_vals ] ;
      for( int i=0 ; i<nb_vals ; ++i )
      {
         valp1[ i ] = val[ i ] + 1 ;
      }
      com->send( next_proc, valp1, nb_vals ) ;
      delete [] valp1 ;
   }

   bool ok = true ;
   for( int i=0 ; i<nb_vals ; ++i )
   {
      ok &= ( rank == 0 ? val[i] == (char) (val0[i]+nbrk-1)
                        : val[i] == (char) (val0[i]+rank-1) ) ;
   }
   notify_one_test_result( "send_receive(char*)", ok ) ;

   delete [] val  ;
   delete [] val0 ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_NB_send_recv_intPtr( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_NB_send_recv_intPtr" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   int const nb_vals = 4 ;
   int* val0 = new int[ nb_vals ] ;
   for( int i=0 ; i<nb_vals ; ++i )
   {
      val0[i] = -2*i ;
   }

   int* val = new int[ nb_vals ] ;
   if( rank == 0 )
   {
      for( int i=0 ; i<nb_vals ; ++i )
      {
         val[i] = val0[i] ;
      }
      com->Isend( next_proc, val, nb_vals ) ;
      void* r_request = com->Ireceive( prev_proc, val, nb_vals ) ;
      com->wait( r_request ) ;
   }
   else
   {
      void* r_request = com->Ireceive( prev_proc, val, nb_vals ) ;
      com->wait( r_request ) ;
      int* valp1 = new int[ nb_vals ] ;
      for( int i=0 ; i<nb_vals ; ++i )
      {
         valp1[ i ] = val[ i ] + 1 ;
      }
      void* s_request = com->Isend( next_proc, valp1, nb_vals ) ;
      com->wait( s_request ) ;
      delete [] valp1 ;
   }

   bool ok = true ;
   for( int i=0 ; i<nb_vals ; ++i )
   {
      ok &= ( rank == 0 ? val[i] == val0[i] + (int)nbrk-1
                        : val[i] == val0[i] + (int)rank-1 ) ;
   }
   notify_one_test_result( "Isend_Ireceive(int*)", ok ) ;

   delete [] val  ;
   delete [] val0 ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_NB_send_recv_doublePtr( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_NB_send_recv_doublePtr" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   if( nbrk == 1 ) return ;
   //              ------

   size_t next_proc = ( rank < (nbrk-1) ? rank+1 : 0  ) ;
   size_t prev_proc = ( rank != 0 ? rank-1 : nbrk-1 ) ;

   int const nb_vals = 4 ;
   double* val0 = new double[ nb_vals ] ;
   for( int i=0 ; i<nb_vals ; ++i )
   {
      val0[i] = 2.0*i + 1.0 ;
   }

   double* val = new double[ nb_vals ] ;
   if( rank == 0 )
   {
      for( int i=0 ; i<nb_vals ; ++i )
      {
         val[i] = val0[i] ;
      }
      com->Isend( next_proc, val, nb_vals ) ;
      void* r_request = com->Ireceive( prev_proc, val, nb_vals ) ;
      com->wait( r_request ) ;
   }
   else
   {
      void* r_request = com->Ireceive( prev_proc, val, nb_vals ) ;
      com->wait( r_request ) ;
      double* valp1 = new double[ nb_vals ] ;
      for( int i=0 ; i<nb_vals ; ++i )
      {
         valp1[ i ] = val[ i ] + 1.0 ;
      }
      void* s_request = com->Isend( next_proc, valp1, nb_vals ) ;
      com->wait( s_request ) ;
      delete [] valp1 ;
   }

   bool ok = true ;
   for( int i=0 ; i<nb_vals ; ++i )
   {
      ok &= ( rank == 0 ?
                 PEL::abs( 1.0 - val[i]/(val0[i]+nbrk-1.0) ) < 1.e-10
            : PEL::abs( 1.0 - val[i]/(val0[i]+rank-1.0) ) < 1.e-10 ) ;
   }
   notify_one_test_result( "Isend_Ireceive(double*)", ok ) ;

   delete [] val  ;
   delete [] val0 ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_broadcast( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_broadcast" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   size_t root = nbrk-1 ;

   // ---------------
   // broadcast(size_t)
   {
      size_t val = rank+10 ;
      com->broadcast( val, root ) ;
      bool ok = ( val == root+10 ) ;
      notify_one_test_result( "broadcast(size_t)", ok ) ;
   }

   // ---------------
   // broadcast(size_t_vector)
   {
      size_t const nb_vals = 10 ;
      size_t_vector vals(0) ;
      if( rank == root )
      {
         vals.re_initialize( nb_vals ) ;
         for( size_t i=0 ; i<nb_vals ; ++i )
         {
            vals(i) = rank+10*i ;
         }
      }
      com->broadcast( vals, root ) ;
      bool ok = ( vals.size() == nb_vals ) ;
      for( size_t i=0 ; ok && i<nb_vals ; ++i )
      {
         ok &= ( vals(i) == root+10*i ) ;
      }
      notify_one_test_result( "broadcast(size_t_vector)", ok ) ;
   }

   // ---------------
   // broadcast(size_t_array2D)
   {
      size_t const nb_vals = 10 ;
      size_t_array2D vals(0,0) ;
      if( rank == root )
      {
         vals.re_initialize( nb_vals, nb_vals ) ;
         for( size_t i=0 ; i<nb_vals ; ++i )
         {
            for( size_t j=0 ; j<nb_vals ; ++j )
            {
               vals(i,j) = rank+10*i+20*j ;
            }
         }
      }
      com->broadcast( vals, root ) ;
      bool ok = ( vals.index_bound(0) == nb_vals &&
                  vals.index_bound(1) == nb_vals ) ;
      for( size_t i=0 ; ok && i<nb_vals ; ++i )
      {
         for( size_t j=0 ; ok && j<nb_vals ; ++j )
         {
            ok &= ( vals(i,j) == root+10*i+20*j ) ;
         }
      }
      notify_one_test_result( "broadcast(size_t_array2D)", ok ) ;
   }

   // ---------------
   // broadcast(int)
   {
      int val = rank+10 ;
      com->broadcast( val, root ) ;
      bool ok = ( val == int( root+10 ) ) ;
      notify_one_test_result( "broadcast(int)", ok ) ;
   }

   // ---------------
   // broadcast(intVector)
   {
      size_t const nb_vals = 10 ;
      intVector vals(0) ;
      if( rank == root )
      {
         vals.re_initialize( nb_vals ) ;
         for( size_t i=0 ; i<nb_vals ; ++i )
         {
            vals(i) = int( rank+10*i ) ;
         }
      }
      com->broadcast( vals, root ) ;
      bool ok = ( vals.size() == nb_vals ) ;
      for( size_t i=0 ; ok && i<nb_vals ; ++i )
      {
         ok &= ( vals(i) == int( root+10*i ) ) ;
      }
      notify_one_test_result( "broadcast(intVector)", ok ) ;
   }

   // ---------------
   // broadcast(intArray2D)
   {
      size_t const nb_vals = 10 ;
      intArray2D vals(0,0) ;
      if( rank == root )
      {
         vals.re_initialize( nb_vals, nb_vals ) ;
         for( size_t i=0 ; i<nb_vals ; ++i )
         {
            for( size_t j=0 ; j<nb_vals ; ++j )
            {
               vals(i,j) = int( rank+10*i-20*j ) ;
            }
         }
      }
      com->broadcast( vals, root ) ;
      bool ok = ( vals.index_bound(0) == nb_vals &&
                  vals.index_bound(1) == nb_vals ) ;
      for( size_t i=0 ; ok && i<nb_vals ; ++i )
      {
         for( size_t j=0 ; ok && j<nb_vals ; ++j )
         {
            ok &= ( vals(i,j) == int( root+10*i-20*j ) ) ;
         }
      }
      notify_one_test_result( "broadcast(intArray2D)", ok ) ;
   }

   // ---------------
   // broadcast(bool)
   {
      bool val = ( rank == root ) ;
      com->broadcast( val, root ) ;
      bool ok = ( val == true ) ;
      notify_one_test_result( "broadcast(bool)", ok ) ;
   }

   // ---------------
   // broadcast(boolVector)
   {
      size_t const nb_vals = 10 ;
      boolVector vals(0) ;
      if( rank == root )
      {
         vals.re_initialize( nb_vals ) ;
         for( size_t i=0 ; i<nb_vals ; ++i )
         {
            vals(i) = ( i%2 == 0 ) ;
         }
      }
      com->broadcast( vals, root ) ;
      bool ok = ( vals.size() == nb_vals ) ;
      for( size_t i=0 ; ok && i<nb_vals ; ++i )
      {
         ok &= ( vals(i) == ( i%2 == 0 ) ) ;
      }
      notify_one_test_result( "broadcast(boolVector)", ok ) ;
   }

   // ---------------
   // broadcast(boolArray2D)
   {
      size_t const nb_vals = 10 ;
      boolArray2D vals(0,0) ;
      if( rank == root )
      {
         vals.re_initialize( nb_vals, nb_vals ) ;
         for( size_t i=0 ; i<nb_vals ; ++i )
         {
            for( size_t j=0 ; j<nb_vals ; ++j )
            {
               vals(i,j) = ( i%2 == 0 && j%3 == 0 ) ;
            }
         }
      }
      com->broadcast( vals, root ) ;
      bool ok = ( vals.index_bound(0) == nb_vals &&
                  vals.index_bound(1) == nb_vals ) ;
      for( size_t i=0 ; ok && i<nb_vals ; ++i )
      {
         for( size_t j=0 ; ok && j<nb_vals ; ++j )
         {
            ok &= ( vals(i,j) == ( i%2 == 0 && j%3 == 0 ) ) ;
         }
      }
      notify_one_test_result( "broadcast(boolArray2D)", ok ) ;
   }

   // ---------------
   // broadcast(double)
   {
      double val = rank+10. ;
      com->broadcast( val, root ) ;
      bool ok = ( val == root+10. ) ;
      notify_one_test_result( "broadcast(double)", ok ) ;
   }

   // ---------------
   // broadcast(doubleVector)
   {
      size_t const nb_vals = 10 ;
      doubleVector vals( nb_vals ) ;
      for( size_t i=0 ; i<nb_vals ; ++i )
      {
         vals( i ) = double( rank+10*i ) ;
      }
      com->broadcast( vals, root ) ;
      bool ok = ( vals.size() == nb_vals ) ;
      for( size_t i=0 ; ok && i<nb_vals ; ++i )
      {
         ok &= ( vals(i) == double( root+10*i ) ) ;
      }
      notify_one_test_result( "broadcast(doubleVector)", ok ) ;
   }

   // ---------------
   // broadcast(doubleArray2D)
   {
      size_t const nb_vals = 10 ;
      doubleArray2D vals(0,0) ;
      if( rank == root )
      {
         vals.re_initialize( nb_vals, nb_vals ) ;
         for( size_t i=0 ; i<nb_vals ; ++i )
         {
            for( size_t j=0 ; j<nb_vals ; ++j )
            {
               vals(i,j) = double( rank+10*i-20*j ) ;
            }
         }
      }
      com->broadcast( vals, root ) ;
      bool ok = ( vals.index_bound(0) == nb_vals &&
                  vals.index_bound(1) == nb_vals ) ;
      for( size_t i=0 ; ok && i<nb_vals ; ++i )
      {
         for( size_t j=0 ; ok && j<nb_vals ; ++j )
         {
            ok &= ( vals(i,j) == double( root+10*i-20*j ) ) ;
         }
      }
      notify_one_test_result( "broadcast(doubleArray2D)", ok ) ;
   }

   // ---------------
   // broadcast(string)
   {
      std::string const string_val0 = "Hello world" ;
      std::string string_val = ( rank == root ? string_val0 : "other string" ) ;
      com->broadcast( string_val, root ) ;
      notify_one_test_result( "broadcast(string)",
                              string_val==string_val0 ) ;
   }

   // ---------------
   // broadcast(stringVector)
   {
      stringVector const stringVec_val0 = "Hello, ,world" ;
      stringVector stringVec_val =
         ( rank == root ? stringVec_val0 : "other, ,stringVector" ) ;
      com->broadcast( stringVec_val, root ) ;
      notify_one_test_result( "broadcast(stringVector)",
                              stringVec_val==stringVec_val0 ) ;
   }
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_gather( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_gather" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   size_t root = 0 ;

   // ---------------
   // gather(double)
   {
      double val = 1.0*rank ;
      doubleVector result( 0 ) ;
      if( rank == root ) result.re_initialize( nbrk ) ;
      com->gather( val, result, root ) ;
      bool ok = true ;
      if( rank == root )
      {
         for( size_t i=0 ; i<nbrk ; ++i )
         {
            ok &= ( result(i) == 1.0*i ) ;
         }
      }
      else
      {
         ok = ( result.size() == 0 ) ;
      }
      notify_one_test_result( "gather(double)", ok ) ;
   }

   // ---------------
   // gather(doubleVector)
   {
      size_t const nb_vals = 10 ;
      doubleVector dvals( nb_vals ) ;
      dvals.set( 1.0*rank ) ;
      doubleVector result( 0 ) ;
      if( rank == root ) result.re_initialize( nbrk * nb_vals ) ;
      com->gather( dvals, result, root ) ;
      bool ok = true ;
      if( rank == root )
      {
         for( size_t i=0 ; i<nbrk*nb_vals ; ++i )
         {
            ok &= ( result(i) == 1.0*(i/nb_vals) ) ;
         }
      }
      else
      {
         ok = ( result.size() == 0 ) ;
      }
      notify_one_test_result( "gather(doubleVector)", ok ) ;
   }

   // ---------------
   // gather(int)
   {
      int val = (int) rank ;
      intVector result( 0 ) ;
      if( rank == root ) result.re_initialize( nbrk ) ;
      com->gather( val, result, root ) ;
      bool ok = true ;
      if( rank == root )
      {
         for( size_t i=0 ; i<nbrk ; ++i )
         {
            ok &= ( result(i) == (int) i ) ;
         }
      }
      else
      {
         ok = ( result.size() == 0 ) ;
      }
      notify_one_test_result( "gather(int)", ok ) ;
   }

   // ---------------
   // gather(doubleVector)
   {
      size_t const nb_vals = 10 ;
      intVector ivals( nb_vals ) ;
      ivals.set( (int) rank ) ;
      intVector result( 0 ) ;
      if( rank == root ) result.re_initialize( nbrk * nb_vals ) ;
      com->gather( ivals, result, root ) ;
      bool ok = true ;
      if( rank == root )
      {
         for( size_t i=0 ; i<nbrk*nb_vals ; ++i )
         {
            ok &= ( result(i) == 1.0*(i/nb_vals) ) ;
         }
      }
      else
      {
         ok = ( result.size() == 0 ) ;
      }
      notify_one_test_result( "gather(doubleVector)", ok ) ;
   }
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_all_gather( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_all_gather" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   // ---------------
   // all_gather(int)
   int ival = rank ;

   intVector result( nbrk ) ;
   com->all_gather( ival, result ) ;

   bool ok = true ;
   for( size_t i=0 ; i<nbrk ; ++i )
   {
      ok &= ( result(i) == (int)i ) ;
   }
   notify_one_test_result( "all_gather(int)", ok ) ;

   // ---------------
   // all_gather(intVector)
   size_t const nb_vals = 10 ;
   intVector ivals( nb_vals ) ;
   for( size_t i=0 ; i<nb_vals ; ++i )
   {
      ivals(i) = nb_vals*rank + i ;
   }

   intVector iresult( nb_vals*nbrk ) ;
   com->all_gather( ivals, iresult ) ;

   ok = true ;
   for( size_t i=0 ; i<nb_vals*nbrk ; ++i )
   {
      ok = ok && ( iresult(i) == (int)i ) ;
   }

   notify_one_test_result( "all_gather(intVector)", ok ) ;

   // ---------------
   // all_gather_v(doubleVector)
   doubleVector dvals( nb_vals ) ;
   for( size_t i=0 ; i<nb_vals ; ++i )
   {
      dvals(i) = nb_vals*rank + i ;
   }
   intVector partition( nbrk ) ;
   intVector start( nbrk ) ;
   for( size_t i=0 ; i<nbrk ; i++ )
   {
      partition( i ) = (int)nb_vals ;
      start( i ) = (int)nb_vals*i ;
   }

   doubleVector dresult( nb_vals*nbrk ) ;
   com->all_gather_v( dvals.data(), nb_vals,
                      const_cast<double*>( dresult.data() ),
                      partition, start ) ;
   ok = true ;
   for( size_t i=0 ; i<nb_vals*nbrk ; i++ )
   {
      ok = ok && ( dresult(i) == 1.0*i ) ;
   }
   notify_one_test_result( "all_gather_v(doubleVector)", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_all_to_all( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_all_to_all" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   size_t const nb_vals = 10 ;
   intVector send_array( nb_vals*nbrk ) ;
   intVector recv_array( nb_vals*nbrk ) ;

   for( size_t p=0 ; p<nbrk ; ++p )
   {
      for( size_t i=0 ; i<nb_vals ; ++i )
      {
         send_array(p*nb_vals+i) = 50*rank+p-12*i ;
      }
   }

   com->all_to_all( send_array, recv_array ) ;

   bool ok = true ;
   for( size_t p=0 ; p<nbrk ; ++p )
   {
      for(  size_t i=0 ; i<nb_vals ; ++i )
      {
         ok = ok && recv_array(p*nb_vals+i) == (int)(50*p+rank-12*i) ;
      }
   }

   notify_one_test_result( "all_to_all(intVector)", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_boolean_and_or( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_boolean_and_or" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   bool btrue  = true ;
   bool bfalse = false ;
   bool bmixed = ( rank==0 ? false : true ) ;

   bool ok = ( !com->boolean_and( bfalse ) ) &&
             ( !com->boolean_and( bmixed ) ) &&
             (  com->boolean_and( btrue  ) ) ;
   notify_one_test_result( "boolean_and", ok ) ;

   ok = ( !com->boolean_or( bfalse ) ) &&
        ( (nbrk==1) || com->boolean_or( bmixed ) ) &&
        ( com->boolean_or( btrue ) ) ;
   notify_one_test_result( "boolean_or", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_sum( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_sum" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   double dval = 1.0*rank ;
   double dsum = 1.0*nbrk*(nbrk-1)/2.0 ;

   bool ok = ( com->sum( dval ) == dsum ) ;
   notify_one_test_result( "sum", ok ) ;

   size_t const nb_vals = 10 ;
   doubleVector dvals( nb_vals ) ;
   dvals.set( dval ) ;
   com->sum_vector( dvals ) ;
   ok = true ;
   for( size_t i=0 ; i<nb_vals ; ++i )
   {
      ok &= ( dvals(i) == dsum ) ;
   }
   notify_one_test_result( "sum_vector", ok ) ;

   size_t const nb_vals2 = 10 ;
   doubleArray2D dvals2( nb_vals, nb_vals2 ) ;
   dvals2.set( dval ) ;
   com->sum_array( dvals2 ) ;
   ok = true ;
   for( size_t i=0 ; i<nb_vals ; ++i )
   {
      for( size_t j=0 ; j<nb_vals2 ; ++j )
      {
         ok &= ( dvals2(i,j) == dsum ) ;
      }
   }
   notify_one_test_result( "sum_array", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_max_min( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_max_min" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   double dval = 1.0*rank ;

   bool ok = ( com->max( dval ) == 1.0*(nbrk-1) ) ;
   ok &= ( com->max( rank ) == nbrk - 1 ) ;
   notify_one_test_result( "max", ok ) ;

   ok = ( com->min( dval ) == 0.0 ) ;
   ok &= ( com->min( rank ) == 0 ) ;
   notify_one_test_result( "min", ok ) ;
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_same_value( PEL_Communicator const* com )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_same_value" ) ;

   size_t rank = com->rank() ;
   size_t nbrk = com->nb_ranks() ;

   {
      double val1 = 3.14 ;
      double val2 = ( rank==0 ? 12.0 : -4.0 ) ;

      bool ok = com->same_value_everywhere( val1 ) ;
      if( nbrk != 1 )
      {
         ok &= !com->same_value_everywhere( val2 ) ;
      }
      else
      {
         ok &= com->same_value_everywhere( val2 ) ;
      }
      notify_one_test_result( "same_value_eveywhere(double)", ok ) ;
   }

   {
      int val1 = 3 ;
      int val2 = ( rank==0 ? 12 : -4 ) ;

      bool ok = com->same_value_everywhere( val1 ) ;
      if( nbrk != 1 )
      {
         ok &= !com->same_value_everywhere( val2 ) ;
      }
      else
      {
         ok &= com->same_value_everywhere( val2 ) ;
      }
      notify_one_test_result( "same_value_eveywhere(int)", ok ) ;
   }
}

//-----------------------------------------------------------------------------
void
PEL_Communicator_TEST:: test_merge( PEL_Communicator const* com,
                                    PEL_ModuleExplorer const* exp )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Communicator_TEST:: test_merge" ) ;

   size_t rank = com->rank() ;

   if( exp->has_module( "merge" ) )
   {
      PEL_DoubleComparator const* dbl_cmp = PEL_DoubleComparatorExact::object() ;

      PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "merge" ) ;
      se->start_module_iterator() ;
      for( ; se->is_valid_module() ; se->go_next_module() )
      {
         bool ok = true ;
         PEL_ModuleExplorer const* sse = se->create_subexplorer( 0 ) ;
         doubleArray2D coord = sse->doubleArray2D_data( "coord" ) ;
         doubleArray2D merge = sse->doubleArray2D_data( "merge" ) ;

         size_t_vector idx( 0 ) ;
         doubleArray2D old_coord = coord ;
         com->merge( dbl_cmp, coord, idx ) ;

         ok &= ( coord.index_bound( 0 ) == old_coord.index_bound( 0 ) ) ;
         if( rank == 0 )
         {
            ok &= ( coord == merge ) ;
         }
         else
         {
            ok &= ( coord == old_coord ) ;
         }
         ok &= ( idx.size() == old_coord.index_bound( 1 ) ) ;
         for( size_t i=0 ; i<old_coord.index_bound( 1 ) ; ++i )
         {
            for( size_t d=0 ; d<old_coord.index_bound( 0 ) ; ++d )
            {
               ok &= ( old_coord( d, i ) == merge( d, idx( i ) ) ) ;
            }
         }
         notify_one_test_result( sse->name(), ok ) ;
         sse->destroy() ;
      }
      se->destroy() ; se = 0 ;
   }
}
