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

#include <PEL_ArrayExp.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>

#include <iostream>

PEL_ArrayExp const* PEL_ArrayExp::PROTOTYPE = new PEL_ArrayExp() ;

//----------------------------------------------------------------------
PEL_ArrayExp:: PEL_ArrayExp( void ) 
//----------------------------------------------------------------------
   : PEL_Expression( "array" )
   , RESULT_D2D( 0, 0 )
   , RESULT_I2D( 0, 0 )
   , RESULT_B2D( 0, 0 )
   , RESULT_S2D( 0, 0 )
   , RESULT_D3D( 0, 0, 0 )
   , RESULT_I3D( 0, 0, 0 )
{
   PEL_LABEL( "PEL_ArrayExp:: PEL_ArrayExp" ) ;
}

//----------------------------------------------------------------------
PEL_ArrayExp:: PEL_ArrayExp( PEL_Object* a_owner,
                             PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, "array", argument_list )
   , RESULT_D2D( 0, 0 )
   , RESULT_I2D( 0, 0 )
   , RESULT_B2D( 0, 0 )
   , RESULT_S2D( 0, 0 )
   , RESULT_D3D( 0, 0, 0 )
   , RESULT_I3D( 0, 0, 0 )
{
   PEL_LABEL( "PEL_ArrayExp:: PEL_ArrayExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_ArrayExp:: ~PEL_ArrayExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArrayExp:: ~PEL_ArrayExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_ArrayExp*
PEL_ArrayExp:: create_replica( PEL_Object* a_owner,
                               PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArrayExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_ArrayExp* result = new PEL_ArrayExp( a_owner, argument_list ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
PEL_ArrayExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result = "array(<list of IV|DV|SV|BV>)" ;
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_ArrayExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArrayExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = some_arguments->count()>0 ;
   bool prem = true ;
   Type k = Undefined ;
   for( size_t idx=0 ; result && idx<some_arguments->count() ; ++idx )
   {
      if( prem )
      {
         k = extract_arg( some_arguments, idx )->data_type() ;
         if( k!=IntVector &&  k!=DoubleVector &&
             k!=IntArray2D &&  k!=DoubleArray2D &&
             k!=BoolVector &&  k!=StringVector )
         {
            result = false ;
         }
         prem = false ;
      }
      else
      {
         
         result &= ( k==extract_arg( some_arguments, idx )->data_type() ) ;
      }
   }
   
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_ArrayExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArrayExp:: data_type" ) ;
   PEL_Data::Type my_type = Undefined ;
   if( arg(0)->data_type()==IntVector )
   {
      my_type = IntArray2D ;
   }
   else if( arg(0)->data_type()==DoubleVector ) 
   {
      my_type = DoubleArray2D ;
   }
   else if( arg(0)->data_type()==IntArray2D ) 
   {
      my_type = IntArray3D ;
   }
   else if( arg(0)->data_type()==DoubleArray2D ) 
   {
      my_type = DoubleArray3D ;
   }
   else if( arg(0)->data_type()==BoolVector ) 
   {
      my_type = BoolArray2D ;
   }
   else if( arg(0)->data_type()==StringVector ) 
   {
      my_type = StringArray2D ;
   }
   else
   {
      PEL_Error::object()->raise_internal( "Bad type " ) ;
   }
   return my_type ;
}

//----------------------------------------------------------------------
doubleArray2D const&
PEL_ArrayExp:: to_double_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArrayExp:: to_double_array2D" ) ;
   PEL_CHECK_PRE( to_double_array2D_PRE( ct ) ) ;
   
   size_t nb_col = 0 ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      nb_col = PEL::max( nb_col, arg(idx)->to_double_vector(ct).size() ) ;
   }
   RESULT_D2D.re_initialize( nb_arguments(), nb_col ) ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      doubleVector const& vec = arg(idx)->to_double_vector( ct ) ;
      for( size_t j=0 ; j<vec.size() ; ++j )
      {
         RESULT_D2D( idx, j ) = vec( j ) ;
      }
   }
   return( RESULT_D2D ) ;
}

//----------------------------------------------------------------------
intArray2D const&
PEL_ArrayExp:: to_int_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArrayExp:: to_int_array2D" ) ;
   PEL_CHECK_PRE( to_int_array2D_PRE( ct ) ) ;
   
   size_t nb_col = 0 ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      nb_col = PEL::max( nb_col, arg(idx)->to_int_vector(ct).size() ) ;
   }
   RESULT_I2D.re_initialize( nb_arguments(), nb_col ) ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      intVector const& vec = arg(idx)->to_int_vector( ct ) ;
      for( size_t j=0 ; j<vec.size() ; ++j )
      {
         RESULT_I2D( idx, j ) = vec( j ) ;
      }
   }
   return( RESULT_I2D ) ;
}


//----------------------------------------------------------------------
boolArray2D const&
PEL_ArrayExp:: to_bool_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArrayExp:: to_bool_array2D" ) ;
   PEL_CHECK_PRE( to_bool_array2D_PRE( ct ) ) ;
   
   size_t nb_col = 0 ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      nb_col = PEL::max( nb_col, arg(idx)->to_bool_vector(ct).size() ) ;
   }
   RESULT_B2D.re_initialize( nb_arguments(), nb_col ) ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      boolVector const& vec = arg(idx)->to_bool_vector( ct ) ;
      for( size_t j=0 ; j<vec.size() ; ++j )
      {
         RESULT_B2D( idx, j ) = vec( j ) ;
      }
   }
   return( RESULT_B2D ) ;
}


//----------------------------------------------------------------------
stringArray2D const&
PEL_ArrayExp:: to_string_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArrayExp:: to_string_array2D" ) ;
   PEL_CHECK_PRE( to_string_array2D_PRE( ct ) ) ;
   
   size_t nb_col = 0 ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      nb_col = PEL::max( nb_col, arg(idx)->to_string_vector(ct).size() ) ;
   }
   RESULT_S2D.re_initialize( nb_arguments(), nb_col ) ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      stringVector const& vec = arg(idx)->to_string_vector( ct ) ;
      for( size_t j=0 ; j<vec.size() ; ++j )
      {
         RESULT_S2D( idx, j ) = vec( j ) ;
      }
   }
   return( RESULT_S2D ) ;
}

//----------------------------------------------------------------------
doubleArray3D const&
PEL_ArrayExp:: to_double_array3D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArrayExp:: to_double_array3D" ) ;
   PEL_CHECK_PRE( to_double_array3D_PRE( ct ) ) ;
   
   size_t dim1 = 0 ;
   size_t dim2 = 0 ;

   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      doubleArray2D const& a = arg(idx)->to_double_array2D(ct) ;
      dim1 = PEL::max( dim1, a.index_bound(0) ) ;
      dim2 = PEL::max( dim2, a.index_bound(1) ) ;
   }
   RESULT_D3D.re_initialize( nb_arguments(), dim1, dim2 ) ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      doubleArray2D const& a = arg(idx)->to_double_array2D(ct) ;
      for( size_t j=0 ; j<a.index_bound(0) ; ++j )
      {
         for( size_t k=0 ; k<a.index_bound(1) ; ++k )
         {
            RESULT_D3D( idx, j, k ) = a( j, k ) ;
         }
      }
   }
   return( RESULT_D3D ) ;
}

//----------------------------------------------------------------------
intArray3D const&
PEL_ArrayExp:: to_int_array3D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArrayExp:: to_int_array3D" ) ;
   PEL_CHECK_PRE( to_int_array3D_PRE( ct ) ) ;
   
   size_t dim1 = 0 ;
   size_t dim2 = 0 ;
   
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      intArray2D const& a = arg(idx)->to_int_array2D(ct) ;
      dim1 = PEL::max( dim1, a.index_bound(0) ) ;
      dim2 = PEL::max( dim2, a.index_bound(1) ) ;
   }
   RESULT_I3D.re_initialize( nb_arguments(), dim1, dim2 ) ;
   for( size_t idx=0 ; idx<nb_arguments() ; ++idx )
   {
      intArray2D const& a = arg(idx)->to_int_array2D(ct) ;
      for( size_t j=0 ; j<a.index_bound(0) ; j++ )
      {
         for( size_t k=0 ; k<a.index_bound(1) ; k++ )
         {
            RESULT_I3D( idx, j, k ) = a( j, k ) ;
         }
      }
   }
   return( RESULT_I3D ) ;
}
