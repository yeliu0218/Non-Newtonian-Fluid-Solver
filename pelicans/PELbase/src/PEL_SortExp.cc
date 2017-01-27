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

#include <PEL_SortExp.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_BalancedBinaryTree.hh>
#include <PEL_Double.hh>
#include <PEL_Error.hh>
#include <PEL_Vector.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL_String.hh>

#include <boolVector.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <iostream>

PEL_SortExp const* PEL_SortExp::PROTOTYPE = new PEL_SortExp() ;

//----------------------------------------------------------------------
PEL_SortExp:: PEL_SortExp( void ) 
//----------------------------------------------------------------------
   : PEL_Expression( "sort" )
   , VECTOR( 0 )
   , GROWING( false )
   , RESULT( 0 )
{
   PEL_LABEL( "PEL_SortExp:: PEL_SortExp" ) ;
}

//----------------------------------------------------------------------
PEL_SortExp*
PEL_SortExp:: create_replica( PEL_Object* a_owner,
                                 PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SortExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_SortExp* result = new PEL_SortExp( a_owner, argument_list ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_SortExp:: PEL_SortExp( PEL_Object* a_owner,
                                 PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, "sort", argument_list )
   , VECTOR( arg(0) )
   , GROWING( arg(1)->to_string() == "<" )
   , RESULT( 0 )
{
   PEL_LABEL( "PEL_SortExp:: PEL_SortExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_SortExp:: ~PEL_SortExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SortExp:: ~PEL_SortExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
std::string const&
PEL_SortExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result = "sort(DV,\">\"|\"<\")" ;
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_SortExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SortExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = some_arguments->count()==2 &&
      extract_arg( some_arguments, 0 )->data_type()==DoubleVector &&
      extract_arg( some_arguments, 1 )->data_type()==String ;
   if( result )
   {
      std::string const& opt = extract_arg( some_arguments, 1 )->to_string() ;
      result = opt=="<" || opt==">" ;
   }
   return result ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_SortExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SortExp:: data_type" ) ;
   return PEL_Data::DoubleVector ;
}

//----------------------------------------------------------------------
doubleVector const&
PEL_SortExp:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SortExp:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ;
   doubleVector const& initial = VECTOR->to_double_vector( ct ) ;
   PEL_BalancedBinaryTree* tree =  PEL_BalancedBinaryTree::create( 0 ) ;
   for( size_t j=0 ; j<initial.size() ; j++ )
   {
      PEL_Double* d = PEL_Double::create( tree, initial(j) ) ;
      tree->extend( d) ;
   }
   doubleVector& result = RESULT ;
   size_t n = tree->count() ;
   result.re_initialize( n ) ;
   PEL_Iterator* it = tree->create_iterator( tree ) ;
   size_t i = 0 ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      size_t idx = ( GROWING ? i : n-1-i ) ;
      result( idx ) = static_cast<PEL_Data*>( it->item() )->to_double() ;
      i++ ;
   }
   tree->destroy() ;
   return result ;
}
