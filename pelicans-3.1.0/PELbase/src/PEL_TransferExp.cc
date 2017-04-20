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

#include <PEL_TransferExp.hh>

#include <PEL_assertions.hh>
#include <PEL_Context.hh>
#include <PEL_Sequence.hh>

//----------------------------------------------------------------------
PEL_TransferExp:: PEL_TransferExp( std::string const& a_name ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
{
}

//----------------------------------------------------------------------
PEL_TransferExp:: PEL_TransferExp( PEL_Object* a_owner,
                                   std::string const& a_name,
                                   PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_TransferExp:: ~PEL_TransferExp( void ) 
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
bool
PEL_TransferExp:: to_bool( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TransferExp:: to_bool" ) ;
   PEL_CHECK_PRE( to_bool_PRE( ct ) ) ;
   return( data( ct )->to_bool( ct ) ) ;
}

//----------------------------------------------------------------------
double
PEL_TransferExp:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TransferExp:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE( ct ) ) ;
   return( data( ct )->to_double( ct ) ) ;
}

//----------------------------------------------------------------------
int
PEL_TransferExp:: to_int( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TransferExp:: to_int" ) ;
   PEL_CHECK_PRE( to_int_PRE( ct ) ) ;
   return( data( ct )->to_int( ct ) ) ;
}

//----------------------------------------------------------------------
std::string const&
PEL_TransferExp:: to_string( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TransferExp:: to_string" ) ;
   PEL_CHECK_PRE( to_string_PRE( ct ) ) ;
   return( data( ct )->to_string( ct ) ) ;
}

//----------------------------------------------------------------------
doubleVector const& 
PEL_TransferExp:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TransferExp:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ;
   return( data( ct )->to_double_vector( ct ) ) ;
}

//----------------------------------------------------------------------
intVector const&
PEL_TransferExp:: to_int_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TransferExp:: to_int_vector" ) ;
   PEL_CHECK_PRE( to_int_vector_PRE( ct ) ) ;
   return( data( ct )->to_int_vector( ct ) ) ;
}

//----------------------------------------------------------------------
stringVector const& 
PEL_TransferExp:: to_string_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TransferExp:: to_string_vector" ) ;
   PEL_CHECK_PRE( to_string_vector_PRE( ct ) ) ;
   return( data( ct )->to_string_vector( ct ) ) ;
}

//----------------------------------------------------------------------
boolVector const&
PEL_TransferExp:: to_bool_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TransferExp:: to_bool_vector" ) ;
   PEL_CHECK_PRE( to_bool_vector_PRE( ct ) ) ;
   return( data( ct )->to_bool_vector( ct ) ) ;
}

//----------------------------------------------------------------------
doubleArray2D const&
PEL_TransferExp:: to_double_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TransferExp:: to_double_array2D" ) ;
   PEL_CHECK_PRE( to_double_array2D_PRE( ct ) ) ;
   return( data( ct )->to_double_array2D( ct ) ) ;
}

//----------------------------------------------------------------------
boolArray2D const&
PEL_TransferExp:: to_bool_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TransferExp:: to_bool_array2D" ) ;
   PEL_CHECK_PRE( to_bool_array2D_PRE( ct ) ) ;
   return( data( ct )->to_bool_array2D( ct ) ) ;
}

//----------------------------------------------------------------------
stringArray2D const&
PEL_TransferExp:: to_string_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TransferExp:: to_string_array2D" ) ;
   PEL_CHECK_PRE( to_string_array2D_PRE( ct ) ) ;
   return( data( ct )->to_string_array2D( ct ) ) ;
}

//----------------------------------------------------------------------
intArray2D const&
PEL_TransferExp:: to_int_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TransferExp:: to_int_array2D" ) ;
   PEL_CHECK_PRE( to_int_array2D_PRE( ct ) ) ;
   return( data( ct )->to_int_array2D( ct ) ) ;
}

//----------------------------------------------------------------------
doubleArray3D const&
PEL_TransferExp:: to_double_array3D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TransferExp:: to_double_array3D" ) ;
   PEL_CHECK_PRE( to_double_array3D_PRE( ct ) ) ;
   return( data( ct )->to_double_array3D( ct ) ) ;
}

//----------------------------------------------------------------------
intArray3D const&
PEL_TransferExp:: to_int_array3D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TransferExp:: to_int_array3D" ) ;
   PEL_CHECK_PRE( to_int_array3D_PRE( ct ) ) ;
   return( data( ct )->to_int_array3D( ct ) ) ;
}

//----------------------------------------------------------------------
bool
PEL_TransferExp:: data_PRE( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( !is_a_prototype() ) ;
   PEL_ASSERT( value_can_be_evaluated( ct ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PEL_TransferExp:: data_POST( PEL_Data const* result,
                             PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->value_can_be_evaluated( ct ) ) ;
   return( true ) ;
}










