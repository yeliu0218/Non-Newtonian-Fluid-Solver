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

#include <PEL_DataWithContext.hh>

#include <PEL_assertions.hh>
#include <PEL_Context.hh>
#include <PEL_ContextPair.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContextExp.hh>
#include <PEL_List.hh>
#include <PEL_Sequence.hh>
#include <PEL_Variable.hh>

#include <iostream>

//----------------------------------------------------------------------------
PEL_DataWithContext*
PEL_DataWithContext:: create( PEL_Object* a_owner,
                              PEL_Data const* data,
                              PEL_Context const* ct )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: create" ) ;
   PEL_CHECK_PRE( data != 0 ) ;
   PEL_CHECK_PRE( ct != 0 ) ;

   PEL_DataWithContext* result = new PEL_DataWithContext( a_owner, data, ct ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->data_type() == data->data_type() ) ;
   PEL_CHECK_POST( result->is_raw_data() == data->is_raw_data() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_DataWithContext:: PEL_DataWithContext( PEL_Object* a_owner,
                                           PEL_Data const* data,
                                           PEL_Context const* ct )
//----------------------------------------------------------------------------
   : PEL_Data( a_owner )
   , DATA( data )
   , CTX( ct->create_clone( this ) )
   , TMP_CTX( PEL_ContextPair::create( this,0,0 ) )
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
PEL_DataWithContext*
PEL_DataWithContext:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: create_clone" ) ;

   PEL_DataWithContext* result = new PEL_DataWithContext( a_owner, DATA, CTX ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->data_type() == data_type() ) ;
   PEL_CHECK_POST( result->is_constant() == is_constant() ) ;
   PEL_CHECK_POST( result->is_raw_data() == is_raw_data() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_DataWithContext:: ~PEL_DataWithContext( void ) 
//----------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
PEL_DataWithContext:: declare( PEL_List* lst ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: declare" ) ;
   PEL_CHECK_PRE( declare_PRE( lst ) ) ;

   PEL_List* ctx = PEL_List::create( 0 ) ;
   DATA->declare( ctx ) ;
   for( size_t i=0 ; i<ctx->count() ; i++ )
   {
      PEL_Variable* var = static_cast<PEL_Variable*>( ctx->at(i) ) ;
      if( !CTX->has_variable( var ) )
      {
         lst->extend( var ) ;
      }
      else
      {
         PEL_DataWithContext* data =
            PEL_DataWithContext::create( 0,
                                         CTX->value( var ),
                                         CTX ) ;
         data->declare( lst ) ;
         data->destroy() ; data = 0 ;
      }
   }
   ctx->destroy() ; ctx=0 ;
   
   PEL_CHECK_POST( declare_POST( lst ) ) ;
}

//----------------------------------------------------------------------------
PEL_Data::Type
PEL_DataWithContext:: data_type( void ) const
//----------------------------------------------------------------------------
{
   return( DATA->data_type() ) ;
}

//----------------------------------------------------------------------
PEL_Context const*
PEL_DataWithContext:: context( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: context" ) ;
   PEL_Context const* result = CTX ;
   if( ct!=0 && ct!=CTX )
   {
      TMP_CTX->re_initialize( CTX, ct ) ;
      result = TMP_CTX ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_DataWithContext:: value_can_be_evaluated( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: value_can_be_evaluated" ) ;
   return( DATA->value_can_be_evaluated( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
stringVector const&
PEL_DataWithContext:: undefined_variables( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: undefined_variables" ) ;
   return( DATA->undefined_variables( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
bool
PEL_DataWithContext:: to_bool( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: to_bool" ) ;
   PEL_CHECK_PRE( to_bool_PRE( ct ) ) ;
   return( DATA->to_bool( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
double
PEL_DataWithContext:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE( ct ) ) ;
   return( DATA->to_double( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
int
PEL_DataWithContext:: to_int( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: to_int" ) ;
   PEL_CHECK_PRE( to_int_PRE( ct ) ) ;
   return( DATA->to_int( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
std::string const&
PEL_DataWithContext:: to_string( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: to_string" ) ;
   PEL_CHECK_PRE( to_string_PRE( ct ) ) ;
   return( DATA->to_string( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
doubleVector const& 
PEL_DataWithContext:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ;
   return( DATA->to_double_vector( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
intVector const&
PEL_DataWithContext:: to_int_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: to_int_vector" ) ;
   PEL_CHECK_PRE( to_int_vector_PRE( ct ) ) ;
   return( DATA->to_int_vector( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
stringVector const& 
PEL_DataWithContext:: to_string_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: to_string_vector" ) ;
   PEL_CHECK_PRE( to_string_vector_PRE( ct ) ) ;
   return( DATA->to_string_vector( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
boolVector const&
PEL_DataWithContext:: to_bool_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: to_bool_vector" ) ;
   PEL_CHECK_PRE( to_bool_vector_PRE( ct ) ) ;
   return( DATA->to_bool_vector( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
doubleArray2D const&
PEL_DataWithContext:: to_double_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: to_double_array2D" ) ;
   PEL_CHECK_PRE( to_double_array2D_PRE( ct ) ) ;
   return( DATA->to_double_array2D( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
boolArray2D const&
PEL_DataWithContext:: to_bool_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: to_bool_array2D" ) ;
   PEL_CHECK_PRE( to_bool_array2D_PRE( ct ) ) ;
   return( DATA->to_bool_array2D( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
stringArray2D const&
PEL_DataWithContext:: to_string_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: to_string_array2D" ) ;
   PEL_CHECK_PRE( to_string_array2D_PRE( ct ) ) ;
   return( DATA->to_string_array2D( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
intArray2D const&
PEL_DataWithContext:: to_int_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: to_int_array2D" ) ;
   PEL_CHECK_PRE( to_int_array2D_PRE( ct ) ) ;
   return( DATA->to_int_array2D( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
doubleArray3D const&
PEL_DataWithContext:: to_double_array3D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: to_double_array3D" ) ;
   PEL_CHECK_PRE( to_double_array3D_PRE( ct ) ) ;
   return( DATA->to_double_array3D( context( ct ) ) ) ;
}

//----------------------------------------------------------------------
intArray3D const&
PEL_DataWithContext:: to_int_array3D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: to_int_array3D" ) ;
   PEL_CHECK_PRE( to_int_array3D_PRE( ct ) ) ;
   return( DATA->to_int_array3D(context( ct ) ) ) ;
}

//----------------------------------------------------------------------
bool
PEL_DataWithContext:: is_raw_data( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: is_raw_data" ) ;
   bool result = DATA->is_raw_data() ;
   PEL_CHECK_POST( is_raw_data_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_DataWithContext:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DataWithContext:: print" ) ;

   if( CTX->nb_variables() != 0 )
   {
      PEL_DataWithContextExp const* exp =
              PEL_DataWithContextExp::create( 0, DATA, CTX ) ;
      exp->print( os, indent_width ) ;
      exp->destroy() ; exp = 0 ;
   }
   else
   {
      DATA->print( os, indent_width ) ;
   }
}
