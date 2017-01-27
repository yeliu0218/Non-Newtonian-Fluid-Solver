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

#include <PEL_VariableExp.hh>

#include <PEL_assertions.hh>
#include <PEL_Context.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>
#include <PEL_Variable.hh>

#include <boolVector.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <stringVector.hh>
#include <boolArray2D.hh>
#include <stringArray2D.hh>
#include <doubleArray2D.hh>
#include <intArray2D.hh>

PEL_VariableExp const*  
PEL_VariableExp::PROTOTYPE_var_def = 
                             new PEL_VariableExp( var_def, "is_defined" ) ;

PEL_VariableExp const*  
PEL_VariableExp::PROTOTYPE_var_value = 
                             new PEL_VariableExp( var_value, "value" ) ;

struct PEL_VariableExp_ERROR
{
   static void n0( std::string const& f_name, std::string const& op_name ) ;
   static void n1( std::string const& a_name ) ;
   static void n2( std::string const& a_name,
                   stringVector const& undef_var ) ;
} ;

//----------------------------------------------------------------------
PEL_VariableExp:: PEL_VariableExp( VarExp exp_id,
                                   std::string const& a_name  ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , OP( exp_id )
{
   PEL_LABEL( "PEL_VariableExp:: PEL_VariableExp" ) ;
}

//----------------------------------------------------------------------
PEL_VariableExp*
PEL_VariableExp:: create_replica( PEL_Object* a_owner,
                                  PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VariableExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_VariableExp* result = new PEL_VariableExp( a_owner,
                                                  OP,
                                                  name(),
                                                  argument_list ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_VariableExp:: PEL_VariableExp( PEL_Object* a_owner,
                                   VarExp exp_id,
                                   std::string const& a_name,
                                   PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , OP( exp_id )
{
   PEL_LABEL( "PEL_VariableExp:: PEL_VariableExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_VariableExp:: ~PEL_VariableExp( void ) 
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
std::string const& 
PEL_VariableExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   switch( OP )
   {
      case var_def :
         result = name() + "(SS)" ;
	 break ;
      case var_value :
         result = name() + "(SS,<default value>)" ;
	 break ;
      default :
         PEL_VariableExp_ERROR::n0( "usage", name() ) ;
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_VariableExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VariableExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   bool result = false ;
   switch( OP )
   {
      case var_def :
         result = some_arguments->count()==1 ;
	 if( result )
	 {
	    Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
	    result = ( t0==String ) ;
	 }
	 break ;
      case var_value :
         result = some_arguments->count()==2 ;
	 if( result )
	 {
	    Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
	    result = ( t0==String ) ;
	 }
	 break ;
      default :
         PEL_VariableExp_ERROR::n0( "valid_arguments", name() ) ;
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_VariableExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VariableExp:: data_type" ) ;
   
   PEL_Data::Type result = PEL_Data::Undefined ;
   switch( OP )
   {
      case var_def :
         result = PEL_Data::Bool ;
         break ;
      case var_value :
         result = arg(1)->data_type() ;
	 break ;
      default :
         PEL_VariableExp_ERROR::n0( "data_type", name() ) ;
         break ;   
   }   
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_VariableExp:: is_constant( void ) const
//----------------------------------------------------------------------
{
   return( false ) ;
}

//----------------------------------------------------------------------
bool
PEL_VariableExp:: to_bool( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VariableExp:: to_bool" ) ;
   PEL_CHECK_PRE( to_bool_PRE( ct ) ) ;

   bool result = false ;
   switch( OP )
   {
      case var_def :
         if( ct != 0 )
         {
            PEL_Variable const* var =
                       PEL_Variable::object( arg(0)->to_string( ct ) ) ;
            result = ct->has_variable( var ) ;
         }
         break ;
      case var_value :
      {
         PEL_Variable const* var =
                       PEL_Variable::object( arg(0)->to_string( ct ) ) ; 
         if( var->data_type() != arg(1)->data_type() )
         {
            PEL_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }   
            result = var->to_bool( ct ) ;
         }
         else
         {
            result = arg(1)->to_bool( ct ) ;
         }
      }
      break ;
      default :
         PEL_VariableExp_ERROR::n0( "to_bool", name() ) ;
         break ;   
   }
   return( result ) ;
}
   
//----------------------------------------------------------------------
double
PEL_VariableExp:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VariableExp:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE( ct ) ) ;

   double result = 0. ;
   switch( OP )
   {
      case var_value :
      {
         PEL_Variable const* var =
                       PEL_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            PEL_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }   
            result = var->to_double( ct ) ;
         }
         else
         {
            result = arg(1)->to_double( ct ) ;
         }
      }
      break ;
      default :
         PEL_VariableExp_ERROR::n0( "to_double", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
int
PEL_VariableExp:: to_int( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VariableExp:: to_int" ) ;
   PEL_CHECK_PRE( to_int_PRE( ct ) ) ;

   int result = 0 ;
   switch( OP )
   {
      case var_value :
      {
         PEL_Variable const* var =
                       PEL_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            PEL_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_int( ct ) ;
         }
         else
         {
            result = arg(1)->to_int( ct ) ;
         }
      }
      break ;
      default :
         PEL_VariableExp_ERROR::n0( "to_int", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
PEL_VariableExp:: to_string( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VariableExp:: to_string" ) ;
   PEL_CHECK_PRE( to_string_PRE( ct ) ) ;

   static std::string result = "unexpected" ;
   switch( OP )
   {
      case var_value :
      {
         PEL_Variable const* var =
                       PEL_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            PEL_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_string( ct ) ;
         }
         else
         {
            result = arg(1)->to_string( ct ) ;
         }
      }
      break ;
      default :
         PEL_VariableExp_ERROR::n0( "to_string", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
boolVector const&
PEL_VariableExp:: to_bool_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VariableExp:: to_bool_vector" ) ;
   PEL_CHECK_PRE( to_bool_vector_PRE( ct ) ) ;

   static boolVector result(0) ;
   switch( OP )
   {
      case var_value :
      {
         PEL_Variable const* var =
                       PEL_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            PEL_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_bool_vector( ct ) ;
         }
         else
         {
            result = arg(1)->to_bool_vector( ct ) ;
         }
      }
      break ;
      default :
         PEL_VariableExp_ERROR::n0( "to_bool_vector", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
doubleVector const&
PEL_VariableExp:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VariableExp:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ;

   static doubleVector result(0) ;
   switch( OP )
   {
      case var_value :
      {
         PEL_Variable const* var =
                       PEL_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            PEL_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_double_vector( ct ) ;
         }
         else
         {
            result = arg(1)->to_double_vector( ct ) ;
         }
      }
      break ;
      default :
         PEL_VariableExp_ERROR::n0( "to_double_vector", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
intVector const&
PEL_VariableExp:: to_int_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VariableExp:: to_int_vector" ) ;
   PEL_CHECK_PRE( to_int_vector_PRE( ct ) ) ;

   static intVector result(0) ;
   switch( OP )
   {
      case var_value :
      {
         PEL_Variable const* var =
                       PEL_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            PEL_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_int_vector( ct ) ;
         }
         else
         {
            result = arg(1)->to_int_vector( ct ) ;
         }
      }
      break ;
      default :
         PEL_VariableExp_ERROR::n0( "to_int_vector", name() ) ;
         break ;   
   }
   return( result ) ;
}
   
//----------------------------------------------------------------------
stringVector const&
PEL_VariableExp:: to_string_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VariableExp:: to_string_vector" ) ;
   PEL_CHECK_PRE( to_string_vector_PRE( ct ) ) ;

   static stringVector result(0) ;
   switch( OP )
   {
      case var_value :
      {
         PEL_Variable const* var =
                       PEL_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            PEL_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_string_vector( ct ) ;
         }
         else
         {
            result = arg(1)->to_string_vector( ct ) ;
         }
      }
      break ;
      default :
         PEL_VariableExp_ERROR::n0( "to_string_vector", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
doubleArray2D const&
PEL_VariableExp:: to_double_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VariableExp:: to_double_array2D" ) ;
   PEL_CHECK_PRE( to_double_array2D_PRE( ct ) ) ;

   static doubleArray2D result(0,0) ;
   switch( OP )
   {
      case var_value :
      {
         PEL_Variable const* var =
                       PEL_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            PEL_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_double_array2D( ct ) ;
         }
         else
         {
            result = arg(1)->to_double_array2D( ct ) ;
         }
      }
      break ;
      default :
         PEL_VariableExp_ERROR::n0( "to_double_array2D", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
stringArray2D const&
PEL_VariableExp:: to_string_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VariableExp:: to_string_array2D" ) ;
   PEL_CHECK_PRE( to_string_array2D_PRE( ct ) ) ;

   static stringArray2D result(0,0) ;
   switch( OP )
   {
      case var_value :
      {
         PEL_Variable const* var =
                       PEL_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            PEL_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_string_array2D( ct ) ;
         }
         else
         {
            result = arg(1)->to_string_array2D( ct ) ;
         }
      }
      break ;
      default :
         PEL_VariableExp_ERROR::n0( "to_string_array2D", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
boolArray2D const&
PEL_VariableExp:: to_bool_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VariableExp:: to_bool_array2D" ) ;
   PEL_CHECK_PRE( to_bool_array2D_PRE( ct ) ) ;

   static boolArray2D result(0,0) ;
   switch( OP )
   {
      case var_value :
      {
         PEL_Variable const* var =
                       PEL_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            PEL_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_bool_array2D( ct ) ;
         }
         else
         {
            result = arg(1)->to_bool_array2D( ct ) ;
         }
      }
      break ;
      default :
         PEL_VariableExp_ERROR::n0( "to_bool_array2D", name() ) ;
         break ;   
   }
   return( result ) ;
}

//----------------------------------------------------------------------
intArray2D const&
PEL_VariableExp:: to_int_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_VariableExp:: to_int_array2D" ) ;
   PEL_CHECK_PRE( to_int_array2D_PRE( ct ) ) ;

   static intArray2D result(0,0) ;
   switch( OP )
   {
      case var_value :
      {
         PEL_Variable const* var =
                       PEL_Variable::object( arg(0)->to_string( ct ) ) ;
         if( var->data_type() != arg(1)->data_type() )
         {
            PEL_VariableExp_ERROR:: n1( arg(0)->to_string( ct ) ) ;
         }
         if( !arg(1)->value_can_be_evaluated( ct ) )
         {
            PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                        arg(1)->undefined_variables( ct ) ) ;
         }
         if( ct != 0 && ct->has_variable( var ) )
         {
            if( !var->value_can_be_evaluated( ct ) )
            {
               PEL_VariableExp_ERROR:: n2( arg(0)->to_string( ct ),
                                           var->undefined_variables( ct ) ) ;
            }
            result = var->to_int_array2D( ct ) ;
         }
         else
         {
            result = arg(1)->to_int_array2D( ct ) ;
         }
      }
      break ;
      default :
         PEL_VariableExp_ERROR::n0( "to_int_array2D", name() ) ;
         break ;   
   }
   return( result ) ;
}

//internal--------------------------------------------------------------
void 
PEL_VariableExp_ERROR:: n0( std::string const& f_name,
                            std::string const& op_name )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** PEL_VariableExp::" + f_name +"\n" ;
   mesg += "    operation " + op_name + " not implemented." ;
   PEL_Error::object()->raise_internal( mesg ) ;
}

//internal--------------------------------------------------------------
void 
PEL_VariableExp_ERROR:: n1( std::string const& a_name )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** PEL_VariableExp:: value error\n" ;
   mesg += "    variable name \""+a_name+"\" and default value should have the same type." ;
   PEL_Error::object()->raise_plain( mesg ) ;
}

//internal--------------------------------------------------------------
void 
PEL_VariableExp_ERROR:: n2( std::string const& a_name,
                            stringVector const& undef_var )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** PEL_VariableExp:: value error\n" ;
   mesg += "    variable name \""+a_name+"\" cannot be evaluated.\n" ;
   if( undef_var.size() > 0 )
   {
      mesg += "    undefined variable(s):\n" ;
      for( size_t i=0 ; i<undef_var.size() ; ++i )
      {
         mesg += "       - \"" + undef_var(i) + "\"\n" ;
      }
   }
   PEL_Error::object()->raise_plain( mesg ) ;
}
