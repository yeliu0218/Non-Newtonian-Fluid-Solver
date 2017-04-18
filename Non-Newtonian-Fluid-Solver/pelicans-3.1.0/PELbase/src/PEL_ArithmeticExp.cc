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

#include <PEL_ArithmeticExp.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Double.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_Sequence.hh>

#include <iostream>

PEL_ArithmeticExp const* 
PEL_ArithmeticExp::PROTOTYPE_M = new PEL_ArithmeticExp( "+", M ) ;

PEL_ArithmeticExp const* 
PEL_ArithmeticExp::PROTOTYPE_L = new PEL_ArithmeticExp( "-", L ) ;

PEL_ArithmeticExp const* 
PEL_ArithmeticExp::PROTOTYPE_T = new PEL_ArithmeticExp( "*", T ) ;

PEL_ArithmeticExp const* 
PEL_ArithmeticExp::PROTOTYPE_D = new PEL_ArithmeticExp( "/", D ) ;

PEL_ArithmeticExp const* 
PEL_ArithmeticExp::PROTOTYPE_MOD = new PEL_ArithmeticExp( "modulo", MOD ) ;

//----------------------------------------------------------------------
PEL_ArithmeticExp:: PEL_ArithmeticExp( std::string const& a_name,
                                       AlgebraicOperator a_op ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , OP( a_op )
   , ARG0( 0 )
   , ARG1( 0 )
{
   PEL_LABEL( "PEL_ArithmeticExp:: PEL_ArithmeticExp" ) ;
}

//----------------------------------------------------------------------
PEL_ArithmeticExp*
PEL_ArithmeticExp:: create_replica( PEL_Object* a_owner,
                                    PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArithmeticExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_ArithmeticExp* result = new PEL_ArithmeticExp( a_owner,
                                                      name(),
                                                      argument_list,
                                                      OP ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ArithmeticExp:: PEL_ArithmeticExp( PEL_Object* a_owner,
                                       std::string const& a_name,
                                       PEL_Sequence const* argument_list,
                                       AlgebraicOperator a_op ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , OP( a_op )
   , ARG0( arg(0) )
   , ARG1( arg(1) )
{
   PEL_LABEL( "PEL_ArithmeticExp:: PEL_ArithmeticExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_ArithmeticExp:: ~PEL_ArithmeticExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArithmeticExp:: ~PEL_ArithmeticExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
std::string const& 
PEL_ArithmeticExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   if( OP==M )
   {
      result = "IS|DS|SS " + name() + " <same type> " ;
   }
   else if( OP==MOD )
   {
      result = "modulo(IS,IS)" ;
   }
   else 
   {
      result = "IS|DS " + name() + " <same type> " ;
   }
       
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_ArithmeticExp:: valid_arguments(
                              PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArithmeticExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   bool result = some_arguments->count()==2 ;
   if( result )
   {
      PEL_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
      PEL_Data::Type k1 = extract_arg( some_arguments, 1 )->data_type() ;
      result = ( k0==k1 ) ;
      if( OP==MOD )
      {
         result &= ( k0 == PEL_Data::Int ) ;
      }
      else if( OP==M )
      {
         result &= ( k0 == PEL_Data::Int ||
                     k0 == PEL_Data::Double ||
                     k0 == PEL_Data::String )  ;
      }
      else
      {
         result &= ( k0 == PEL_Data::Int ||
                     k0 == PEL_Data::Double )  ;
      }
   }
   return result ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_ArithmeticExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   return arg(0)->data_type() ;
}

//----------------------------------------------------------------------
double
PEL_ArithmeticExp:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArithmeticExp:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE( ct ) ) ;
   
   double result = PEL::max_double() ;
   double v1 = ARG0->to_double( ct ) ;
   double v2 = ARG1->to_double( ct ) ;
   switch( OP )
   {
      case M : result = v1 + v2 ;
         break ;
      case L : result = v1 - v2 ;
         break ;
      case T : result = v1 * v2 ;
         break ;
      case D : result = v1 / v2 ;
         break ;
      default :
         PEL_Error::object()->raise_plain( "Internal error" ) ;
   }
   
   return result ;
}

//----------------------------------------------------------------------
int
PEL_ArithmeticExp:: to_int( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArithmeticExp:: to_int" ) ;
   PEL_CHECK_PRE( to_int_PRE(ct) ) ;
   
   int result = PEL::max_int() ;
   int v1 = ARG0->to_int( ct ) ;
   int v2 = ARG1->to_int( ct ) ;

   switch( OP )
   {
      case M : result = v1 + v2 ;
         break ;
      case L : result = v1 - v2 ;
         break ;
      case T : result = v1 * v2 ;
         break ;
      case D : result = v1 / v2 ;
         break ;
      case MOD : result = v1 % v2 ;
         break ;
      default :
         PEL_Error::object()->raise_plain( "Internal error" ) ;
   }
   
   return result ;
}

//----------------------------------------------------------------------
std::string const&
PEL_ArithmeticExp:: to_string( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArithmeticExp:: to_string" ) ;
   PEL_CHECK_PRE( to_string_PRE(ct) ) ;
   PEL_ASSERT( name()=="+" ) ;
   
   std::string const& v1 = ARG0->to_string( ct ) ;
   std::string const& v2 = ARG1->to_string( ct ) ;
   RESULT_STR = v1 + v2 ;
   return RESULT_STR ;
}

//----------------------------------------------------------------------
PEL_Data*
PEL_ArithmeticExp:: create_derivative( PEL_Object* a_owner,
                                       PEL_Variable const* var,
                                       PEL_Context const* ct  ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArithmeticExp:: create_derivative" ) ;
   PEL_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;
   
   PEL_Data* result = 0 ;
   PEL_List* list = PEL_List::create( 0 ) ;
   PEL_Data* d1 = ARG0->create_derivative( list, var, ct ) ;
   PEL_Data* d2 = ARG1->create_derivative( list, var, ct ) ;
   
   if( OP==M || OP==L )
   {
      list->append(d1) ;
      list->append(d2) ;
      
      result = PEL_Expression::create( a_owner, name(), list ) ;
      list->set_owner( result ) ;
   }
   else if( OP==T )
   {
      list->append(ARG0->create_clone(list)) ;
      list->append(d2) ; 
      PEL_Expression* m1 = PEL_Expression::create( 0, "*", list ) ;
      list->set_owner( m1 ) ;
      
      list = PEL_List::create( 0 ) ;
      list->append(d1) ;
      list->append(ARG1->create_clone(list)) ; 
      PEL_Expression* m2 = PEL_Expression::create( 0, "*", list ) ;
      list->set_owner( m2 ) ;
      
      list = PEL_List::create( 0 )  ;
      list->append(m1) ; m1->set_owner( list ) ;
      list->append(m2) ; m2->set_owner( list ) ;
      result = PEL_Expression::create( a_owner, "+", list ) ;
      list->set_owner( result ) ;
   }
   else if( OP==D )
   {
      list->append(d1) ;
      list->append(ARG1->create_clone(list)) ; 
      PEL_Expression* m1 = PEL_Expression::create( 0, "*", list ) ;
      list->set_owner( m1 ) ;
      
      list = PEL_List::create( 0 ) ;
      list->append(ARG0->create_clone(list)) ;
      list->append(d2) ; 
      PEL_Expression* m2 = PEL_Expression::create( 0, "*", list ) ;
      list->set_owner( m2 ) ;
      
      list = PEL_List::create( 0 ) ;
      list->append(m1) ; m1->set_owner( list ) ;
      list->append(m2) ; m2->set_owner( list ) ;
      PEL_Expression* num = PEL_Expression::create( 0, "-", list ) ;
      list->set_owner( num ) ;
      
      list = PEL_List::create( 0 ) ;
      list->append(ARG1->create_clone(list)) ;
      PEL_Expression* den = PEL_Expression::create( 0, "sqr", list ) ;
      list->set_owner( den ) ;

      list = PEL_List::create( 0 ) ;
      list->append(num) ; num->set_owner( list ) ;
      list->append(den) ; den->set_owner( list ) ;
      result = PEL_Expression::create( a_owner, "/", list ) ;
      list->set_owner( result ) ;
   }
   
   PEL_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_ArithmeticExp:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string space( indent_width, ' ' ) ;
   os << space ;
   if( external_brackets_are_set() ) os << "(" ;
   ARG0->print( os, 0 ) ;
   os << name() ;
   ARG1->print( os, 0 ) ;
   if( external_brackets_are_set() ) os << ")" ;
}

//----------------------------------------------------------------------
PEL_Data*
PEL_ArithmeticExp:: create_operator_simplification( PEL_Object* a_owner )  
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ArithmeticExp:: create_operator_simplification" ) ;

   PEL_Data* result = this ;
   
   bool first_const = ARG0->is_constant() ;
   bool second_const = ARG1->is_constant() ;
   PEL_CHECK( ! ( first_const && second_const ) ) ;
   if( first_const || second_const )
   {
      double v = ( first_const ? ARG0->to_double() : ARG1->to_double() ) ;
      PEL_Data const* non_const = ( first_const ? ARG1 : ARG0 ) ;
      bool null = v==0.0 ;
      bool unity = v==1.0 ;
      if( OP==M && null )
      {
         result = non_const->create_clone( a_owner ) ;
      }
      else if( OP==L && null && second_const )
      {
         result = non_const->create_clone( a_owner ) ;   
      }
      else if( OP==T && null )
      {
         result = PEL_Double::create( a_owner, 0.0 ) ;
      }
      else if( OP==T && unity )
      {
         result = non_const->create_clone( a_owner ) ;
      }
      else if( OP==D && null )
      {
         if( first_const )
         {
            result = PEL_Double::create( a_owner, 0.0 ) ;
         }
         else
         {
            PEL_Error::object()->raise_plain(
               "When simplifiing expression, null dividend found" ) ;
         }
      }
      else if( OP==D && unity &&  second_const )
      {
         
         result = non_const->create_clone( a_owner ) ;
      }
   }

   PEL_CHECK_POST( create_operator_simplification_POST( a_owner, result ) ) ;
   return result ;
}
