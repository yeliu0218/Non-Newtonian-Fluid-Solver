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

#include <PEL_DerivativeExp.hh>

#include <PEL_assertions.hh>
#include <PEL_Double.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL_Variable.hh>
#include <PEL.hh>

#include <iostream>

PEL_DerivativeExp const* 
PEL_DerivativeExp:: PROTOTYPE_d = new PEL_DerivativeExp( "d", d ) ;

PEL_DerivativeExp const* 
PEL_DerivativeExp:: PROTOTYPE_dnum = new PEL_DerivativeExp( "dnum", dnum ) ;

//----------------------------------------------------------------------
PEL_DerivativeExp:: PEL_DerivativeExp( std::string const& a_name,
                                       OP_TYPE a_op ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , OP( a_op )
   , EXP( 0 )
   , VAR( 0 )
   , DERIVATIVE( 0 )
{
   PEL_LABEL( "PEL_DerivativeExp:: PEL_DerivativeExp" ) ;
}

//----------------------------------------------------------------------
PEL_DerivativeExp*
PEL_DerivativeExp:: create_replica( PEL_Object* a_owner,
                                    PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DerivativeExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_DerivativeExp* result = new PEL_DerivativeExp( a_owner,
                                                      name(),
                                                      OP,
                                                      argument_list ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_DerivativeExp:: PEL_DerivativeExp( PEL_Object* a_owner,
                                       std::string const& a_name,
                                       OP_TYPE a_op,
                                       PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , OP( a_op )
   , EXP( arg(0) )
   , VAR( PEL_Variable::object( arg(1)->to_string() ) )
   , DERIVATIVE( 0 )
{
   PEL_LABEL( "PEL_DerivativeExp:: PEL_DerivativeExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_DerivativeExp:: ~PEL_DerivativeExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DerivativeExp:: ~PEL_DerivativeExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_DerivativeExp*
PEL_DerivativeExp:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DerivativeExp:: create_clone" ) ;
   PEL_CHECK_PRE( !is_a_prototype() ) ;
   
   PEL_DerivativeExp* result =
      static_cast<PEL_DerivativeExp*>(
         PEL_Expression::create_clone( a_owner ) ) ;

   if( !is_a_prototype() && DERIVATIVE!=0 )
      result->DERIVATIVE = DERIVATIVE->create_clone( result ) ;
   
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
PEL_DerivativeExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   result = name() + "(<expression>,SS)" ;
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_DerivativeExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DerivativeExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   bool result = some_arguments->count()==2 ;
   if( result )
   {
      PEL_Data::Type k0 =  extract_arg( some_arguments, 0 )->data_type() ;
      result = result && ( k0==Double || k0==DoubleVector ) ;
      PEL_Data::Type k1 =  extract_arg( some_arguments, 1 )->data_type() ;
      result = result && ( k1==String ) ;
   }
   return result ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_DerivativeExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DerivativeExp:: data_type" ) ;
   
   return EXP->data_type() ;
}

//----------------------------------------------------------------------
double
PEL_DerivativeExp:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DerivativeExp:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE(ct) ) ;

   double result ;
   
   if( OP==d )
   {
      result = derivative(ct)->to_double(ct) ;
   }
   else
   {
      PEL_ASSERT( OP==dnum ) ;
      static double const eps = 1.0e-8 ;
      
      double val = ct->value( VAR )->to_double(ct) ;
      PEL_ContextSimple* ctx = PEL_ContextSimple::create( 0 ) ;
      ctx->extend( ct ) ;
      double dx = PEL::max( eps, val*eps ) ;
      
      ctx->set_value_of( VAR, PEL_Double::create( ctx, val+dx ) ) ;
      result = ( EXP->to_double( ctx ) - EXP->to_double( ct ) )/dx ;
      ctx->destroy() ; ctx=0 ;
   }
   return result ;
}

//----------------------------------------------------------------------
doubleVector const&
PEL_DerivativeExp:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DerivativeExp:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ;

   PEL_ASSERT( OP==d ) ;

   doubleVector const& result = derivative(ct)->to_double_vector(ct) ;
   return result ;
}

//----------------------------------------------------------------------
PEL_Data*
PEL_DerivativeExp:: create_derivative( PEL_Object* a_owner,
                                       PEL_Variable const* var,
                                       PEL_Context const* ct ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DerivativeExp:: create_derivative" ) ;
   PEL_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;

   if( OP!=d )
   {
      PEL_Error::object()->raise_plain( "Unable to differentiate "+name()+" operator" ) ;
   }
   
   PEL_Data* result = derivative(ct)->create_derivative( a_owner, var, ct ) ;
   
   PEL_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_DerivativeExp:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string space( indent_width, ' ' ) ;
   os << space << name() << "( " ;
   EXP->print( os, 0 ) ;
   os << ", \"" << VAR->name() << "\" )" ;
}

//----------------------------------------------------------------------
PEL_Data*
PEL_DerivativeExp:: create_non_const_simplification(
                                              PEL_Object* a_owner ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DerivativeExp:: create_non_const_simplification" ) ;
   PEL_CHECK( create_non_const_simplification_PRE( a_owner ) ) ;
   PEL_Data* result = 0 ;
   
   if( OP==d )
   {
      PEL_ASSERT( DERIVATIVE!=0 ) ;
      result = DERIVATIVE->create_simplification( a_owner ) ;
   }
   
   PEL_CHECK( create_non_const_simplification_POST( a_owner, result ) ) ;
   return result ;
}

//----------------------------------------------------------------------
PEL_Data const*
PEL_DerivativeExp:: derivative( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DerivativeExp:: derivative" ) ;
   PEL_CHECK( OP==d ) ;
   
   if( DERIVATIVE==0 )
   {
      PEL_Data* der = EXP->create_derivative( 0, VAR, ct ) ;
//       std::cout << "Derivative of " ;
//       EX->print( std::cout, 1 ) ;
//       std::cout << std::endl << " with respect to "  ;
//       VAR->print( std::cout, 1 ) ;
//       std::cout << std::endl << " is : " ;
//       der->print( std::cout, 0 ) ;
//       std::cout << std::endl << " Simplified in " ;
      DERIVATIVE = der->create_simplification(
         const_cast<PEL_DerivativeExp*>(this) ) ;
//       DERIVATIVE->print( std::cout, 0 ) ;
//      std::cout << std::endl ;
      
      der->destroy() ;
   }
   return DERIVATIVE ;
}
