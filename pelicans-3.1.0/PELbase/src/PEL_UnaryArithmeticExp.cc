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

#include <PEL_UnaryArithmeticExp.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <iostream>

//----------------------------------------------------------------------
PEL_UnaryArithmeticExp const* 
PEL_UnaryArithmeticExp:: PROTOTYPE_Minus = 
                          new PEL_UnaryArithmeticExp( "unary_minus", Minus ) ;
//----------------------------------------------------------------------

//----------------------------------------------------------------------
PEL_UnaryArithmeticExp:: PEL_UnaryArithmeticExp( std::string const& a_name,
                                                 UnaryOperator a_op ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , OP( a_op )
   , FIRST( 0 )
{
   PEL_LABEL( "PEL_UnaryArithmeticExp:: PEL_UnaryArithmeticExp" ) ;
}

//----------------------------------------------------------------------
PEL_UnaryArithmeticExp*
PEL_UnaryArithmeticExp:: create_replica( PEL_Object* a_owner,
                             PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_UnaryArithmeticExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_UnaryArithmeticExp* result = new PEL_UnaryArithmeticExp( a_owner, 
                                                          name(), 
                                                          argument_list, 
                                                          OP ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_UnaryArithmeticExp:: PEL_UnaryArithmeticExp( PEL_Object* a_owner,
                                     std::string const& a_name,
                                     PEL_Sequence const* argument_list,
                                     UnaryOperator a_op ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , OP( a_op )
   , FIRST( arg(0) )
{
   PEL_LABEL( "PEL_UnaryArithmeticExp:: PEL_UnaryArithmeticExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_UnaryArithmeticExp:: ~PEL_UnaryArithmeticExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_UnaryArithmeticExp:: ~PEL_UnaryArithmeticExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_UnaryArithmeticExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_UnaryArithmeticExp:: data_type" ) ;
   
   return FIRST->data_type() ;
}

//----------------------------------------------------------------------
std::string const& 
PEL_UnaryArithmeticExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result = name() + "<double|integer>" ;
   
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_UnaryArithmeticExp:: valid_arguments(
                              PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_UnaryArithmeticExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = some_arguments->count()==1 ;
   if( result )
   {
      PEL_Data::Type k0 =  extract_arg( some_arguments, 0 )->data_type() ;
      result = result && ( k0==Double || k0==Int ) ;
   }
   return result ;
}

//----------------------------------------------------------------------
double
PEL_UnaryArithmeticExp:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_UnaryArithmeticExp:: to_double" ) ;
   
   PEL_CHECK_PRE( to_double_PRE(ct) ) ;
   double result = PEL::max_double() ;
   double v = FIRST->to_double(ct) ;
   switch(OP)
   {
      case Minus : result = -v ;
         break ;
      default :
         PEL_Error::object()->raise_plain( "Internal error" ) ;
   }
   
   return result ;
}

//----------------------------------------------------------------------
int
PEL_UnaryArithmeticExp:: to_int( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_UnaryArithmeticExp:: to_int" ) ;
   
   PEL_CHECK_PRE( to_int_PRE(ct) ) ;
   int result = PEL::max_int() ;
   int v = FIRST->to_int(ct) ;
   switch(OP)
   {
      case Minus : result = -v ;
         break ;
      default :
         PEL_Error::object()->raise_plain( "Internal error" ) ;
   }
   
   return result ;
}

//----------------------------------------------------------------------
PEL_Data*
PEL_UnaryArithmeticExp:: create_derivative( PEL_Object* a_owner,
                                            PEL_Variable const* var,
                                            PEL_Context const* ct  ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_UnaryArithmeticExp:: create_derivative" ) ;
   PEL_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;

   PEL_List* list = PEL_List::create( 0 ) ;
   list->append( FIRST->create_derivative( list, var, ct ) ) ;
   
   PEL_Data* result = create_replica( a_owner, list ) ;
   list->set_owner( result ) ;
   
   PEL_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_UnaryArithmeticExp:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string space( indent_width, ' ' ) ;
   os << space << "-" ;
   if( external_brackets_are_set() ) os << "(" ;
   FIRST->print( os, 0 ) ;
   if( external_brackets_are_set() ) os << ")" ;
}
