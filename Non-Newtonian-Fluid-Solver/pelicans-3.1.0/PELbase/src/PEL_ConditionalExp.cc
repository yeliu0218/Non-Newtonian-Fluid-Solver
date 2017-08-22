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

#include <PEL_ConditionalExp.hh>

#include <PEL_assertions.hh>
#include <PEL_List.hh>
#include <PEL_Vector.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL_String.hh>

#include <boolVector.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <iostream>

PEL_ConditionalExp const*
PEL_ConditionalExp::PROTOTYPE = new PEL_ConditionalExp() ;

//----------------------------------------------------------------------
PEL_ConditionalExp:: PEL_ConditionalExp( void ) 
//----------------------------------------------------------------------
      : PEL_TransferExp( "(?:)" )
      , DEFAULT(0)
{
   PEL_LABEL( "PEL_ConditionalExp:: PEL_ConditionalExp" ) ;
}

//----------------------------------------------------------------------
PEL_ConditionalExp:: PEL_ConditionalExp( PEL_Object* a_owner,
                                         PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
      : PEL_TransferExp( a_owner, "(?:)", argument_list )
      , DEFAULT( arg( nb_arguments()-1 ) )
{
   PEL_LABEL( "PEL_ConditionalExp:: PEL_ConditionalExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_ConditionalExp:: ~PEL_ConditionalExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConditionalExp:: ~PEL_ConditionalExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_ConditionalExp*
PEL_ConditionalExp:: create_replica( PEL_Object* a_owner,
                             PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConditionalExp:: create_replica" ) ;
   PEL_CHECK_PRE( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_ConditionalExp* result = new PEL_ConditionalExp( a_owner, 
                                                      argument_list ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Data const*
PEL_ConditionalExp:: data( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConditionalExp:: data" ) ;
   PEL_CHECK( data_PRE( ct ) ) ;
   
   PEL_Data const* result = 0 ;
   size_t nb_tests = (nb_arguments()-1)/2 ;
   for( size_t i=0 ; i<nb_tests ; i++ )
   {
      size_t idx = 2*i ;
      if( arg(idx)->to_bool(ct) ) result=arg(idx+1) ;
   }
   if( result==0 ) result=DEFAULT ;
   
   PEL_CHECK( data_POST( result, ct ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Data*
PEL_ConditionalExp:: create_derivative( PEL_Object* a_owner,
                                        PEL_Variable const* var,
                                        PEL_Context const* ct  ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConditionalExp:: create_derivative" ) ;
   PEL_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;
   
   PEL_Data* result = 0 ;
   PEL_List* list = PEL_List::create( 0 ) ;
   size_t nb_tests = (nb_arguments()-1)/2 ;
   for( size_t i=0 ; i<nb_tests ; i++ )
   {
      size_t idx = 2*i ;
      list->append( arg(idx)->create_clone( list ) ) ;
      list->append( arg(idx+1)->create_derivative( list, var, ct ) ) ;
   }
   list->append( DEFAULT->create_derivative( list, var, ct ) ) ;
      
   result = PEL_Expression::create( a_owner, name(), list ) ;
   list->set_owner( result ) ;
   
   PEL_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
PEL_ConditionalExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result =
      "( BS ? <val1> : ... BS ? <valN> : <defaultVal> )" ;
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_ConditionalExp:: valid_arguments(
                              PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConditionalExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = some_arguments->count()>=3 &&
      ( some_arguments->count()%2 ) == 1;
   if( result )
   {
      PEL_Data::Type k_default =
         extract_arg( some_arguments, some_arguments->count()-1 )->data_type() ;
      size_t nb_tests = (some_arguments->count()-1)/2 ;
      for( size_t i=0 ; i<nb_tests ; i++ )
      {
         size_t idx = 2*i ;
         PEL_Data::Type k0 =  extract_arg( some_arguments, idx )->data_type() ;
         PEL_Data::Type k1 =  extract_arg( some_arguments, idx+1 )->data_type() ;
         result = result && k0==Bool && k1==k_default ;
      }
   }
   return result ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_ConditionalExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   Type result = DEFAULT->data_type() ;
   return result ;
}

//----------------------------------------------------------------------
void
PEL_ConditionalExp:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string space( indent_width, ' ' ) ;
   os << space << "( " ;
   size_t nb_tests = (nb_arguments()-1)/2 ;
   for( size_t i=0 ; i<nb_tests ; i++ )
   {
      size_t idx = 2*i ;
      arg(idx)->print(os,0) ;
      os << " ? " ;
      arg(idx+1)->print(os,0) ;
      os << " : " ;
   }
   DEFAULT->print( os, 0 ) ;
   os << " ) " ;
}
