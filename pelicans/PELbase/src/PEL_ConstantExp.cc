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

#include <PEL_ConstantExp.hh>

#include <PEL.hh>

#include <PEL_assertions.hh>
#include <PEL_Double.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>

PEL_ConstantExp const* 
PEL_ConstantExp::PROTOTYPE_pi = new PEL_ConstantExp( "pi", PEL::pi() ) ;

PEL_ConstantExp const* 
PEL_ConstantExp::PROTOTYPE_e = new PEL_ConstantExp( "e", PEL::e() ) ;

PEL_ConstantExp const* 
PEL_ConstantExp::PROTOTYPE_eu = new PEL_ConstantExp( "euler", PEL::euler() ) ;

//----------------------------------------------------------------------
PEL_ConstantExp:: PEL_ConstantExp( std::string const& a_name,
                                   double a_val) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , VAL( a_val )
{
   PEL_LABEL( "PEL_ConstantExp:: PEL_ConstantExp" ) ;
}

//----------------------------------------------------------------------
PEL_ConstantExp*
PEL_ConstantExp:: create_replica( PEL_Object* a_owner,
                                  PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConstantExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_ConstantExp* result =
                  new PEL_ConstantExp( a_owner, name(), VAL, argument_list ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ConstantExp:: PEL_ConstantExp( PEL_Object* a_owner,
                                   std::string const& a_name,
                                   double a_val,
                                   PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , VAL( a_val )
{
   PEL_LABEL( "PEL_ConstantExp:: PEL_ConstantExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_ConstantExp:: ~PEL_ConstantExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConstantExp:: ~PEL_ConstantExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_ConstantExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConstantExp:: data_type" ) ;
   PEL_Data::Type result = Double ;
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_ConstantExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConstantExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   bool result = some_arguments->count()==0 ;
   return result ;
}

//----------------------------------------------------------------------
std::string const&
PEL_ConstantExp:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConstantExp:: usage" ) ;
   static std::string result ;
   result = name()+"()" ;
   return result ;
}

//----------------------------------------------------------------------
double
PEL_ConstantExp:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConstantExp:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE(ct) ) ;
   return VAL ;
}

//----------------------------------------------------------------------
PEL_Data*
PEL_ConstantExp:: create_derivative( PEL_Object* a_owner,
                                     PEL_Variable const* var,
                                     PEL_Context const* ct ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConstantExp:: create_derivative" ) ;
   PEL_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;
   
   PEL_Data* result = PEL_Double::create( a_owner, 0.0 ) ;
   
   PEL_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

