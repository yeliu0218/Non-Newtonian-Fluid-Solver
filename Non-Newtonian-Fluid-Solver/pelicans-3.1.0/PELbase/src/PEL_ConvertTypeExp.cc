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

#include <PEL_ConvertTypeExp.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <iostream>
#include <cmath>

PEL_ConvertTypeExp const* 
PEL_ConvertTypeExp::PROTOTYPE_DOUBLE = new PEL_ConvertTypeExp( "double" ) ;

PEL_ConvertTypeExp const* 
PEL_ConvertTypeExp::PROTOTYPE_INT = new PEL_ConvertTypeExp( "int" ) ;


//----------------------------------------------------------------------
PEL_ConvertTypeExp:: PEL_ConvertTypeExp( std::string const& a_name ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , P0( 0 )
{
   PEL_LABEL( "PEL_ConvertTypeExp:: PEL_ConvertTypeExp" ) ;
}

//----------------------------------------------------------------------
PEL_ConvertTypeExp*
PEL_ConvertTypeExp:: create_replica( PEL_Object* a_owner,
                             PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConvertTypeExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_ConvertTypeExp* result = new PEL_ConvertTypeExp( a_owner,
                                                        name(),
                                                        argument_list ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ConvertTypeExp:: PEL_ConvertTypeExp( PEL_Object* a_owner,
                                         std::string const& a_name,
                                         PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , P0( arg(0) )
{
   PEL_LABEL( "PEL_ConvertTypeExp:: PEL_ConvertTypeExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_ConvertTypeExp:: ~PEL_ConvertTypeExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConvertTypeExp:: ~PEL_ConvertTypeExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_ConvertTypeExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConvertTypeExp:: data_type" ) ;
   
   return ( name()=="double" ? PEL_Data::Double : PEL_Data::Int ) ;
}

//----------------------------------------------------------------------
std::string const& 
PEL_ConvertTypeExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   result = ( (name()=="double") ? "double(IS)" : "int(DS)" ) ;
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_ConvertTypeExp:: valid_arguments(
                              PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConvertTypeExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   bool result = some_arguments->count()==1 ;
   if( result )
   {
      PEL_Data::Type k0 =  extract_arg( some_arguments, 0 )->data_type() ;
      result = result && ( ( name() == "double" && k0 == Int     ) ||
                           ( name() == "int"    && k0 == Double  ) ) ;
   }
   return result ;
}

//----------------------------------------------------------------------
double
PEL_ConvertTypeExp:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConvertTypeExp:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE(ct) ) ;

   double result = (double) P0->to_int( ct ) ;

   return result ;
}

//----------------------------------------------------------------------
int
PEL_ConvertTypeExp:: to_int( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_ConvertTypeExp:: to_int" ) ;
   PEL_CHECK_PRE( to_int_PRE(ct) ) ;

   int result = (int) P0->to_double(ct) ;

   return result ;
}
