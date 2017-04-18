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

#include <GE_ColorExp.hh>

#include <GE_Color.hh>
#include <PEL_assertions.hh>
#include <PEL_Sequence.hh>
#include <doubleVector.hh>
#include <intVector.hh>

GE_ColorExp const*  
GE_ColorExp::PROTOTYPE_matching = 
                             new GE_ColorExp( matching, "matching_color" ) ;

GE_ColorExp const*  
GE_ColorExp::PROTOTYPE_valid = 
                             new GE_ColorExp( valid, "valid_color" ) ;

//----------------------------------------------------------------------
GE_ColorExp:: GE_ColorExp( MembExp exp_id, 
                           std::string const& a_name  ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , OP( exp_id )
   , ARG0( 0 )
   , ARG1( 0 )
{
   PEL_LABEL( "GE_ColorExp:: GE_ColorExp" ) ;
}

//----------------------------------------------------------------------
GE_ColorExp*
GE_ColorExp:: create_replica( PEL_Object* a_owner,
                              PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_ColorExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   GE_ColorExp* result = new GE_ColorExp( a_owner, 
                                          OP,
                                          name(),
                                          argument_list ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_ColorExp:: GE_ColorExp( PEL_Object* a_owner,
                           MembExp exp_id,
                           std::string const& a_name,
                           PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list ),
     OP( exp_id ),
     ARG0( arg( 0 ) ),
     ARG1( exp_id == matching ? arg( 1 ) : 0 )
{
   PEL_LABEL( "GE_ColorExp:: GE_ColorExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_ColorExp:: ~GE_ColorExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_ColorExp:: ~GE_ColorExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
std::string const& 
GE_ColorExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   switch( OP )
   {
      case valid :
         result =
           "valid_color(SS)" ;
	 break ;
      case matching :
         result = 
           "matching_color(SS,SS)" ;
	 break ;
   }
   return result ;
}

//----------------------------------------------------------------------
PEL_Data::Type
GE_ColorExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   return PEL_Data::Bool ;
}

//----------------------------------------------------------------------
bool
GE_ColorExp:: to_bool( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_ColorExp:: to_bool" ) ;
   PEL_CHECK_PRE( to_bool_PRE( ct ) ) ;

   bool result = false ;
   std::string color_name = ARG0->to_string(ct) ;
   switch( OP )
   {
      case valid :
         result = GE_Color::exist( color_name ) ;
         break ;
      case matching :
         GE_Color const* col1 = GE_Color::object( color_name ) ;
         GE_Color const* col2 = GE_Color::object( ARG1->to_string(ct) ) ;
         result = col1->is_matching( col2 ) ;
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
bool
GE_ColorExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_ColorExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = false ;
   switch( OP )
   {
      case valid :
         result = some_arguments->count()==1 &&
            extract_arg( some_arguments, 0 )->data_type() == String ;
	 break ;
      case matching :
         result = ( some_arguments->count() == 2 ) &&
            extract_arg( some_arguments, 0 )->data_type() == String &&
            extract_arg( some_arguments, 1 )->data_type() == String ;
	 break ;
      default : result = false ;
   }
   return result ;
}
