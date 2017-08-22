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

#include <PEL_SystemExp.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Module.hh>
#include <PEL_Sequence.hh>
#include <PEL_System.hh>

#include <iostream>
#include <sstream>


PEL_SystemExp const*
PEL_SystemExp::PROTOTYPE_PWD = new PEL_SystemExp( "getcwd" ) ;

PEL_SystemExp const* 
PEL_SystemExp::PROTOTYPE_GETENV = new PEL_SystemExp( "getenv" ) ;

PEL_SystemExp const* 
PEL_SystemExp::PROTOTYPE_JOIN = new PEL_SystemExp( "join" ) ;

PEL_SystemExp const* 
PEL_SystemExp::PROTOTYPE_SEPARATOR = new PEL_SystemExp( "path_name_separator" ) ;

PEL_SystemExp const* 
PEL_SystemExp::PROTOTYPE_DIRNAME = new PEL_SystemExp( "dirname" ) ;

PEL_SystemExp const* 
PEL_SystemExp::PROTOTYPE_BASENAME = new PEL_SystemExp( "basename" ) ;

PEL_SystemExp const* 
PEL_SystemExp::PROTOTYPE_GETPID = new PEL_SystemExp( "getpid" ) ;

PEL_SystemExp const* 
PEL_SystemExp::PROTOTYPE_UNAME = new PEL_SystemExp( "uname" ) ;

PEL_SystemExp const* 
PEL_SystemExp::PROTOTYPE_HOSTNAME = new PEL_SystemExp( "host_name" ) ;

//----------------------------------------------------------------------
PEL_SystemExp:: PEL_SystemExp( std::string const& a_name ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
{
   PEL_LABEL( "PEL_SystemExp:: PEL_SystemExp" ) ;
   PEL_CHECK( a_name=="getcwd" || a_name=="getenv"
              || a_name=="join" || a_name=="getpid"
              || a_name=="uname" || a_name=="host_name"
              || a_name=="path_name_separator"
              || a_name=="dirname" || a_name=="basename" ) ;
}

//----------------------------------------------------------------------
PEL_SystemExp:: PEL_SystemExp( PEL_Object* a_owner,
                               std::string const& a_name,
                               PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
{
   PEL_LABEL( "PEL_SystemExp:: PEL_SystemExp" ) ;
   PEL_CHECK( a_name=="getcwd" || a_name=="getenv"
              || a_name=="join" || a_name=="getpid"
              || a_name=="uname" || a_name=="host_name"
              || a_name=="path_name_separator"
              || a_name=="dirname" || a_name=="basename"  ) ;   
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_SystemExp:: ~PEL_SystemExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SystemExp:: ~PEL_SystemExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_SystemExp*
PEL_SystemExp:: create_replica( PEL_Object* a_owner,
                                PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SystemExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_SystemExp* result = new PEL_SystemExp( a_owner, name(), argument_list ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_SystemExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SystemExp:: data_type" ) ;
   PEL_Data::Type result = String ;
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_SystemExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SystemExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   bool result = false ;
   if( name()=="getcwd" || name()=="getpid" || name()=="uname" 
       || name()=="host_name" || name()=="path_name_separator" )
   {
      result = some_arguments->count()==0 ;
   }
   else if( name()=="getenv" || name()=="dirname" || name()=="basename" )
   {
      result = some_arguments->count()==1 &&
         extract_arg( some_arguments, 0 )->data_type() == String ;
   }
   else 
   {
      PEL_ASSERT( name()=="join" ) ;
      result = some_arguments->count()>1 ;
      for( size_t i=0 ; i<some_arguments->index_limit() ; i++ )
      {
         result = result &&
            extract_arg( some_arguments, i )->data_type() == String ;
      }
   }

   return result ;
}

//----------------------------------------------------------------------
std::string const&
PEL_SystemExp:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SystemExp:: usage" ) ;
   static std::string result ;
   if( name()=="getcwd" || name()=="getpid" || name()=="path_name_separator" ||
       name()=="uname" || name()=="host_name" )
   {
      result = name()+"()" ;
   }
   else if( name()=="getenv" || name()=="dirname" || name()=="basename" )
   {
      result = name() + "(SS)" ;
   }
   else if( name()=="join" )
   {
      result = "join(<list of SS>)" ;
   }
   return result ;
}

//----------------------------------------------------------------------
std::string
const&
PEL_SystemExp:: to_string( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SystemExp:: to_string" ) ;
   PEL_CHECK_PRE( to_string_PRE(ct) ) ;
   
   static char sep = PEL_System::path_name_separator() ;
   
   RESULT_STR = "" ;
  
   if( name()=="getcwd" )
   {
      RESULT_STR = PEL_System::working_directory() ;
   }
   else if( name()=="getpid" )
   {
      std::ostringstream os ;
      os << PEL_System::process_id() ;
      RESULT_STR = os.str() ;
   }
   else if( name()=="getenv" )
   {
      RESULT_STR = PEL_System::getenv( arg(0)->to_string( ct ) ) ;
   }
   else if( name()=="dirname" )
   {
      RESULT_STR = PEL_System::dirname( arg(0)->to_string( ct ) ) ;
   }
   else if( name()=="basename" )
   {
      RESULT_STR = PEL_System::basename( arg(0)->to_string( ct ) ) ;
   }
   else if( name()=="path_name_separator" )
   {
      RESULT_STR = PEL_System::path_name_separator() ;
   }
   else if( name()=="join" )
   {
      RESULT_STR = "" ;
      for( size_t i=0 ; i<nb_arguments() ; i++ )
      {
         std::string const& item = arg(i)->to_string( ct ) ;
         if( i!=0 && item.length()>0 && item[0]!=sep ) RESULT_STR += sep ;
         RESULT_STR += item ;
      }
   }
   else if( name()=="uname" )
   {
      RESULT_STR = PEL_System::sysname() ;
   }
   else if( name()=="host_name" )
   {
      RESULT_STR = PEL_System::host_name() ;
   }
   
   return RESULT_STR ;   
}
