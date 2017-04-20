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

#include <PEL_Error.hh>

#include <PEL_Exec.hh>
#include <PEL_Exceptions.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Int.hh>
#include <PEL_Map.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>
#include <PEL_String.hh>
#include <PEL_StringVector.hh>
#include <PEL_assertions.hh>
#include <stringVector.hh>

#include <iostream>
#include <sstream>

using std::endl ;
using std::string ;

//---------------------------------------------------------------------------
PEL_Error*
PEL_Error:: object( void )
//---------------------------------------------------------------------------
{
   static PEL_Error UNIQUE_INSTANCE ;
   return( &UNIQUE_INSTANCE ) ;
}

//---------------------------------------------------------------------------
PEL_Error:: PEL_Error( void )
//---------------------------------------------------------------------------
   : os( PEL_Exec::out() )
{
}

//---------------------------------------------------------------------------
PEL_Error:: ~PEL_Error( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_internal( std::string message )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "Internal Error" << endl << endl ;
   os << message << endl << endl ;
   print_invocated_methods() ;
   end() ;
   exit_with_internal_error() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_plain( std::string message )
//---------------------------------------------------------------------------
{
   begin() ;
   os << message << endl << endl ;
   end() ;
   // print_invocated_methods() ;
   exit() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: display_info( std::string message )
//---------------------------------------------------------------------------
{
   os << endl << hline() << endl ;
   os << message << endl ;
   os << hline() << endl ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_not_tested( std::string file, int line, std::string test )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "Situation never tested : " << test << endl << endl ;
   os << "in file : " << file << endl ;
   os << "at line : " << line << endl << endl ;
   print_invocated_methods() ;
   end() ;
   exit_with_internal_error() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_not_implemented( PEL_Object const* oo,
                                   std::string method )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "Call to not implemented method : " << method << endl ;
   os << "         for an object of type : " << oo->type_name()
      << endl << endl ;
   print_invocated_methods() ;
   end() ;
   exit_with_internal_error() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_bad_types( PEL_Object const* oo,
                             std::string method,
                             PEL_Object const* argument )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "Uncovered combination of run-time types" << endl ;
   os << "in call to method     : " << method << endl ;
   os << "for an object of type : " << oo->type_name() << endl ;
   os << "with argument attached to object of type : " << endl ;
   os << "   " << argument->type_name() << endl << endl ;
   print_invocated_methods() ;
   end() ;
   exit_with_internal_error() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_bad_types( PEL_Object const* oo,
                             std::string method,
                             PEL_Object const* argument_1,
                             PEL_Object const* argument_2 )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "Uncovered combination of run-time types" << endl ;
   os << "in call to method     : " << method << endl ;
   os << "for an object of type : " << oo->type_name() << endl ;
   os << "with arguments attached to objects of type : " << endl ;
   os << "   " << argument_1->type_name() << endl ;
   os << "   " << argument_2->type_name() << endl << endl ;
   print_invocated_methods() ;
   end() ;
   exit_with_internal_error() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_read_syntax_error( std::string const& file,
                                    int line,
                                    std::string const& last,
                                    std::string const& nature )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "Error while reading data file." << endl << endl ;
   os << nature << endl << endl ;
   os << "Last line number " << line << " read in file " << file << " was :" << endl ;
   os << ">> " << last << " <<" << endl ;
   end() ;
   exit() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_missing_keyword( PEL_ModuleExplorer const* exp,
                                   std::string keyword )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "In module: " << endl << "   " << exp->absolute_path_name() << endl ;
   os << "the following keyword is requested: " << endl
      << "   \"" << keyword << "\"" << endl ;
   end() ;
   exit() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_missing_module( PEL_ModuleExplorer const* exp,
                                  std::string path_and_name )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "In module: " << endl << "   " << exp->absolute_path_name() << endl ;
   os << "the following MODULE is requested: " << endl
      << "   \"" << path_and_name << "\"" << endl ;
   end() ;
   exit() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_module_error( PEL_ModuleExplorer const* exp,
                                std::string error )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "In module: " << endl << "   " << exp->absolute_path_name() << endl ;
   os << error << endl ;
   end() ;
   exit() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_bad_data_type( PEL_ModuleExplorer const* exp,
                                 std::string keyword,
                                 PEL_Data::Type query_kind )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "In module : " << exp->absolute_path_name() << endl ;
   os << "the data of keyword : " << keyword << endl ;
   os << "should be of type : "
      << PEL_Data::type_name(query_kind) << endl ;
   end() ;
   exit() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_bad_file( PEL_ModuleExplorer const* exp,
                            std::string const& filename,
                            std::string const& access )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "In module : " << exp->absolute_path_name() << endl ;
   os << " file : " << exp->string_data( filename ) << endl ;

   os << " specified by the data of keyword : " << filename << endl ;

   os << " can't be accessed in mode : " << access << endl ;
   end() ;
   exit() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_not_evaluable( PEL_ModuleExplorer const* exp,
                                 std::string keyword,
                                 stringVector const& undefined_variables )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "In module : " << exp->absolute_path_name() << endl ;
   os << "the data of keyword : " << keyword << endl ;
   os << "cannot be evaluated " << endl ;
   if( undefined_variables.size() > 0 )
   {
      os << "  undefined variable(s): " << endl ;
      for( size_t i=0 ; i<undefined_variables.size() ; ++i )
      {
         os << "      - \"" << undefined_variables(i) << "\"" << endl ;
      }
   }
   end() ;
   exit() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_bad_data_value( PEL_ModuleExplorer const* exp,
                                  std::string keyword,
                                  std::string allowed_values )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "In module: " << endl << "   " << exp->absolute_path_name() << endl ;
   os << "the data of keyword: " << keyword << endl ;
   os << "should have one of the following values: " << endl ;
   os << allowed_values << endl ;
   end() ;
   exit() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_data_error( PEL_ModuleExplorer const* exp,
                              std::string keyword,
                              std::string error )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "In module: " << endl << "   " << exp->absolute_path_name() << endl ;
   os << "error in the data of keyword: " << keyword << endl ;
   os << error << endl ;
   end() ;
   exit() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_bad_object( std::string message,
                              PEL_Object const* a_object )
//---------------------------------------------------------------------------
{
   begin() ;
   os << message << endl ;
   a_object->print( os, 1 ) ;
   end() ;
   exit() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_precondition_violation( std::string test )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "PRECONDITION assertion violation :" << endl ;
   os << "   " << test << endl << endl ;
   print_invocated_methods() ;
   os << endl ;
   print_client_and_server() ;
   os << endl ;
   os << "HINTS : " << endl ;
   os << "   - Meeting a PRECONDITION is the CLIENT\'s responsibility."
      << endl ;
   os << "   - A PRECONDITION violation is the manifestation of " << endl
      << "     a bug in the CLIENT." << endl ;
   end() ;
   exit_with_internal_error() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_postcondition_violation( std::string file,
                                           int line,
                                           std::string test )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "POSTCONDITION assertion violation :" << endl ;
   os << "   " << test << endl ;
   os << "in file : " << file << endl ;
   os << "at line : " << line << endl << endl ;
   print_invocated_methods() ;
   os << endl ;
   print_client_and_server() ;
   os << endl ;
   os << "HINTS : " << endl ;
   os << "   - Meeting a POSTCONDITION is the SUPPLIER\'s responsibility."
      << endl ;
   os << "   - A POSTCONDITION violation is the manifestation of " << endl
      << "     a bug in the SUPPLIER." << endl ;
   end() ;
   exit_with_internal_error() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_invariant_violation( std::string file,
                                       int line,
                                       std::string test )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "INVARIANT assertion violation :" << endl ;
   os << "   " << test << endl ;
   os << "in file : " << file << endl ;
   os << "at line : " << line << endl << endl ;
   print_invocated_methods() ;
   end() ;
   exit_with_internal_error() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_assertion_violation( std::string file,
                                       int line,
                                       std::string test )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "ASSERTION violation :" << endl ;
   os << "   " << test << endl  ;
   os << "in file : " << file << endl ;
   os << "at line : " << line << endl << endl ;
   print_invocated_methods() ;
   end() ;
   exit_with_internal_error() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: raise_file_handling( std::string file,
                                 std::string operation )
//---------------------------------------------------------------------------
{
   begin() ;
   os << "Failure of operation : " << operation << endl ;
   os << "for file : " << file << endl << endl ;
   print_invocated_methods() ;
   end() ;
   exit_with_internal_error() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: begin( void )
//---------------------------------------------------------------------------
{
   os << endl << endl << hline() << endl ;
   os << "              FATAL ERROR" << endl << hline() << endl ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: end( void )
//---------------------------------------------------------------------------
{
   os << hline() << endl ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: print_invocated_methods( void )
//---------------------------------------------------------------------------
{
   if( PEL_Marker::nb_labels() > 0 )
   {
      os << "Stack of invocated methods :" << endl ;
      for( size_t i=PEL_Marker::nb_labels()-1 ; i<PEL_Marker::nb_labels() ; i-- )
      {
         os << "   " << PEL_Marker::label(i) << endl ;
      }
      os << "WARNING : only the methods labelled with \'PEL_LABEL\' "
         << "are listed." << endl ;
   }
   else
   {
      os << "Stack of invocated methods unavailable" << endl ;
      os << "(none of them have been labelled with \'PEL_LABEL\')."
         << endl ;
   }
}

//---------------------------------------------------------------------------
void
PEL_Error:: print_client_and_server( void )
//---------------------------------------------------------------------------
{
   if( PEL_Marker::nb_labels() > 1 )
   {
      os << "CLIENT   : " << PEL_Marker::label( PEL_Marker::nb_labels()-2 ) << endl ;
   }
   if( PEL_Marker::nb_labels() > 0 )
   {
      os << "SUPPLIER : " << PEL_Marker::label( PEL_Marker::nb_labels()-1 ) << endl ;
      os << "(inferred from the above stack of invocated methods)."
         << endl ;
   }
}

//---------------------------------------------------------------------------
std::string const&
PEL_Error:: hline( void ) const
//---------------------------------------------------------------------------
{
   static std::string const h = string( 50, '-' ) ;
   return( h ) ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: trace( std::string const& event )
//---------------------------------------------------------------------------
{
   os << endl << hline() << endl ;
   os << " EVENT TRACKING " ;
   os << event << endl << hline() << endl ;
   print_invocated_methods() ;
   os << hline() << endl << endl ;
}

//-------------------------------------------------------------------------
void
PEL_Error:: display_new_syntax( PEL_ModuleExplorer const* older_module,
                                PEL_ModuleExplorer const* result_module ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Error:: display_new_syntax" ) ;

   os << " ********* PELICANS information " << endl ;
   os << " Syntax for module : " << endl ;
   os << "<==============================" << endl ;
   older_module->print( os, 0 ) ;
   os << "<==============================" << endl ;
   os << " won't be supported anymore. " << endl ;
   os << " Please let replace with following :" << endl ;
   os << ">==============================" << endl ;
   result_module->print( os, 0 ) ;
   os << ">==============================" << endl ;
   os << endl << endl ;
}


//---------------------------------------------------------------------------
void
PEL_Error:: notify( PEL_Module* list,
                    std::string const& message,
                    std::string const& where,
                    stringVector const& choices )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Error:: notify list" ) ;
   PEL_CHECK_PRE( list!=0 ) ;
   PEL_CHECK_PRE( !message.empty() ) ;
   PEL_CHECK_PRE( !where.empty() ) ;

   PEL_Module* result = notify( list, message, where ) ;
   if(choices.size()>0)
      result->add_entry( "valid_choices",
                         PEL_StringVector::create( result, choices ) ) ;

}


//---------------------------------------------------------------------------
PEL_Module*
PEL_Error:: notify( PEL_Module* list,
                    std::string const& message,
                    std::string const& where )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Error:: notify simple" ) ;
   PEL_CHECK_PRE( list!=0 ) ;
   PEL_CHECK_PRE( !message.empty() ) ;
   PEL_CHECK_PRE( !where.empty() ) ;

   size_t nb = 0 ;
   if( list->has_entry( "nb_items" ) )
      nb = list->data_of_entry( "nb_items" )->to_int() ;
   else
      list->add_entry( "nb_items", PEL_Int::create( list, nb ) ) ;

   std::ostringstream is ;
   is << "Error" << nb ;
   PEL_Module* result = PEL_Module::create( list, is.str() ) ;
   result->add_entry( "message", PEL_String::create( result, message ) ) ;
   result->add_entry( "where", PEL_String::create( result, where ) ) ;

   nb++ ;

   const_cast<PEL_Int*>(
      static_cast<PEL_Int const*>(list->data_of_entry( "nb_items" )))
      ->set( (int)nb ) ;
   list->add_module( result ) ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->has_entry( "message" ) ) ;
   PEL_CHECK_POST( result->has_entry( "where" ) ) ;

   return result ;
}


//---------------------------------------------------------------------------
void
PEL_Error:: display_data_checking( PEL_ModuleExplorer* report ) const
//---------------------------------------------------------------------------
{
   if( !report->is_empty() )
   {
      os << "Data deck checking failed : " << endl ;
      os << "*************************" << endl ;
      os << endl ;

      for( report->start_module_iterator() ;
           report->is_valid_module() ;
           report->go_next_module() )
      {
         PEL_ModuleExplorer const* error =
            report->create_subexplorer( 0 ) ;
         os << "   " << error->string_data( "message" ) << endl ;
         os << "   in module " << error->string_data( "where" ) << endl ;
         if( error->has_entry( "valid_choices" ) )
         {
            stringVector val = error->stringVector_data( "valid_choices" ) ;
            val.sort() ;
            os << "   Valid ones are :" << endl ;
            for( size_t j=0 ; j<val.size() ; ++j )
            {
               os << "      - \"" << val(j) << "\"" << endl ;
            }
         }
         os << endl ;
         error->destroy() ; error = 0 ;
      }
   }
}

//---------------------------------------------------------------------------
void
PEL_Error:: exit_with_internal_error( void )
//---------------------------------------------------------------------------
{
   PEL_Exec::set_exit_code( -1 ) ;
   throw PEL_Exceptions::InternalError() ;
}

//---------------------------------------------------------------------------
void
PEL_Error:: exit( size_t status )
//---------------------------------------------------------------------------
{
   PEL_Exec::set_exit_code( status ) ;
   throw PEL_Exceptions::UserError() ;
}



