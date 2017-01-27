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

#include <PEL_NativeFileToModule.hh>

#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Int.hh>
#include <PEL_Module.hh>
#include <PEL_StringVector.hh>
#include <PEL_TICio.hh>
#include <PEL_assertions.hh>
#include <stringVector.hh>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>

using std::endl ;
using std::string ;

PEL_NativeFileToModule* 
PEL_NativeFileToModule:: PEL_OBJ = new PEL_NativeFileToModule( "PEL", 
                                                               ".pel" ) ;

PEL_NativeFileToModule* 
PEL_NativeFileToModule:: GENE_OBJ = new PEL_NativeFileToModule( "GENE", 
                                                                ".gene" ) ;

PEL_NativeFileToModule* 
PEL_NativeFileToModule:: CSV_OBJ = new PEL_NativeFileToModule( "CSV", 
                                                               ".csv" ) ;

struct PEL_NativeFileToModule_ERROR
{
   static void n0( std::string const& file_name ) ;
   static void n1( std::string const& file_name, size_t nb_cols,
                   std::string const& value, std::string const& separator ) ;
} ;

//----------------------------------------------------------------------
PEL_NativeFileToModule:: PEL_NativeFileToModule( std::string const& a_format,
                                           std::string const& a_default_motif )
//----------------------------------------------------------------------
   : PEL_FileToModule( a_format, a_default_motif )
{
}

//----------------------------------------------------------------------
PEL_NativeFileToModule:: ~PEL_NativeFileToModule( void )
//----------------------------------------------------------------------
{
   PEL_OBJ = 0 ;
   GENE_OBJ = 0 ; 
   CSV_OBJ = 0 ;
}

//----------------------------------------------------------------------
PEL_Module*
PEL_NativeFileToModule:: create_from_file( 
                                 PEL_Object* a_owner,
                                 std::string const& module_name,
                                 std::string const& file_name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_NativeFileToModule:: create_from_file" ) ;
   PEL_CHECK_PRE( create_from_file_PRE( a_owner, module_name, file_name ) ) ;
   
   PEL_Module* result = 0 ;
   if( format() == "PEL" )
   {
      result = PEL_Module::create( a_owner, module_name, file_name, 
                                   PEL_Exec::execution_context() ) ;
   }
   else if( format() == "GENE" )
   {
      result = PEL_TICio::create_from_gene_file( a_owner, 
                                                 module_name, 
                                                 file_name ) ;      
   }
   else if( format() == "CSV" )
   {
      result = create_from_multicolumns( a_owner, module_name, file_name,
                                         "," ) ;
   }
   
   PEL_CHECK_POST( create_from_file_POST( result, 
                                          a_owner, module_name, file_name ) ) ;
   return( result) ;
}

//----------------------------------------------------------------------------
PEL_Module*
PEL_NativeFileToModule:: create_from_multicolumns( 
                                       PEL_Object* a_owner,
                                       std::string const& name,
                                       std::string const& file_name,
                                       std::string const& separator ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_NativeFileToModule:: create_from_multicolumns" ) ;

   typedef string::const_iterator iter ;

   std::ifstream in( file_name.c_str() ) ;
   if( !in ) PEL_Error::object()->raise_plain( "unable to open file \"" + 
                                               file_name + "\"" ) ;

   PEL_Module* result = PEL_Module::create( a_owner, name ) ;
   
   std::vector< doubleVector > double_vals ;
   std::vector< stringVector > string_vals ;

   bool first_line = true ;
   size_t nb_columns = 0 ;

   string line ;
   while( getline( in, line ) )
   {
      stringVector chaines( 0 ) ;

      size_t i_col = 0 ;

      iter it_b = line.begin() ;
      iter it_e = line.end() ;
      iter it_i = it_b ;
      while( it_b != it_e )
      {
         it_i = std::search( it_i, it_e, separator.begin(), separator.end() ) ;
         chaines.append( string( it_b, it_i ) ) ;

         if( it_i != it_e ) it_i += separator.size() ;
         it_b = it_i ;

         i_col++ ;
      }

      if( first_line )
      {
         nb_columns = i_col ;
         double_vals.resize( nb_columns, doubleVector( 0 ) ) ;
         string_vals.resize( nb_columns, stringVector( 0 ) ) ;
         first_line = false ;
      }
      else
      {
         if( i_col != nb_columns ) 
            PEL_NativeFileToModule_ERROR::n0( file_name ) ;
      }
      for( size_t j=0 ; j<nb_columns ; ++j )
      {
         std::istringstream is( chaines( j ) ) ;
         double val ;
         is >> val ;
         if( is )
         {
            std::string str ;
            is >> str ;
            if( str.empty() ) 
            {  
               double_vals[ j ].append( val ) ;
               string_vals[ j ].append( ""  ) ;
            }
            else
            {
               PEL_NativeFileToModule_ERROR::n1( file_name, nb_columns,
                                                 chaines( j ), separator ) ;
            }
            
         }
         else
         {
            double_vals[ j ].append( 0.0 ) ;
            string_vals[ j ].append( chaines( j ) ) ;
         }
      }
   }

   result->add_entry( "nb_columns", PEL_Int::create( result, nb_columns ) ) ;
   for( size_t j=0 ; j<nb_columns ; j++ )
   {
      std::ostringstream nn ;
      nn << "column_" << j ;
      PEL_Module* mod = PEL_Module::create( result, nn.str() ) ;
      mod->add_entry( "double_values", 
                      PEL_DoubleVector::create( mod, double_vals[ j ] ) ) ;
      mod->add_entry( "string_values", 
                      PEL_StringVector::create( mod, string_vals[ j ] ) ) ;
      result->add_module( mod ) ;
   }   
   
   return result ;
}

//internal--------------------------------------------------------------
void 
PEL_NativeFileToModule_ERROR:: n0( std::string const& file_name )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** PEL_Comparator error:" << endl << endl ;
   mesg << "    error in CSV file \"" <<  file_name << "\"" << endl ;
   mesg << "    each line should contain the same number of fields" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}


//internal--------------------------------------------------------------
void 
PEL_NativeFileToModule_ERROR:: n1( std::string const& file_name,
                                   size_t nbcols,
                                   std::string const& value,
                                   std::string const& separator )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** PEL_Comparator error:" << endl << endl ;
   mesg << "    error in CSV file \"" <<  file_name << "\"" << endl ;
   mesg << "    (number of columns found: " << nbcols << " )" << endl ;
   mesg << "    unable to convert the value" << endl ;
   mesg << "       \"" << value << "\"" << endl ;
   mesg << "    into double or string" << endl << endl ;
   mesg << "    hint: the separator may be incorrect" << endl ;
   mesg << "          (it is supposed to be \"" << separator << "\")" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

