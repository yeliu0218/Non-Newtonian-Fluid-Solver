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

#include <PEL_Comparator.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Data.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_FileToModule.hh>
#include <PEL_Int.hh>
#include <PEL_KeywordDataPair.hh>
#include <PEL_KeywordDataIterator.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ModuleIterator.hh>
#include <PEL_StringVector.hh>
#include <PEL_System.hh>
#include <PEL_TICio.hh>
#include <PEL_String.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>
#include <doubleVector.hh>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

using std::cout ;
using std::endl ;
using std::string ;

PEL_Comparator const* PEL_Comparator::PROTOTYPE = new PEL_Comparator() ;

struct PEL_Comparator_ERROR
{
   static void n0( std::string const& filename ) ;
   static void n1( std::string const& file1, std::string const& format1,
                   std::string const& file2, std::string const& format2 ) ;
   static void n2( PEL_ModuleExplorer const* exp, 
                   std::string const& choices ) ;
} ;

//----------------------------------------------------------------------------
PEL_Comparator:: PEL_Comparator( void )
//----------------------------------------------------------------------------
   : PEL_Application( "pelcmp" )
   , FILE1()
   , FILE2()
   , FILE_OUT()
   , EXP( 0 )
   , IGNORE( 0 )
{
}

//----------------------------------------------------------------------------
PEL_Comparator*
PEL_Comparator:: create_replica( PEL_Object* a_owner,
                                 PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Comparator:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;
   
   PEL_Comparator* result = new PEL_Comparator( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_Comparator:: PEL_Comparator( PEL_Object* a_owner, 
                                 PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , FILE1( exp->string_data( "left_file" ) )
   , FILE2( exp->string_data( "right_file" ) )
   , FILE_OUT()
   , EXP( 0 )
   , VERBOSE( false )
   , MY_DBL_EPS( PEL::bad_double() )
   , MY_DBL_MIN( PEL::bad_double() )
   , IGNORE( 0 )
   , FORMAT( "UNKNOWN" )
{
   if( exp->has_entry( "output_file" ) )
      FILE_OUT = exp->string_data( "output_file" ) ;
   if( exp->has_entry( "verbose" ) )
      VERBOSE = exp->bool_data( "verbose" ) ;
   if( exp->has_module( "PEL_ModuleComparator" ) )
      EXP = exp->create_subexplorer( this, "PEL_ModuleComparator" ) ;
   if( exp->has_module( "tolerance" ) )
   {
      PEL_ModuleExplorer const* ee = 
                         exp->create_subexplorer( 0 ,"tolerance" ) ;
      MY_DBL_EPS = ee->double_data( "dbl_eps" ) ;
      MY_DBL_MIN = ee->double_data( "dbl_min" ) ;
      ee->destroy() ;
   }
   if( exp->has_entry( "ignore_data" ) )
      IGNORE = exp->stringVector_data( "ignore_data" ) ;
   if( exp->has_entry( "format" ) )
   {
      FORMAT = exp->string_data( "format" ) ;
   }
}

//----------------------------------------------------------------------------
PEL_Comparator*
PEL_Comparator:: create_replica_from_args( PEL_Object* a_owner,
                                           stringVector& args ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Comparator:: create_replica_from_args" ) ;
   
   PEL_Comparator* result = new PEL_Comparator( a_owner, args ) ;

   PEL_CHECK( create_replica_from_args_POST( result, a_owner, args ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_Comparator:: PEL_Comparator( PEL_Object* a_owner, stringVector& args )
//----------------------------------------------------------------------------
   : PEL_Application( a_owner, 0 )
   , FILE1()
   , FILE2()
   , FILE_OUT()
   , EXP( 0 )
   , VERBOSE( false )
   , MY_DBL_EPS( PEL::bad_double() )
   , MY_DBL_MIN( PEL::bad_double() )
   , IGNORE( 0 )
   , FORMAT( "UNKNOWN" )
{
   size_t i_non_opt = 0 ;
   while( args.size() != 0 )
   {
      string arg = args( 0 ) ;
      args.remove_at( 0 ) ;

      // no options
      if( arg[0] == '-' )
      {         
         if( arg == "-verbose" )
         {
            VERBOSE = true ;
         }
         else if( arg=="-dbl_eps" && args.size()>0 )
         {
            std::istringstream is( args(0) ) ;
            is >> MY_DBL_EPS ;
            args.remove_at( 0 ) ;
         }
         else if( arg=="-dbl_min" && args.size()>0 )
         {
            std::istringstream is( args(0) ) ;
            is >> MY_DBL_MIN ;
            args.remove_at( 0 ) ;
         }
         else if( arg=="-ignore_data" && args.size()>0 )
         {
            IGNORE.append( args( 0 ) ) ;
            args.remove_at( 0 ) ;
         }
         else if( arg=="-format" && args.size()>0 )
         {
            FORMAT = args( 0 ) ;
            args.remove_at( 0 ) ;
         }
         else if( arg=="-output_file" && args.size()>0 )
         {
            FILE_OUT = args( 0 ) ;
            args.remove_at( 0 ) ;
         }
         else
         {         
            notify_error_in_arguments() ;
         }
      }
      else
      {
         i_non_opt++ ;
         if( i_non_opt==1 ) 
            FILE1 = arg ;
         else if( i_non_opt==2 )
            FILE2 = arg ;
         else if( i_non_opt==3 )
            FILE_OUT = arg ;
         else
            notify_error_in_arguments() ;
      }
   }
   if( ! ( ( MY_DBL_EPS==PEL::bad_double() && MY_DBL_MIN==PEL::bad_double() ) ||
           ( MY_DBL_EPS!=PEL::bad_double() && MY_DBL_MIN!=PEL::bad_double() ) ) )
       notify_error_in_arguments() ;
   if( FILE1.empty() || FILE2.empty() ) notify_error_in_arguments() ;
}

//----------------------------------------------------------------------------
PEL_Comparator:: ~PEL_Comparator( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
PEL_Comparator:: set_preferred_motifs_formats( 
                                          PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Comparator:: set_preferred_motifs_formats" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   
   pref_motifs()  = exp->stringVector_data( "motifs" ) ;
   pref_formats() = exp->stringVector_data( "formats" ) ;
   exp->test_data( "formats", "size(formats)=size(motifs)" ) ;
   
   for( size_t i=0 ; i<pref_formats().size() ; i++ )
   {
      std::string const& frmt = pref_formats()( i ) ;
      if( ! PEL_FileToModule::has( frmt ) ) 
         PEL_Comparator_ERROR::n2( exp, PEL_FileToModule::list_of_formats() ) ;
   }
}

//----------------------------------------------------------------------------
void
PEL_Comparator:: detect_file_format( PEL_ModuleExplorer const* exp,
                                     std::string const& file_name,
                                     std::string& format )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Comparator:: detect_file_format" ) ;

   if( exp!=0 && exp->has_entry( "format" ) )
   {
      format = exp->string_data( "format" ) ;
      exp->test_data_in( "format", PEL_FileToModule::list_of_formats() ) ;
   }
   else
   {
      bool found = false ;
      for( size_t n=0 ; n<pref_motifs().size() ; n++ )
      {
         if( file_name.find( pref_motifs()(n) ) < file_name.length() )
         {
            found = true ;
            format = pref_formats()(n) ;
            break ;
         }
      }
      if( !found ) PEL_FileToModule::find_file_format( file_name, format ) ;
   }
}

//----------------------------------------------------------------------------
void
PEL_Comparator:: run( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Comparator::run" ) ;

   std::ofstream ofs ;
   if( !FILE_OUT.empty() )
   {
      ofs.open( FILE_OUT.c_str(), std::ios::out|std::ios::trunc ) ;
      if( ! ofs )
      {
         string mess = "PEL_Comparator: unable to open file \"" ;
         mess += FILE_OUT ;
         mess += "\"" ;
         PEL_Error::object()->raise_plain( mess ) ;
      }
   }
   
   if( FORMAT == "UNKNOWN" )
   {
      std::string fmt = FORMAT ;
      detect_file_format( 0, FILE1, fmt ) ;
      detect_file_format( 0, FILE2, FORMAT ) ;
      if( FORMAT != fmt ) PEL_Comparator_ERROR::n1( FILE1, fmt, 
                                                    FILE2, FORMAT ) ;
   }

   if( VERBOSE ) PEL::out() << "|    format: " << FORMAT << endl ;
   if( FORMAT != "UNKNOWN" )
   {
      if( FILE_OUT.empty() )
         builtin_compare( PEL::out() ) ;
      else
         builtin_compare( ofs ) ;
   }
   else
   {
      if( FILE_OUT.empty() )
         do_diff( PEL::out() ) ;
      else
         builtin_compare( ofs ) ;
   }
   
   if( !FILE_OUT.empty() ) ofs.close() ;
}

//-----------------------------------------------------------------------------
void
PEL_Comparator:: builtin_compare( std::ostream& os )
//-----------------------------------------------------------------------------
{  
   PEL_LABEL( "PEL_Comparator:: builtin_compare" ) ;
   
   PEL_FileToModule const* oo = PEL_FileToModule::object( FORMAT ) ;
   PEL_Module* mod1 = oo->create_from_file( 0, "left_file",  FILE1 ) ;
   PEL_Module* mod2 = oo->create_from_file( 0, "right_file", FILE2 ) ;

   if( EXP==0 ) 
   {
      PEL_Module* comp = PEL_Module::create( this, "PEL_ModuleComparator" ) ;
      
      if( IGNORE.size()>0  )
      {
         comp->add_entry( "ignore_data", 
                          PEL_StringVector::create( comp, IGNORE ) ) ;
      }
      if( MY_DBL_EPS != PEL::bad_double()  )
      {
         comp->add_entry( "dbl_eps", PEL_Double::create( comp, MY_DBL_EPS ) ) ;
      }
      if( MY_DBL_MIN != PEL::bad_double()  )
      {
         comp->add_entry( "dbl_min", PEL_Double::create( comp, MY_DBL_MIN ) ) ;
      }
      if( comp!=0 ) EXP = PEL_ModuleExplorer::create( this, comp ) ;     
   }
   
   if( EXP!=0 && VERBOSE )
   {
      PEL::out() << "PEL_ModuleComparator configuration : " << std::endl ;
      EXP->print( PEL::out(), 2 ) ;
      PEL::out() << std::endl ;
   }
   
   PEL_Module* mdif = PEL_Module::create_as_difference( 
                                            this, "result", mod1, mod2, EXP );

   mdif->add_entry( "left_file",  PEL_String::create( mdif, FILE1 ) ) ;
   mdif->add_entry( "right_file", PEL_String::create( mdif, FILE2 ) ) ;

   int nb_errors = mdif->data_of_entry( "nb_differences" )->to_int() ;   
   if( nb_errors>0 )
   {
      if( FILE_OUT.empty() ) // for consistance with the old versions
      {
         os << "---> " << nb_errors << " difference" ;
         if( nb_errors > 1 ) os << "s" ;
         os << endl;
      }
      display_glob_info_diff( os, mdif ) ;
      display_submodule_diff( os, mdif ) ;
      PEL_Exec::set_exit_code( nb_errors ) ;
   }
   
   mod1->destroy() ; mod1 = 0 ;
   mod2->destroy() ; mod2 = 0 ;
}

//----------------------------------------------------------------------------
void
PEL_Comparator:: do_diff( std::ostream& os ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Comparator:: do_diff" ) ;
   if( VERBOSE ) PEL::out() << "|    do_diff " << FILE1 
                       << " and " << FILE2 << std::endl ;
   
   std::ifstream ifs1( FILE1.c_str() ) ;
   std::ifstream ifs2( FILE2.c_str() ) ;

   bool ok = ifs1.good() && ifs2.good() ;
   std::string fline ;
   std::string line1, line2 ;
   size_t iline = 0 ;
   
   while( ifs1.good() && ifs2.good() )
   {
      iline++ ;
      getline( ifs1, line1 ) ;
      getline( ifs2, line2 ) ;
      if( line1 != line2 ) 
      {
         // For linefeed in files written under windows
         PEL::replace( line1, "\r", "" ) ;
         PEL::replace( line2, "\r", "" ) ;
         
         // because windows notation is not the same than linux one...
         PEL::replace( line1, "E+00", "E+0" ) ;
         PEL::replace( line2, "E+00", "E+0" ) ;
         PEL::replace( line1, "E-00", "E-0" ) ;
         PEL::replace( line2, "E-00", "E-0" ) ;
         PEL::replace( line1, "e+00", "e+0" ) ;
         PEL::replace( line2, "e+00", "e+0" ) ;
         PEL::replace( line1, "e-00", "e-0" ) ;
         PEL::replace( line2, "e-00", "e-0" ) ;
         PEL::replace( line1, "  ", " " ) ;
         PEL::replace( line2, "  ", " " ) ;

         if( line1 != line2 )
         {
            os << "   first difference at line " << iline << endl ;
            if( line1.length() < 100 && line2.length() < 100 )
            {
               os << "   \"" << line1 << "\"" << endl ;
               os << "   \"" << line2 << "\"" << endl ;
            }
            PEL_Exec::set_exit_code( 12 ) ;
            ok = false ;
            break ;
         }
         
      }
   }
   if( ok && ( ifs1.good() || ifs2.good() ) ) 
   {
      PEL_Exec::set_exit_code( 12 ) ;
   } 
}

//---------------------------------------------------------------------------
void
PEL_Comparator:: print_usage( void ) const
//---------------------------------------------------------------------------
{
   PEL::out() << usage_title( "pelcmp" )  ;
   PEL::out() << "[options...] file1 file2 [file3]" << endl << endl ;
   PEL::out() << "     Compare and display differences between pairs of files"
              << endl ;
}

//---------------------------------------------------------------------------
void
PEL_Comparator:: print_operands( void ) const
//---------------------------------------------------------------------------
{
   PEL::out() << "     -dbl_eps <value>" << endl
	<< "          argument a_dbl_eps in calls to PEL::double_equality"
        << endl 
        << "          when comparing values of type double" 
        << endl << endl ;
   PEL::out() << "     -dbl_min <value>" << endl
	<< "          argument a_dbl_min in calls to PEL::double_equality"
        << endl 
        << "          when comparing values of type double" 
        << endl << endl ;
   PEL::out() << "     -ignore_data <keyword>" << endl
	<< "          Ignore entries of keyword <keyword>" 
        << endl << endl ;
   PEL::out() << "     -format [" << PEL_FileToModule::list_of_formats() 
                                  << "]" << endl
	<< "          Specify the format of file1 and file2" 
        << endl << endl ;
   PEL::out() << operands_title() ;
   PEL::out() << "     file1, file2" << endl
	<< "          Path names of the two files to be compared" 
        << endl << endl ;
   PEL::out() << "     file3" << endl
        << "          Path names of the file containing the difference on exit"
	<< endl << endl ;
}

//---------------------------------------------------------------------------
void
PEL_Comparator:: print_exit_status( void ) const
//---------------------------------------------------------------------------
{
   PEL::out() << exit_status_title() ;
   PEL::out() << "     0    The files are identical" << endl ;
   PEL::out() << "    >0    Number of differences found" << endl ;
}

//----------------------------------------------------------------------------
void
PEL_Comparator:: display_glob_info_diff( std::ostream& os,
                                         PEL_Module const* mod ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Comparator:: display_glob_info_diff" ) ;

   PEL_KeywordDataIterator* it = mod->create_entry_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      string const& keywd  = it->item()->keyword() ;
      PEL_Data const* data = it->item()->data() ;
      if( keywd!="left_file" && keywd!="right_file" &&
          keywd!="nb_differences" )
      {
         PEL_ASSERT( data->data_type() == PEL_Data::String ) ;
         os << keywd << ": " << data->to_string() << endl ;
      }
   }
   it->destroy() ; it = 0 ;
}

//----------------------------------------------------------------------------
void
PEL_Comparator:: display_submodule_diff( std::ostream& os,
                                         PEL_Module const* mod ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Comparator:: display_submodule_diff" ) ;

   PEL_ModuleIterator* it = mod->create_module_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PEL_Module const* mm = it->item() ;
      string modname =  mm->absolute_path_name() ;
      modname.erase( modname.begin(), modname.begin()+8 ) ; 
      PEL_KeywordDataIterator* kwit = mm->create_entry_iterator( 0 ) ;
      bool first_diff = true ;
      for( kwit->start() ; kwit->is_valid() ; kwit->go_next() )
      {
         string const& str = kwit->item()->keyword() ;
         PEL_Data const* data = kwit->item()->data() ;
         if( first_diff && data->data_type()!=PEL_Data::DoubleVector )
         {
            os << modname << endl ;
            first_diff=false ;
         }
         std::ios_base::fmtflags original_flags = os.flags() ;
         os.setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
         std::streamsize p = os.precision() ;
         os << std::setprecision( 5 ) ;
         switch( data->data_type() )
         {
            case PEL_Data::DoubleArray2D :
            {
               // the compared data is of type PEL_Data::DoubleVector
               double max = -PEL::max_double() ;
               doubleArray2D const& stat = data->to_double_array2D() ;
               size_t i_max = PEL::bad_index() ;
               size_t i_other = PEL::bad_index() ;
               for( size_t i=0 ; i<stat.index_bound( 0 ) ; ++i )
               {
                  if( ( stat( i, 0 ) < 0.0 ) && ( i_other == PEL::bad_index() ) )
                  {
                     i_other = i ;
                  }
                  if( stat( i,0 ) > max )
                  {
                     i_max = i ;
                     max = stat( i, 0 ) ;
                  }
               }
               if( i_max != PEL::bad_index() &&  max > 0.0 )
               {  
                  os << "    relative error max : " << max << endl ;
                  os << "    obtained for -> " << stat( i_max, 1 ) << " <!!> ";
                  os << stat( i_max, 2 ) << " <-" << endl ;
               }
               else
               {
                  os << "    first difference at index " << i_other << endl ;
                  os << "    -> " << stat( i_other, 1 ) << " <!!> ";
                  os << stat( i_other, 2 ) << " <-" << endl ;
               }
            }
            break ;
            case PEL_Data::DoubleArray3D :
            {
               // the compared data is of type PEL_Data::DoubleArray2D
               double max = -PEL::max_double() ;
               doubleArray3D const& stat = data->to_double_array3D() ;
               size_t i_max = PEL::bad_index() ;
               size_t j_max = PEL::bad_index() ;
               size_t i_other = PEL::bad_index() ;
               size_t j_other = PEL::bad_index() ;
               for( size_t i=0 ; i<stat.index_bound( 0 ) ; ++i )
               {
                  for( size_t j=0 ; j<stat.index_bound( 1 ) ; ++j )
                  {
                     if( ( stat( i, j, 0 ) < 0.0 ) && 
                         ( i_other == PEL::bad_index() ) )
                     {
                        i_other = i ;
                        j_other = j ;
                     }
                     if( stat( i, j , 0 ) > max ) 
                     {
                        i_max = i ;
                        j_max = j ;
                        max = stat( i, j, 0 ) ;
                     }
                  }
               }
               if( i_max != PEL::bad_index() &&  max > 0.0 )
               {  
                  os << "    relative error max : " << max << endl ;
                  os << "    obtained for -> " << stat( i_max, j_max, 1 ) 
                     << " <!!> " << stat( i_max, j_max, 2 ) << " <-" << endl ;
               }
               else
               {
                  os << "    first difference at indices (" 
                     << i_other << "," << j_other << ")" << endl ;
                  os << "    -> " << stat( i_other, j_other, 1 ) << " <!!> ";
                  os << stat( i_other, j_other, 2 ) << " <-" << endl ;
               }
            }
            break ;
            case PEL_Data::String :
            {
               os << " " << str << " : " << data->to_string() << endl ;
            }
            default:
            {
            }
         }
         os << std::setprecision(p) ;
         os.flags( original_flags ) ;
      }
      kwit->destroy() ;
      display_submodule_diff( os, mm ) ; 
   }
   it->destroy() ;
}

//----------------------------------------------------------------------
stringVector&
PEL_Comparator:: pref_motifs( void )
//----------------------------------------------------------------------
{
   static stringVector result( 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
stringVector&
PEL_Comparator:: pref_formats( void )
//----------------------------------------------------------------------
{
   static stringVector result( 0 ) ;
   return( result ) ;
}

//internal--------------------------------------------------------------
void 
PEL_Comparator_ERROR:: n0( std::string const& filename )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** PEL_Comparator error:" << endl << endl ;
   mesg << "    unable to open file \"" <<  filename << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PEL_Comparator_ERROR:: n1( std::string const& file1,
                           std::string const& format1,
                           std::string const& file2,
                           std::string const& format2 )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** PEL_Comparator error:" << endl << endl ;
   mesg << "the two files should have the same format" << endl << endl ;
   mesg << file1 << endl << "   detected format: " << format1 << endl << endl ;
   mesg << file2 << endl << "   detected format: " << format2 << endl ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PEL_Comparator_ERROR:: n2( PEL_ModuleExplorer const* exp,
                           std::string const& choices )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "In module: " << endl << "   " << exp->absolute_path_name() << endl ;
   mesg << "the data of keyword \"formats\"" << endl ;
   mesg << "should contains strings with one of the following values: " ;
   stringVector liste( 0 ) ;
   liste.build_from_string( choices, ',' ) ;
   for( size_t i=0 ; i<liste.size() ; i++ )
   {
      mesg << endl << "   - \"" << liste(i) << "\"" ;
   }
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}



