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

#include <PEL_BinStored_TEST.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_BinStored.hh>
#include <PEL_Error.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Double.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>
#include <PEL_System.hh>

#include <intVector.hh>
#include <intArray2D.hh>
#include <intArray3D.hh>
#include <doubleVector.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>
#include <stringVector.hh>

#include <fstream>

PEL_BinStored_TEST*
PEL_BinStored_TEST::registered_test = new PEL_BinStored_TEST() ;

//-------------------------------------------------------------------------
void
PEL_BinStored_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BinStored_TEST:: process_all_tests" ) ;
   
   std::string const file = "test.bin" ;
   PEL_ModuleExplorer* datum_exp =
      data_deck_explorer()->create_subexplorer( 0, "DATA" ) ;
   size_t record = 0 ;
   PEL_BinStored::init_binary_file( file ) ;
   PEL_Module* mod = PEL_Module::create( 0, "BINARY" ) ;
   stringVector keywords(0) ;
   PEL_Vector * vector = PEL_Vector::create( 0, 0 ) ;
   
   // First saving binary data
   for( datum_exp->start_module_iterator() ;
        datum_exp->is_valid_module() ;
        datum_exp->go_next_module() )
   {
      PEL_ModuleExplorer* under = datum_exp->create_subexplorer( datum_exp ) ;
      std::string const& t_name = under->string_data( "type" ) ;
      std::string const& a_name = under->string_data( "key" ) ;
      keywords.append( a_name ) ;
      PEL_Data const* entry = under->abstract_data( vector, a_name ) ;
      if( !entry->value_can_be_evaluated(0) )
      {
         PEL_Error::object()->raise_not_evaluable(
               under, a_name, entry->undefined_variables(0) ) ;
      }
      notify_one_test_result( a_name+" type",
                           PEL_Data::type_name(entry->data_type())==t_name ) ;
      vector->append( const_cast<PEL_Data*>( entry ) ) ;
      PEL_BinStored const* expression =
         PEL_BinStored::create_reference( mod, entry, file ) ;
      // Verifying binary data saved as internal expression generated
      bool ok = equal_data( entry, expression ) ;
      notify_one_test_result( a_name+" internal", ok ) ;
      mod->add_entry( a_name, expression ) ;
   }

   // Verify binary data saved as external way
   PEL_ModuleExplorer* binary_exp = 
      data_deck_explorer()->create_subexplorer( 0, "BINARY_EXP" ) ;
   record = 0 ;
   for( size_t i=0 ; i<keywords.size() ; i++ )
   {
      PEL_Data const* recovered =
                           binary_exp->abstract_data( vector, keywords(i) ) ;
      if( !recovered->value_can_be_evaluated(0) )
      {
         PEL_Error::object()->raise_not_evaluable(
               binary_exp, keywords(i), recovered->undefined_variables(0) ) ;
      }
      PEL_Data const* entry = static_cast<PEL_Data const*>( vector->at(i) ) ;
      bool ok = equal_data( entry, recovered ) ;
      notify_one_test_result( keywords(i)+" external", ok ) ;
   }
   // Append data to another file
   std::string const another_file = "anothertest.bin" ;
   PEL_BinStored::init_binary_file( another_file ) ;

   // First saving binary data
   for( size_t i=0 ; i<keywords.size() ; i++ )
   {
      std::string const& a_name = keywords(i) ;
      PEL_Data const* entry = static_cast<PEL_Data const*>( vector->at(i) ) ;
      PEL_BinStored const* expression =
         PEL_BinStored::create_reference( vector, entry, another_file ) ;
      // Verifying binary data saved as internal expression generated
      bool ok = equal_data( entry, expression ) ;
      notify_one_test_result( a_name+" another internal", ok ) ;
   }

   // Testing expression generating
   std::ofstream ostr( "binary.pel" , std::ios::out ) ;
   mod->print( ostr, 0 ) ;
   ostr.close() ;

   PEL_Module* recov = PEL_Module::create( 0, "ROOT", "binary.pel" ) ;
   PEL_ModuleExplorer* root_exp = PEL_ModuleExplorer::create( recov, recov ) ;
   PEL_ModuleExplorer* exp = root_exp->create_subexplorer( recov, "BINARY" ) ;
   size_t j = 0 ;
   for( exp->start_entry_iterator() ;
        exp->is_valid_entry() ;
        exp->go_next_entry() )
   {
      std::string const& a_name = exp->keyword() ;
      PEL_Data const* entry = exp->abstract_data( recov, a_name ) ;
      if( !entry->value_can_be_evaluated(0) )
      {
         PEL_Error::object()->raise_not_evaluable(
               exp, a_name, entry->undefined_variables(0) ) ;
      }
      PEL_Data const* expression =
         static_cast<PEL_Data const*>( vector->at(j) ) ;
      // Verifying binary data saved as internal expression read
      bool ok = equal_data( entry, expression ) ;
      notify_one_test_result( a_name+" throught file", ok ) ;
      j++;
   }

   PEL_System::erase( "test.bin" ) ;
   PEL_System::erase( "anothertest.bin" ) ;
   PEL_System::erase( "binary.pel" ) ;

   datum_exp->destroy() ; datum_exp = 0 ;
   mod->destroy() ; mod = 0 ;
   vector->destroy() ; vector = 0 ;
   binary_exp->destroy() ; binary_exp = 0 ;
   recov->destroy() ; recov = 0 ;
}

//-------------------------------------------------------------------------
bool
PEL_BinStored_TEST:: equal_data( PEL_Data const* entry,
                                 PEL_Data const* recovered ) 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BinStored_TEST:: equal_data" ) ;
   PEL_CHECK( entry != 0 ) ;
   PEL_CHECK( recovered != 0 ) ;
   
   bool result = entry->data_type()==recovered->data_type() ;
   if( result )
   {
      if( entry->data_type()==PEL_Data::Double )
      {
         result = entry->to_double()==recovered->to_double() ;
      }
      else if( entry->data_type()==PEL_Data::Int )
      {
         result = entry->to_int()==recovered->to_int() ;
      }
      else if( entry->data_type()==PEL_Data::Bool )
      {
         result = entry->to_bool()==recovered->to_bool() ;
      }
      else if( entry->data_type()==PEL_Data::String )
      {
         result = entry->to_string()==recovered->to_string() ;
      }
      else if( entry->data_type()==PEL_Data::DoubleVector )
      {
         doubleVector const& vec1 = entry->to_double_vector() ;
         doubleVector const& vec2 = recovered->to_double_vector() ;
         size_t const nb = vec1.size() ;
         result = ( vec1.size() == vec2.size() ) ;
         for( size_t i=0 ; result && i<nb ; i++ )
         {
            result &= vec1(i)==vec2(i) ;
         }
      }
      else if( entry->data_type()==PEL_Data::IntVector )
      {
         intVector const& vec1 = entry->to_int_vector() ;
         intVector const& vec2 = recovered->to_int_vector() ;
         size_t const nb = vec1.size() ;
         result = ( vec1.size() == vec2.size() ) ;
         for( size_t i=0 ; result && i<nb ; i++ )
         {
            result &= vec1(i)==vec2(i) ;
         }
      }
      else if( entry->data_type()==PEL_Data::BoolVector )
      {
         boolVector const& vec1 = entry->to_bool_vector() ;
         boolVector const& vec2 = recovered->to_bool_vector() ;
         size_t const nb = vec1.size() ;
         result = ( vec1.size() == vec2.size() ) ;
         for( size_t i=0 ; result &&  i<nb ; i++ )
         {
            result &= vec1(i)==vec2(i) ;
         }
      }
      else if( entry->data_type()==PEL_Data::StringVector )
      {
         stringVector const& vec1 = entry->to_string_vector() ;
         stringVector const& vec2 = recovered->to_string_vector() ;
         size_t const nb = vec1.size() ;
         result = ( vec1.size() == vec2.size() ) ;
         for( size_t i=0 ; result &&  i<nb ; i++ )
         {
            result &= vec1(i)==vec2(i) ;
         }
      }
      else if( entry->data_type()==PEL_Data::DoubleArray2D )
      {
         doubleArray2D const& vec1 = entry->to_double_array2D() ;
         doubleArray2D const& vec2 = recovered->to_double_array2D() ;
         size_t const nb1 = vec1.index_bound(0) ;
         size_t const nb2 = vec1.index_bound(1) ;
         result = ( vec1.index_bound(0) == vec2.index_bound(0) &&
                    vec1.index_bound(1) == vec2.index_bound(1) ) ;
         for( size_t i=0 ; result && i<nb1 ; i++ )
         {
            for( size_t j=0 ; result && j<nb2 ; j++ )
            {
               result &= vec1(i,j)==vec2(i,j) ;
            }
         }
      }
      else if( entry->data_type()==PEL_Data::IntArray2D )
      {
         intArray2D const& vec1 = entry->to_int_array2D() ;
         intArray2D const& vec2 = recovered->to_int_array2D() ;
         size_t const nb1 = vec1.index_bound(0) ;
         size_t const nb2 = vec1.index_bound(1) ;
         result = ( vec1.index_bound(0) == vec2.index_bound(0) &&
                    vec1.index_bound(1) == vec2.index_bound(1) ) ;
         for( size_t i=0 ; result && i<nb1 ; i++ )
         {
            for( size_t j=0 ; result && j<nb2 ; j++ )
            {
               result &= vec1(i,j)==vec2(i,j) ;
            }
         }
      }
      else if( entry->data_type()==PEL_Data::DoubleArray3D )
      {
         doubleArray3D const& vec1 = entry->to_double_array3D() ;
         doubleArray3D const& vec2 = recovered->to_double_array3D() ;
         size_t const nb1 = vec1.index_bound(0) ;
         size_t const nb2 = vec1.index_bound(1) ;
         size_t const nb3 = vec1.index_bound(2) ;
         result = ( vec1.index_bound(0) == vec2.index_bound(0) &&
                    vec1.index_bound(1) == vec2.index_bound(1) &&
                    vec1.index_bound(2) == vec2.index_bound(2) ) ;
         for( size_t i=0 ; result && i<nb1 ; i++ )
         {
            for( size_t j=0 ; result && j<nb2 ; j++ )
            {
               for( size_t k=0 ; result && k<nb3 ; k++ )
               {
                  result &= vec1(i,j,k)==vec2(i,j,k) ;
               }
            }
         }
      }
      else if( entry->data_type()==PEL_Data::IntArray3D )
      {
         intArray3D const& vec1 = entry->to_int_array3D() ;
         intArray3D const& vec2 = recovered->to_int_array3D() ;
         size_t const nb1 = vec1.index_bound(0) ;
         size_t const nb2 = vec1.index_bound(1) ;
         size_t const nb3 = vec1.index_bound(2) ;
         result = ( vec1.index_bound(0) == vec2.index_bound(0) &&
                    vec1.index_bound(1) == vec2.index_bound(1) &&
                    vec1.index_bound(2) == vec2.index_bound(2) ) ;
         for( size_t i=0 ; result && i<nb1 ; i++ )
         {
            for( size_t j=0 ; result && j<nb2 ; j++ )
            {
               for( size_t k=0 ; result && k<nb3 ; k++ )
               {
                  result &= vec1(i,j,k)==vec2(i,j,k) ;
               }
            }
         }
      }
      else
      {
         PEL_Error::object()->raise_internal( "Type not implemented" ) ;
      }
   }
   return result ;
}

//-------------------------------------------------------------------------
PEL_BinStored_TEST:: PEL_BinStored_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "PEL_BinStored", "PEL_BinStored_TEST" )
{
}

//-------------------------------------------------------------------------
PEL_BinStored_TEST:: ~PEL_BinStored_TEST( void )
//-------------------------------------------------------------------------
{
}
