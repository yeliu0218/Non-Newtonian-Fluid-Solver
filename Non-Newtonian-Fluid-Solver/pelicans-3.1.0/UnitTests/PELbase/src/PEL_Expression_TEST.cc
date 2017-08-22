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

#include <PEL_Expression_TEST.hh>

#include <PEL.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Exceptions.hh>
#include <PEL_Error.hh>
#include <PEL_String.hh>
#include <PEL_Context.hh>
#include <PEL_Map.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>

#include <boolVector.hh>
#include <boolArray2D.hh>
#include <doubleVector.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>
#include <intVector.hh>
#include <intArray2D.hh>
#include <intArray3D.hh>
#include <stringVector.hh>
#include <stringArray2D.hh>


#include <iostream>

using std::string ;
using std::cerr ;
using std::endl ;

//-------------------------------------------------------------------------
PEL_Expression_TEST*
PEL_Expression_TEST::unique_instance = new PEL_Expression_TEST() ;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
PEL_Expression_TEST:: PEL_Expression_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "PEL_Expression", "PEL_Expression_TEST" ),
     modexp( 0 )
{
}



//-------------------------------------------------------------------------
PEL_Expression_TEST:: ~PEL_Expression_TEST( void )
//-------------------------------------------------------------------------
{
}


//-------------------------------------------------------------------------
void
PEL_Expression_TEST:: process_one_test( PEL_ModuleExplorer const* texp )
//-------------------------------------------------------------------------
{
   bool request_exception =
      texp->has_entry( "exception" ) && texp->bool_data( "exception" ) ;
   bool catched = false ;
   try
   {
      process_one_comparison( texp ) ;
   }
   catch( PEL_Exceptions::UserError )
   {
      catched = true ;
   }
   string nn = texp->name() ;
   if( request_exception || catched )
      notify_one_test_result(nn+" exception", catched==request_exception ) ;
}



//-------------------------------------------------------------------------
void
PEL_Expression_TEST:: process_one_comparison( PEL_ModuleExplorer const* texp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Expression_TEST:: process_one_comparison" ) ;
   
   string nn = texp->name() ;
   std::string const& eval = "exp_to_eval" ;
   PEL_Data const* to_eval = texp->abstract_data( this , eval ) ;
   PEL_DataWithContext const* res = texp->abstract_data( this , "result" ) ;
   if( !res->value_can_be_evaluated() )
   {
      PEL_Error::object()->raise_not_evaluable(
                             texp, "result", res->undefined_variables() ) ;
   }
   std::string const& t_name = texp->string_data( "type" ) ;
   bool ok = to_eval->data_type() == res->data_type() &&
      PEL_Data::type_name(to_eval->data_type())== t_name ;
   notify_one_test_result(nn+" kind", ok ) ;
   double epsilon = 1.0e-14 ;
   double min_d = 1.0e-14 ;
   if( texp->has_entry( "relative_precision" ) )
      epsilon = texp->double_data( "relative_precision" ) ;
   
   if( ok )
   {
      if( to_eval->data_type()==PEL_Data::Double )
      {
         ok = PEL::double_equality( texp->double_data(eval), res->to_double(),
                                    epsilon, min_d ) ;
      }
      else if( to_eval->data_type()==PEL_Data::Int )
      {
         ok = texp->int_data(eval) == res->to_int() ;
      }
      else if( to_eval->data_type()==PEL_Data::String )
      {
         ok = texp->string_data(eval) == res->to_string(res->context()) ;
      }
      else if( to_eval->data_type()==PEL_Data::Bool )
      {
         ok = texp->bool_data(eval) == res->to_bool() ;
      }
      else  if( to_eval->data_type()==PEL_Data::DoubleVector )
      {
         ok = ( texp->doubleVector_data(eval).size()==
                                         res->to_double_vector().size() ) ;
         for( size_t i=0 ; ok && i<texp->doubleVector_data(eval).size() ; i++ )
            ok = PEL::double_equality( texp->doubleVector_data(eval)(i),
                             res->to_double_vector()(i),
                                    epsilon, min_d ) ;
      }
      else  if( to_eval->data_type()==PEL_Data::IntVector )
      {
         ok = ( texp->intVector_data(eval).size()==
                                         res->to_int_vector().size() ) ;
         for( size_t i=0 ; ok && i<texp->intVector_data(eval).size() ; i++ )
            ok = ( texp->intVector_data(eval)(i)==res->to_int_vector()(i) ) ;
      }
      else  if( to_eval->data_type()==PEL_Data::BoolVector )
      {
         ok = ( texp->boolVector_data(eval).size()==
                                         res->to_bool_vector().size() ) ;
         for( size_t i=0 ; ok && i<texp->boolVector_data(eval).size() ; i++ )
            ok = ( texp->boolVector_data(eval)(i)==res->to_bool_vector()(i) ) ;
      }
      else  if( to_eval->data_type()==PEL_Data::StringVector )
      {
         ok = ( texp->stringVector_data(eval).size()==
                                         res->to_string_vector().size() ) ;
         for( size_t i=0 ; ok && i<texp->stringVector_data(eval).size() ; i++ )
            ok = ( texp->stringVector_data(eval)(i)==res->to_string_vector()(i) ) ;
      }
      else  if( to_eval->data_type()==PEL_Data::DoubleArray2D )
      {
         ok = ( texp->doubleArray2D_data(eval).index_bound(0)==
                             res->to_double_array2D().index_bound(0) &&
                texp->doubleArray2D_data(eval).index_bound(1)==
                             res->to_double_array2D().index_bound(1) ) ;
         for( size_t i=0 ; ok && i<texp->doubleArray2D_data(eval).index_bound(0) ; i++ )
            for( size_t j=0 ; ok && j<texp->doubleArray2D_data(eval).index_bound(1) ; j++ )
            ok = PEL::double_equality( texp->doubleArray2D_data(eval)(i,j),
                             res->to_double_array2D()(i,j),
                                    epsilon, min_d ) ;
      }
      else  if( to_eval->data_type()==PEL_Data::BoolArray2D )
      {
         ok = ( texp->boolArray2D_data(eval).index_bound(0)==
                             res->to_bool_array2D().index_bound(0) &&
                texp->boolArray2D_data(eval).index_bound(1)==
                             res->to_bool_array2D().index_bound(1) ) ;
         for( size_t i=0 ; ok && i<texp->boolArray2D_data(eval).index_bound(0) ; i++ )
            for( size_t j=0 ; ok && j<texp->boolArray2D_data(eval).index_bound(1) ; j++ )
               ok = texp->boolArray2D_data(eval)(i,j) == res->to_bool_array2D()(i,j) ;         
      }
      else  if( to_eval->data_type()==PEL_Data::StringArray2D )
      {
         ok = ( texp->stringArray2D_data(eval).index_bound(0)==
                             res->to_string_array2D().index_bound(0) &&
                texp->stringArray2D_data(eval).index_bound(1)==
                             res->to_string_array2D().index_bound(1) ) ;
         for( size_t i=0 ; ok && i<texp->stringArray2D_data(eval).index_bound(0) ; i++ )
            for( size_t j=0 ; ok && j<texp->stringArray2D_data(eval).index_bound(1) ; j++ )
            ok = texp->stringArray2D_data(eval)(i,j) == res->to_string_array2D()(i,j) ;
      }
      else  if( to_eval->data_type()==PEL_Data::DoubleArray3D )
      {
         ok = ( texp->doubleArray3D_data(eval).index_bound(0)==
                             res->to_double_array3D().index_bound(0) &&
                texp->doubleArray3D_data(eval).index_bound(1)==
                             res->to_double_array3D().index_bound(1) &&
                texp->doubleArray3D_data(eval).index_bound(2)==
                             res->to_double_array3D().index_bound(2) )
            ;
         for( size_t i=0 ; ok && i<texp->doubleArray3D_data(eval).index_bound(0) ; i++ )
            for( size_t j=0 ; ok && j<texp->doubleArray3D_data(eval).index_bound(1) ; j++ )
               for( size_t k=0 ; ok && k<texp->doubleArray3D_data(eval).index_bound(2) ; k++ )
            ok = PEL::double_equality( texp->doubleArray3D_data(eval)(i,j,k),
                             res->to_double_array3D()(i,j,k),
                                    epsilon, min_d  ) ;
      }
      else  if( to_eval->data_type()==PEL_Data::IntArray2D )
      {
         ok = ( texp->intArray2D_data(eval).index_bound(0)==
                             res->to_int_array2D().index_bound(0) &&
                texp->intArray2D_data(eval).index_bound(1)==
                             res->to_int_array2D().index_bound(1) ) ;
         for( size_t i=0 ; ok && i<texp->intArray2D_data(eval).index_bound(0) ; i++ )
            for( size_t j=0 ; ok && j<texp->intArray2D_data(eval).index_bound(1) ; j++ )
            ok = texp->intArray2D_data(eval)(i,j)==res->to_int_array2D()(i,j) ;
      }
      else  if( to_eval->data_type()==PEL_Data::IntArray3D )
      {
         ok = ( texp->intArray3D_data(eval).index_bound(0)==
                             res->to_int_array3D().index_bound(0) &&
                texp->intArray3D_data(eval).index_bound(1)==
                             res->to_int_array3D().index_bound(1) &&
                texp->intArray3D_data(eval).index_bound(2)==
                             res->to_int_array3D().index_bound(2) )
            ;
         for( size_t i=0 ; ok && i<texp->intArray3D_data(eval).index_bound(0) ; i++ )
            for( size_t j=0 ; ok && j<texp->intArray3D_data(eval).index_bound(1) ; j++ )
               for( size_t k=0 ; ok && k<texp->intArray3D_data(eval).index_bound(2) ; k++ )
            ok = texp->intArray3D_data(eval)(i,j,k)==res->to_int_array3D()(i,j,k) ;
      }
      else
      {
         PEL_Error::object()->raise_internal( "Unexpected kind of expression" ) ;
      }
      if( !ok && to_eval->data_type()==PEL_Data::Double )
      {
         std::cout.precision( 15 ) ;
         std::cout << "Expected : " << texp->double_data(eval) ;
         std::cout << std::endl ;
         std::cout << "Result : " << res->to_double() ;
         std::cout << std::endl ;
      }
      else if( !ok && to_eval->data_type()==PEL_Data::Int )
      {
         std::cout << "Expected : " << texp->int_data(eval) ;
         std::cout << std::endl ;
         std::cout << "Result : " << res->to_int() ;
         std::cout << std::endl ;
      }
      else if( !ok && to_eval->data_type()==PEL_Data::String )
      {
         std::cout << "Expected : " << texp->string_data(eval) ;
         std::cout << std::endl ;
         std::cout << "Result : " << res->to_string() ;
         std::cout << std::endl ;
      }
      else if( !ok && to_eval->data_type()==PEL_Data::Bool )
      {
         std::cout << "Expected : " << texp->bool_data(eval) ;
         std::cout << std::endl ;
         std::cout << "Result : " << res->to_bool() ;
         std::cout << std::endl ;
      }
      else if( !ok && to_eval->data_type()==PEL_Data::DoubleVector )
      {
         std::cout.precision( 15 ) ;
         std::cout << "Expected : " << texp->doubleVector_data(eval) ;
         std::cout << std::endl ;
         std::cout << "Result : " << res->to_double_vector() ;
         std::cout << std::endl ;
      }
      else if( !ok && to_eval->data_type()==PEL_Data::IntVector )
      {
         std::cout << "Expected : " << texp->intVector_data(eval) ;
         std::cout << std::endl ;
         std::cout << "Result : " << res->to_int_vector() ;
         std::cout << std::endl ;
      }
      else if( !ok && to_eval->data_type()==PEL_Data::StringVector )
      {
         std::cout << "Expected : " << texp->stringVector_data(eval) ;
         std::cout << std::endl ;
         std::cout << "Result : " << res->to_string_vector() ;
         std::cout << std::endl ;
      }
      else if( !ok && to_eval->data_type()==PEL_Data::BoolVector )
      {
         std::cout << "Expected : " << texp->boolVector_data(eval) ;
         std::cout << std::endl ;
         std::cout << "Result : " << res->to_bool_vector() ;
         std::cout << std::endl ;
      }
      else if( !ok && to_eval->data_type()==PEL_Data::DoubleArray2D )
      {
         std::cout.precision( 15 ) ;
         std::cout << "Expected : " << texp->doubleArray2D_data(eval) ;
         std::cout << std::endl ;
         std::cout << "Result : " << res->to_double_array2D() ;
         std::cout << std::endl ;
      }
      else if( !ok && to_eval->data_type()==PEL_Data::DoubleArray3D )
      {
         std::cout.precision( 15 ) ;
         std::cout << "Expected : " << texp->doubleArray3D_data(eval) ;
         std::cout << std::endl ;
         std::cout << "Result : " << res->to_double_array3D() ;
         std::cout << std::endl ;
      }
      else if( !ok && to_eval->data_type()==PEL_Data::IntArray2D )
      {
         std::cout << "Expected : " << texp->intArray2D_data(eval) ;
         std::cout << std::endl ;
         std::cout << "Result : " << res->to_int_array2D() ;
         std::cout << std::endl ;
      }
      else if( !ok && to_eval->data_type()==PEL_Data::IntArray3D )
      {
         std::cout << "Expected : " << texp->intArray3D_data(eval) ;
         std::cout << std::endl ;
         std::cout << "Result : " << res->to_int_array3D() ;
         std::cout << std::endl ;
      }
      notify_one_test_result(nn+" evaluation", ok ) ;
   }
}



