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

#include <vectorIO_TEST.hh>

#include <PEL.hh>
#include <PEL_Exceptions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <boolVector.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <size_t_vector.hh>
#include <stringVector.hh>

#include <iostream>
#include <sstream>

vectorIO_TEST const*
vectorIO_TEST::PROTOTYPE = new vectorIO_TEST() ;

//-------------------------------------------------------------------------
vectorIO_TEST:: vectorIO_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "vector", "vectorIO_TEST" )
{
}

//-------------------------------------------------------------------------
vectorIO_TEST:: ~vectorIO_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
vectorIO_TEST:: process_one_test( PEL_ModuleExplorer const* texp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "vectorIO_TEST:: process_one_test" ) ;
   
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
   if( request_exception || catched )
   {
      std::string const& nn = texp->name() ;
      notify_one_test_result(nn+" exception", catched==request_exception ) ;
   }
}

//-------------------------------------------------------------------------
void
vectorIO_TEST:: process_one_comparison( PEL_ModuleExplorer const* texp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "vectorIO_TEST:: process_one_comparison" ) ;
   
   std::string const& nn = texp->name() ;

   std::string const& t = texp->string_data( "type" ) ;
   
   std::string const& exp_to_eval = texp->string_data( "exp_to_eval" ) ;
   
   std::stringstream m ;
   m << exp_to_eval ;
   
   if( t == "intVector" )
   {
      intVector v(0) ;
      m >> v ;
      intVector const& result = texp->intVector_data( "result" ) ;
      bool ok = ( v == result ) ;
      notify_one_test_result( nn+" read", ok ) ;
      if( !ok )
      {
         out() << "   read:" << std::endl << v << std::endl ;
         out() << "   expected:" << std::endl << result << std::endl ;
      }

      std::stringstream m2 ;
      m2 << result ;
      m2 >> v ;
      ok = ( v == result ) ;
      notify_one_test_result( nn+" write_read", ok ) ;
      if( !ok )
      {
         out() << "   read:" << std::endl << v << std::endl ;
         out() << "   expected:" << std::endl << result << std::endl ;
      }
   }
   else if( t == "size_t_vector" )
   {
      size_t_vector v(0) ;
      m >> v ;
      intVector const& r = texp->intVector_data( "result" ) ;
      size_t_vector result( r ) ;
      bool ok = ( v == result ) ;
      notify_one_test_result( nn+" read", ok ) ;
      if( !ok )
      {
         out() << "   read:" << std::endl << v << std::endl ;
         out() << "   expected:" << std::endl << result << std::endl ;
      }

      std::stringstream m2 ;
      m2 << result ;
      m2 >> v ;
      ok = ( v == result ) ;
      notify_one_test_result( nn+" write_read", ok ) ;
      if( !ok )
      {
         out() << "   read:" << std::endl << v << std::endl ;
         out() << "   expected: << std::endl" << result << std::endl ;
      }
   }
   else if( t == "boolVector" )
   {
      boolVector v(0) ;
      m >> v ;
      boolVector const& result = texp->boolVector_data( "result" ) ;
      bool ok = ( v == result ) ;
      notify_one_test_result( nn+" read", ok ) ;
      if( !ok )
      {
         out() << "   read:" << std::endl << v << std::endl ;
         out() << "   expected:" << std::endl << result << std::endl ;
      }

      std::stringstream m2 ;
      m2 << result ;
      m2 >> v ;
      ok = ( v == result ) ;
      notify_one_test_result( nn+" write_read", ok ) ;
      if( !ok )
      {
         out() << "   read:" << std::endl << v << std::endl ;
         out() << "   expected:" << std::endl << result << std::endl ;
      }
   }
   else if( t == "stringVector" )
   {
      stringVector v(0) ;
      m >> v ;
      stringVector const& result = texp->stringVector_data( "result" ) ;
      bool ok = ( v == result ) ;
      notify_one_test_result( nn+" read", ok ) ;
      if( !ok )
      {
         out() << "   read:" << std::endl << v << std::endl ;
         out() << "   expected:" << std::endl << result << std::endl ;
      }

      std::stringstream m2 ;
      m2 << result ;
      m2 >> v ;
      ok = ( v == result ) ;
      notify_one_test_result( nn+" write_read", ok ) ;
      if( !ok )
      {
         out() << "   read:" << std::endl << v << std::endl ;
         out() << "   expected:" << std::endl << result << std::endl ;
      }
   }
   else if( t == "doubleVector" )
   {
      double dbl_eps = 1.0e-14 ;
      double dbl_min = 1.0e-14 ;
      if( texp->has_entry( "dbl_epsilon" ) )
      {
         dbl_eps = texp->double_data( "dbl_epsilon" ) ;
         texp->test_data( "dbl_epsilon", "dbl_epsilon>0." ) ;
      }
      if( texp->has_entry( "dbl_minimal" ) )
      {
         dbl_min = texp->double_data( "dbl_minimal" ) ;
         texp->test_data( "dbl_minimal", "dbl_minimal>0." ) ;
      }
      
      doubleVector v(0) ;
      m >> v ;
      doubleVector const& result = texp->doubleVector_data( "result" ) ;
      bool ok = doubleVector_equality( v, result, dbl_eps, dbl_min ) ;
      notify_one_test_result( nn+" read", ok ) ;
      if( !ok )
      {
         out() << "   read:" << std::endl << v << std::endl ;
         out() << "   expected:" << std::endl << result << std::endl ;
      }

      std::stringstream m2 ;
      m2 << result ;
      m2 >> v ;
      ok = doubleVector_equality( v, result, dbl_eps, dbl_min ) ;
      notify_one_test_result( nn+" write_read", ok ) ;
      if( !ok )
      {
         out() << "   read:" << std::endl << v << std::endl ;
         out() << "   expected:" << std::endl << result << std::endl ;
      }
   }
   else
   {
      PEL_Error::object()->raise_internal( "Unexpected kind of expression" ) ;
   }
}

//-------------------------------------------------------------------------
bool
vectorIO_TEST:: doubleVector_equality( doubleVector const& v1,
                                             doubleVector const& v2,
                                             double dbl_eps,
                                             double dbl_min ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "vectorIO_TEST:: doubleVector_equality" ) ;

   bool result = ( v1.size() == v2.size() ) ;
   for( size_t i=0 ; result && i<v1.size() ; ++i )
   {
      result = PEL::double_equality( v1(i), v2(i), dbl_eps, dbl_min ) ;
   }
   return( result ) ;
}
