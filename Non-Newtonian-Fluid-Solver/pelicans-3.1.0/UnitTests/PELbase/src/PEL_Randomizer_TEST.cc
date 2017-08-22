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

#include <PEL_Randomizer_TEST.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <PEL_Randomizer.hh>
#include <intVector.hh>
#include <doubleVector.hh>

#include <iostream>

PEL_Randomizer_TEST const*
PEL_Randomizer_TEST::PROTOTYPE = new PEL_Randomizer_TEST() ;

//-------------------------------------------------------------------------
PEL_Randomizer_TEST:: PEL_Randomizer_TEST( void )
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "PEL_Randomizer", "PEL_Randomizer_TEST" )
{
}

//-------------------------------------------------------------------------
PEL_Randomizer_TEST:: ~PEL_Randomizer_TEST( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
PEL_Randomizer_TEST:: process_one_test( PEL_ModuleExplorer const* texp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Randomizer_TEST:: process_one_test" ) ;

   std::string const& t = texp->string_data( "type" ) ;
   std::string const& nn = texp->name() ;
   if( t == "randomizer" )
   {
      int const serie = texp->int_data( "serie" ) ;
      doubleVector const& resu = texp->doubleVector_data( "result" ) ;
      bool const verbose = texp->bool_data( "verbose" ) ;
      PEL_Randomizer* r = PEL_Randomizer::create( 0, serie ) ;
      bool ok = true ;
      for( size_t k=0 ; k<resu.size() ; k++ )
      {
         ok &= PEL::toler( r->item()-resu(k), 1.0e-6 ) ;
         if( verbose ) out() <<"Rand( "<<k<<" )=" << r->item() << std::endl ;
         r->go_next() ;
      }
      notify_one_test_result( nn+":: randomizer", ok ) ;
      r->destroy() ; r = 0 ;
   }
   else if( t == "random_double" )
   {
      size_t const nb_calls = texp->int_data( "N" ) ;
      doubleVector val( nb_calls ) ;
      texp->test_data( "N", "N>0" ) ;
      for( size_t i=0 ; i<nb_calls ; ++i )
      {
         val(i) = texp->double_data( "value" ) ;
      }
      bool const verbose = texp->bool_data( "verbose" ) ;
      if( verbose ) out() << val << std::endl ;
      bool ok = true ;
      for( size_t i=0 ; ok && i<nb_calls ; ++i )
      {
         ok &= ( val(i)>0. && val(i)<1. ) ;
         for( size_t j=i+1 ; ok && j<nb_calls ; ++j )
         {
            ok &= ( val(i) != val(j) ) ;
         }
      }
      notify_one_test_result( nn+":: random_double", ok ) ;
   }
   else if( t == "rand" )
   {
      size_t const nb_calls = texp->int_data( "N" ) ;
      intVector val( nb_calls ) ;
      texp->test_data( "N", "N>0" ) ;
      for( size_t i=0 ; i<nb_calls ; ++i )
      {
         val(i) = texp->int_data( "value" ) ;
      }
      bool const verbose = texp->bool_data( "verbose" ) ;
      if( verbose ) out() << val << std::endl ;
      bool ok = true ;
      for( size_t i=0 ; ok && i<nb_calls ; ++i )
      {
         ok &= ( val(i)>0 ) ;
         for( size_t j=i+1 ; ok && j<nb_calls ; ++j )
         {
            ok &= ( val(i) != val(j) ) ;
         }
      }
      notify_one_test_result( nn+":: rand", ok ) ;
   }
}
