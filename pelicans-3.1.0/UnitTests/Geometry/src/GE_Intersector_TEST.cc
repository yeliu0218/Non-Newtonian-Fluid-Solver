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

#include <GE_Intersector_TEST.hh>

#include <GE_Point.hh>
#include <GE_PointPoint_INT.hh>
#include <GE_PointSegment_INT.hh>
#include <GE_SegmentSegment_INT.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <doubleVector.hh>

#include <iostream>
#include <string>
#include <sstream>

//-------------------------------------------------------------------------
GE_Intersector_TEST*
GE_Intersector_TEST::unique_instance = new GE_Intersector_TEST() ;
//-------------------------------------------------------------------------



//-------------------------------------------------------------------------
GE_Intersector_TEST:: GE_Intersector_TEST( void ) 
//-------------------------------------------------------------------------
   : PEL_ObjectTest( "GE_Intersector", "GE_Intersector_TEST" )
{
   PEL_LABEL( "GE_Intersector_TEST:: GE_Intersector_TEST" ) ;
   PEL_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
GE_Intersector_TEST:: ~GE_Intersector_TEST( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Intersector_TEST:: ~GE_Intersector_TEST" ) ;
   PEL_CHECK_INV( invariant() ) ;
}





//-------------------------------------------------------------------------
void
GE_Intersector_TEST:: process_one_test( PEL_ModuleExplorer const* t_exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Intersector_TEST:: process_all_tests" ) ;

   // Build intersectors :
   GE_PointPoint_INT* pt_pt_algo = 0 ;
   GE_PointSegment_INT* pt_seg_algo = 0 ;
   GE_SegmentSegment_INT* seg_seg_algo = 0 ;
   if( t_exp->has_module( "GE_PointPoint_INT" ) )
   {
      PEL_ModuleExplorer const* pt_pt_exp =
         t_exp->create_subexplorer( this, "GE_PointPoint_INT" ) ;
      std::string const& pt_pt_name =
         pt_pt_exp->string_data( "concrete_name" ) ;
      pt_pt_algo = GE_PointPoint_INT::create( this,
                                              pt_pt_name,
                                              pt_pt_exp ) ;
   }
   if( t_exp->has_module( "GE_PointSegment_INT" ) )
   {
      PEL_ModuleExplorer const* pt_seg_exp =
         t_exp->create_subexplorer( this, "GE_PointSegment_INT" ) ;
      std::string const& pt_seg_name =
         pt_seg_exp->string_data( "concrete_name" ) ;
      pt_seg_algo = GE_PointSegment_INT::create( this,
                                                 pt_seg_name,
                                                 pt_seg_exp,
                                                 pt_pt_algo ) ;
   }
   if( t_exp->has_module( "GE_SegmentSegment_INT" ) )
   {
      PEL_ModuleExplorer const* seg_seg_exp =
         t_exp->create_subexplorer( this, "GE_SegmentSegment_INT" ) ;
      std::string const& seg_seg_name =
         seg_seg_exp->string_data( "concrete_name" ) ;
      seg_seg_algo = GE_SegmentSegment_INT::create( this,
                                                    seg_seg_name,
                                                    seg_seg_exp,
                                                    pt_pt_algo,
                                                    pt_seg_algo ) ;
   }

   // Tests explorer :
   PEL_ModuleExplorer* tests_exp =
      t_exp->create_subexplorer( this, "TESTS" ) ;

   // Perform tests :
   for( tests_exp->start_module_iterator() ;
        tests_exp->is_valid_module() ;
        tests_exp->go_next_module() )
   {
      PEL_ModuleExplorer* exp = tests_exp->create_subexplorer( this ) ;
      std::string const& test_type = exp->string_data( "type" ) ;
      if( test_type=="point_point_intersection" )
      {
         test_pt_pt_intersector( pt_pt_algo,
                                 t_exp->name()+"/"+exp->name(),
                                 exp ) ;
      }
      else if( test_type=="point_segment_intersection" )
      {
         test_pt_seg_intersector( pt_seg_algo,
                                  t_exp->name()+"/"+exp->name(),
                                  exp ) ;
      }
      else if( test_type=="segment_segment_intersection" )
      {
         test_seg_seg_intersector( seg_seg_algo,
                                   t_exp->name()+"/"+exp->name(),
                                   exp ) ;
      }
      else
      {
         std::ostringstream allowed_values ;
         allowed_values << "  \"point_point_intersection\"\n"
                        << "  \"point_segment_intersection\"\n"
                        << "  \"segment_segment_intersection\"\n" ;
         PEL_Error::object()->raise_bad_data_value( exp,
                                                    "type",
                                                    allowed_values.str() ) ;
      }
   }
      

   PEL_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
void
GE_Intersector_TEST:: test_pt_pt_intersector(
                                   GE_PointPoint_INT* pt_pt_algo,
                                   std::string const& test_name,
                                   PEL_ModuleExplorer const* test_exp )
//-------------------------------------------------------------------------
{
   bool const result = test_exp->bool_data( "result" ) ;
   doubleVector p1_coord =  test_exp->doubleVector_data( "P1" ) ;
   GE_Point const* P1 = GE_Point::create( 0, p1_coord ) ;
   doubleVector p2_coord =  test_exp->doubleVector_data( "P2" ) ;
   GE_Point const* P2 = GE_Point::create( 0, p2_coord ) ;
   
   bool ok_test = ( pt_pt_algo->points_are_close( P1, P2 )==result ) ;
   
   notify_one_test_result( test_name, ok_test ) ;

   P1->destroy() ; P1 = 0 ;
   P2->destroy() ; P2 = 0 ;
}



//-------------------------------------------------------------------------
void
GE_Intersector_TEST:: test_pt_seg_intersector(
                                   GE_PointSegment_INT* pt_seg_algo,
                                   std::string const& test_name,
                                   PEL_ModuleExplorer const* test_exp )
//-------------------------------------------------------------------------
{
   bool const result = test_exp->bool_data( "result" ) ;
   doubleVector p_coord =  test_exp->doubleVector_data( "P" ) ;
   GE_Point const* P = GE_Point::create( 0, p_coord ) ;
   doubleVector q1_coord =  test_exp->doubleVector_data( "Q1" ) ;
   GE_Point const* Q1 = GE_Point::create( 0, q1_coord ) ;
   doubleVector q2_coord =  test_exp->doubleVector_data( "Q2" ) ;
   GE_Point const* Q2 = GE_Point::create( 0, q2_coord ) ;
   
   bool ok_test = ( pt_seg_algo->point_in_segment( P, Q1, Q2 )==result ) ;
   
   notify_one_test_result( test_name, ok_test ) ;

   P->destroy() ; P = 0 ;
   Q1->destroy() ; Q1 = 0 ;
   Q2->destroy() ; Q2 = 0 ;
}



//-------------------------------------------------------------------------
void
GE_Intersector_TEST:: test_seg_seg_intersector(
                                   GE_SegmentSegment_INT* seg_seg_algo,
                                   std::string const& test_name,
                                   PEL_ModuleExplorer const* test_exp )
//-------------------------------------------------------------------------
{
   PEL_ModuleExplorer* rexp = test_exp->create_subexplorer( this, "result" ) ;
   
   std::string const& result_s = rexp->string_data( "type" ) ;
   GE_SegmentSegment_INT::IntersectionType result =
                                   GE_SegmentSegment_INT::no_intersection ;
   if( result_s=="no_intersection" )
   {
      result = GE_SegmentSegment_INT::no_intersection ;
   }
   else if( result_s=="parallel" )
   {
      result = GE_SegmentSegment_INT::parallel ;
   }
   else if( result_s=="colinear_disjoint" )
   {
      result = GE_SegmentSegment_INT::colinear_disjoint ;
   }
   else if( result_s=="colinear_one_in_the_other" )
   {
      result = GE_SegmentSegment_INT::colinear_one_in_the_other ;
   }
   else if( result_s=="colinear" )
   {
      result = GE_SegmentSegment_INT::colinear ;
   }
   else if( result_s=="one_intersection" )
   {
      result = GE_SegmentSegment_INT::one_intersection ;
   }
   else
   {
      std::ostringstream allowed_values ;
      allowed_values << "  \"no_intersection\"\n"
                     << "  \"parallel\"\n"
                     << "  \"colinear_disjoint\"\n"
                     << "  \"colinear_one_in_the_other\"\n"
                     << "  \"colinear\"\n"
                     << "  \"one_intersection\n" ;
      PEL_Error::object()->raise_bad_data_value( rexp,
                                                 "type",
                                                 allowed_values.str() ) ;
   }

   doubleVector p1_coord =  test_exp->doubleVector_data( "P1" ) ;
   GE_Point const* P1 = GE_Point::create( 0, p1_coord ) ;
   doubleVector p2_coord =  test_exp->doubleVector_data( "P2" ) ;
   GE_Point const* P2 = GE_Point::create( 0, p2_coord ) ;
   doubleVector q1_coord =  test_exp->doubleVector_data( "Q1" ) ;
   GE_Point const* Q1 = GE_Point::create( 0, q1_coord ) ;
   doubleVector q2_coord =  test_exp->doubleVector_data( "Q2" ) ;
   GE_Point const* Q2 = GE_Point::create( 0, q2_coord ) ;

   // has intersection :
   bool has_int_test =
      EQUIVALENT( seg_seg_algo->has_intersection( P1, P2, Q1, Q2 ),
                  result==GE_SegmentSegment_INT::one_intersection ) ;
   notify_one_test_result( test_name+"/has_intersection", has_int_test ) ;

   // intersection type :
   seg_seg_algo->compute_intersection( P1, P2, Q1, Q2 ) ;
   bool int_type_test = seg_seg_algo->intersection_type()==result ;
   notify_one_test_result( test_name+"/intersection_type", int_type_test ) ;

   // intersection :
   if( result==GE_SegmentSegment_INT::one_intersection )
   {
      double const alpha = rexp->double_data( "alpha" ) ;
      bool alpha_test = PEL::equal( alpha, seg_seg_algo->alpha() ) ;
      if( !alpha_test )
      {
         std::cout << "   alpha computed : " << seg_seg_algo->alpha()
                   << std::endl ;
         std::cout << "   alpha expected : " << alpha
                   << std::endl ;
      }
      notify_one_test_result( test_name+"/alpha", alpha_test ) ;
      double const beta = rexp->double_data( "beta" ) ;
      bool beta_test = PEL::equal( beta, seg_seg_algo->beta() ) ;
      if( !beta_test )
      {
         std::cout << "   beta computed : " << seg_seg_algo->beta()
                   << std::endl ;
         std::cout << "   beta expected : " << beta
                   << std::endl ;
      }
      notify_one_test_result( test_name+"/beta", beta_test ) ;         
   }

   P1->destroy() ; P1 = 0 ;
   P2->destroy() ; P2 = 0 ;
   Q1->destroy() ; Q1 = 0 ;
   Q2->destroy() ; Q2 = 0 ;   
}

