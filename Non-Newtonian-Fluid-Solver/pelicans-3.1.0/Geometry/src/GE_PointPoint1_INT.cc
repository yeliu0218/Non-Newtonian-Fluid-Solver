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

#include <GE_PointPoint1_INT.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_Point.hh>

GE_PointPoint1_INT const*
GE_PointPoint1_INT::PROTOTYPE = new GE_PointPoint1_INT() ;

//----------------------------------------------------------------------
GE_PointPoint1_INT*
GE_PointPoint1_INT:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointPoint1_INT:: create_clone" ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   GE_PointPoint1_INT* result = new GE_PointPoint1_INT( a_owner,
                                                        EPSILON ) ;
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_PointPoint1_INT*
GE_PointPoint1_INT:: create( PEL_Object* a_owner, double a_epsilon )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointPoint1_INT:: create" ) ;
   PEL_CHECK_PRE( a_epsilon>0. ) ;
   
   GE_PointPoint1_INT* result = new GE_PointPoint1_INT( a_owner,
                                                        a_epsilon ) ;

   PEL_CHECK( result!=0 ) ;
   PEL_CHECK( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_PointPoint1_INT:: GE_PointPoint1_INT( void )
//----------------------------------------------------------------------
   : GE_PointPoint_INT( "GE_PointPoint1_INT" )
   , EPSILON( -PEL::max_double() )
{
   PEL_LABEL( "GE_PointPoint1_INT:: GE_PointPoint1_INT" ) ;
   PEL_CHECK( is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_PointPoint1_INT:: GE_PointPoint1_INT( PEL_Object* a_owner,
                                         double a_epsilon )
//----------------------------------------------------------------------
   : GE_PointPoint_INT( a_owner )
   , EPSILON( a_epsilon )
{
   PEL_LABEL( "GE_PointPoint1_INT:: GE_PointPoint1_INT" ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
GE_PointPoint1_INT:: ~GE_PointPoint1_INT( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointPoint1_INT:: ~GE_PointPoint1_INT" ) ;
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
bool
GE_PointPoint1_INT:: points_are_close( GE_Point const* P1,
                                        GE_Point const* P2 ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointPoint1_INT:: points_are_close" ) ;
   PEL_CHECK_PRE( points_are_close_PRE( P1, P2 ) ) ;
   PEL_CHECK( !is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = true ;
   size_t const nb_coords = P1->nb_coordinates() ;
   for( size_t i = 0 ; result && i<nb_coords ; ++i )
   {
      result = ( PEL::abs( P1->coordinate(i)-P2->coordinate(i) )<EPSILON ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
GE_PointPoint_INT*
GE_PointPoint1_INT:: create_replica(
                        PEL_Object* a_owner,
                        PEL_ModuleExplorer const* a_mod_exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PointPoint1_INT:: create_replica" ) ;
   PEL_CHECK( a_mod_exp!=0 ) ;
   PEL_CHECK( is_a_prototype() ) ;
   PEL_CHECK_INV( invariant() ) ;

   double const a_epsilon = a_mod_exp->double_data( "epsilon" ) ;
   a_mod_exp->test_data( "epsilon", "epsilon>0." ) ;
   GE_PointPoint_INT* result = new GE_PointPoint1_INT( a_owner,
                                                       a_epsilon ) ;

   PEL_CHECK( result!=0 ) ;
   PEL_CHECK( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
GE_PointPoint1_INT:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( GE_PointPoint_INT::invariant() ) ;
   return( true )  ;
}
