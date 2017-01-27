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

#include <GE_RefinedSegment_QR.hh>

#include <PEL_assertions.hh>
#include <PEL_Root.hh>
#include <PEL.hh>

#include <GE_Point.hh>
#include <GE_ReferenceSegment.hh>

//----------------------------------------------------------------------
GE_RefinedSegment_QR const*
GE_RefinedSegment_QR:: create( std::string a_name,
                               GE_QuadratureRule const* tria_rule )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_RefinedSegment_QR:: create" ) ;
   PEL_CHECK_PRE( tria_rule->reference_polyhedron() == 
                  GE_ReferenceSegment::object() ) ;

   GE_RefinedSegment_QR const* result = 
                         new GE_RefinedSegment_QR( a_name, tria_rule ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == PEL_Root::object() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_RefinedSegment_QR:: GE_RefinedSegment_QR( 
                               std::string a_name,
                               GE_QuadratureRule const* tria_rule )
//----------------------------------------------------------------------
   : GE_QuadratureRule( a_name, 
                        GE_ReferenceSegment::object(), 
                        tria_rule->order() )
{
   // left segment
   // ------------
   for( size_t ip = 0 ; ip<tria_rule->nb_points() ; ++ip )
   {
      double x = tria_rule->point( ip )->coordinate( 0 ) ;
      PEL_ASSERT( (x != 0.0) && (x != 1.0) ) ;

      GE_Point* p = GE_Point::create( this, 0.5*x ) ;

      append_point( p, 0.5*tria_rule->weight( ip ) ) ;
   }

   // right segment
   // ------------
   for( size_t ip = 0 ; ip<tria_rule->nb_points() ; ++ip )
   {
      double x = tria_rule->point( ip )->coordinate( 0 ) ;
      PEL_ASSERT( (x != 0.0) && (x != 1.0) ) ;

      GE_Point* p = GE_Point::create( this, 0.5*(1.0+x) ) ;

      append_point( p, 0.5*tria_rule->weight( ip ) ) ;
   }

   set_sum_of_weights( tria_rule->sum_of_weights() ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
GE_RefinedSegment_QR:: ~GE_RefinedSegment_QR( void )
//----------------------------------------------------------------------------
{
}
