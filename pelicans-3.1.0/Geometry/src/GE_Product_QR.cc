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

#include <GE_Product_QR.hh>

#include <PEL_assertions.hh>
#include <PEL_Root.hh>
#include <PEL.hh>

#include <GE_Point.hh>
#include <GE_ReferenceSegment.hh>
#include <GE_ReferenceSquare.hh>
#include <GE_ReferenceCube.hh>

//----------------------------------------------------------------------
GE_Product_QR const*
GE_Product_QR:: create( std::string a_name,
                        GE_QuadratureRule const* rule1D,
                        GE_ReferencePolyhedron const* poly )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Product_QR:: create" ) ;
   PEL_CHECK_PRE( rule1D->reference_polyhedron() == 
                  GE_ReferenceSegment::object() ) ;
   PEL_CHECK_PRE( poly == GE_ReferenceSquare::object() ||
                  poly == GE_ReferenceCube::object() ) ;

   GE_Product_QR const* result = new GE_Product_QR( a_name, rule1D, poly ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == PEL_Root::object() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Product_QR:: GE_Product_QR( std::string a_name,
                               GE_QuadratureRule const* rule1D,
                               GE_ReferencePolyhedron const* poly )
//----------------------------------------------------------------------
   : GE_QuadratureRule( a_name, poly, rule1D->order() )
{
   size_t nb_pt = rule1D->nb_points() ;
   double sum = rule1D->sum_of_weights() ;

   if( poly == GE_ReferenceSquare::object() )
   {
      for( size_t i=0; i<nb_pt; ++i )
      {
         for( size_t j=0; j<nb_pt; ++j )
         {
            GE_Point* p = GE_Point::create( this,
     	                                    rule1D->point(i)->coordinate(0),
     	                                    rule1D->point(j)->coordinate(0) ) ;
            double w = rule1D->weight(i) * rule1D->weight(j) ;
            append_point( p, w ) ;
         }
      }
      set_sum_of_weights( sum*sum ) ;
   }
   else if( poly == GE_ReferenceCube::object() )
   {
      for( size_t i=0; i<nb_pt; ++i )
      {
         for( size_t j=0; j<nb_pt; ++j )
         {
            for( size_t k=0; k<nb_pt; ++k )
            {
               GE_Point* p = GE_Point::create( this,
     	                                    rule1D->point(i)->coordinate(0),
     	                                    rule1D->point(j)->coordinate(0),
     	                                    rule1D->point(k)->coordinate(0) ) ;
               double w = rule1D->weight(i) * 
                          rule1D->weight(j) * 
                          rule1D->weight(k) ;
               append_point( p, w ) ;
            }
         }
      }
      set_sum_of_weights( sum*sum*sum ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
GE_Product_QR:: ~GE_Product_QR( void )
//----------------------------------------------------------------------------
{
}
