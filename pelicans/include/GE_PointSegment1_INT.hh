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

#ifndef GE_POINT_SEGMENT_1_INT_HH
#define GE_POINT_SEGMENT_1_INT_HH

#include <GE_PointSegment_INT.hh>

class GE_PointPoint_INT ;

/*

Concrete implementation of the intersection of a point and a segment.

A point P is in the segment [Q1,Q2] if
     - P is "close" to Q1 or "close" to "Q2" (GE_POINT_POINT_INT algorithm)
     - or PEL::abs( dQ1Q2-(dPQ1+dPQ2) )<epsilon
       with :
            dQ1Q2 = Q1->distance( Q2 ) ;
            dPQ1  = P->distance( Q1 ) ;
            dPQ2  = P->distance( Q2 ) ;
   
Parameters needed in the module explorer :
   - epsilon (positive double value) :
        defined previously
   - GE_PointPoint_INT (module) :
        module defining a `::GE_PointPoint_INT' algorithm

PUBLISHED
*/

class PEL_EXPORT GE_PointSegment1_INT : public GE_PointSegment_INT
{

   public: //------------------------------------------------------------------

   //-- Instance delivery and initialization

      static GE_PointSegment1_INT* create(
                        PEL_Object* a_owner,
                        GE_PointPoint_INT const* a_pt_pt_intersector,
                        double a_epsilon ) ;
      
      virtual  GE_PointSegment1_INT* create_clone( PEL_Object* a_owner ) const ;
      
   //-- Point and segment intersection

      // IMPLEMENTATION :
      //  A point P is in the segment [Q1,Q2] if
      //       - P is "closed" to Q1 or "closed" to "Q2"
      //       - or PEL::abs( dQ1Q2-(dPQ1+dPQ2) )<epsilon
      //            with :
      //                  dQ1Q2 = Q1->distance( Q2 ) ;
      //                  dPQ1  = P->distance( Q1 ) ;
      //                  dPQ2  = P->distance( Q2 ) ;
      virtual bool point_in_segment( GE_Point const* P,
                                     GE_Point const* Q1,
                                     GE_Point const* Q2 ) const ;
      
   protected: //---------------------------------------------------------------
      
   private: //-----------------------------------------------------------------

      GE_PointSegment1_INT( void ) ;
     ~GE_PointSegment1_INT( void ) ;
      GE_PointSegment1_INT( GE_PointSegment1_INT const& other ) ;
      GE_PointSegment1_INT const& operator=(
                            GE_PointSegment1_INT const& other ) ;

      GE_PointSegment1_INT( PEL_Object* a_owner,
                            GE_PointPoint_INT const* a_pt_pt_intersector,
                            double a_epsilon ) ;

      virtual GE_PointSegment_INT* create_replica( 
                     PEL_Object* a_owner,
                     PEL_ModuleExplorer const* a_mod_exp,
                     GE_PointPoint_INT const* a_pt_pt_intersector ) const ;
      
   //-- Preconditions, Postconditions, Invariant
      
      virtual bool invariant( void ) const ;
      
   //-- Class Attributes

      static GE_PointSegment1_INT const* PROTOTYPE ;
      
   //-- Attributes

      double const EPSILON ;
      GE_PointPoint_INT const* PT_PT_INTERSECTOR ;
      
} ;

#endif
