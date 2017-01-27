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

#ifndef GE_SEGMENT_SEGMENT_1_INT_HH
#define GE_SEGMENT_SEGMENT_1_INT_HH

#include <GE_SegmentSegment_INT.hh>

class GE_PointSegment_INT ;

/*
Concrete implementation of the intersection of two segment.

Parameters needed in the module explorer :
   - epsilon (positive double value) :
        absolute value of vectorial product of two segments under which it
        is set nil
   - GE_PointSegment (module) :
        module defining a `GE_PointSegment::' algorithm

PUBLISHED
*/

class PEL_EXPORT GE_SegmentSegment1_INT : public GE_SegmentSegment_INT
{

   public: //------------------------------------------------------------------

   //-- Instance delivery and initialization

      static GE_SegmentSegment1_INT* create(
                              PEL_Object* a_owner,
                              GE_PointSegment_INT const* a_pt_seg_intersector,
                              double a_epsilon ) ;
      
      virtual GE_SegmentSegment1_INT* create_clone( PEL_Object* a_owner ) const ;

   //-- Segments intersection

      virtual bool has_intersection( GE_Point const* P1,
                                     GE_Point const* P2,
                                     GE_Point const* Q1,
                                     GE_Point const* Q2 ) const ;
      
      virtual void compute_intersection( GE_Point const* P1,
                                         GE_Point const* P2,
                                         GE_Point const* Q1,
                                         GE_Point const* Q2 ) ;

      virtual GE_SegmentSegment_INT::IntersectionType
                                             intersection_type( void ) const ;

      virtual double alpha( void ) const ;

      virtual double beta( void ) const ;
      
   protected: //--------------------------------------------------------------
      
   private: //----------------------------------------------------------------

      
      GE_SegmentSegment1_INT( void ) ;
     ~GE_SegmentSegment1_INT( void ) ;
      GE_SegmentSegment1_INT( GE_SegmentSegment1_INT const& other ) ;
      GE_SegmentSegment1_INT const& operator=(
                              GE_SegmentSegment1_INT const& other ) ;

      GE_SegmentSegment1_INT( PEL_Object* a_owner,
                              GE_PointSegment_INT const* a_pt_seg_intersector,
                              double a_epsilon ) ;

      virtual GE_SegmentSegment_INT* create_replica( 
                     PEL_Object* a_owner,
                     PEL_ModuleExplorer const* a_mod_exp,
                     GE_PointPoint_INT const* a_pt_pt_intersector,
                     GE_PointSegment_INT const* a_pt_seg_intersector ) const ;
      
   //-- Preconditions, Postconditions, Invariant
      
      virtual bool invariant( void ) const ;
      
   //-- Class Attributes

      static GE_SegmentSegment1_INT const* PROTOTYPE ;
      
   //-- Attributes

      double const EPSILON ;
      GE_PointSegment_INT const* const PT_SEG_INT ;
      
      GE_SegmentSegment_INT::IntersectionType INTER_TYPE ;
      double ALPHA ;
      double BETA ;
} ;

#endif
