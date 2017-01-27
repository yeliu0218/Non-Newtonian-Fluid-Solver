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

#ifndef GE_SEGMENT_SEGMENT_INT_HH
#define GE_SEGMENT_SEGMENT_INT_HH

#include <PEL_Object.hh>

class GE_Point ;
class GE_PointPoint_INT ;
class GE_PointSegment_INT ;

class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;

/*
Algorithms investigating the intersection of two segments.

FRAMEWORK INSTANTIATION
   1. Derive a concrete subclass.
   2. Create a static class pointer which defines the model of the concrete
      class (prototype). It is build calling the prototype constructor :
          `GE_SegmentSegment_INT( std::string const& )'
   3. Implement the function
          `create_replica( PEL_Object*, PEL_ModuleExplorer const* ) const'
      which create a new instance of the concrete class, calling the constructor :
          `GE_SegmentSegment_INT( PEL_Object* a_owner )'
   4. Implement a destructor
   5. Implement the virtual function `create_clone'
   6. Implement the virtual functions which compute the intersection of two segments

PUBLISHED
*/

class PEL_EXPORT GE_SegmentSegment_INT : public PEL_Object
{
   public: //------------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return a new instance of type `a_name'.
      // If a `GE_PointPoint_INT::' object is needed, it is built with
      // the module GE_PointPoint_INT of `a_mod_exp' if any, else
      // `a_pt_pt_intersector' is taken.
      // If a `GE_PointSegment_INT::' object is needed, it is built with
      // the module GE_PointSegment_INT of `a_mod_exp' if any, else
      // `a_pt_seg_intersector' is taken.      
      static GE_SegmentSegment_INT* create(
                            PEL_Object* a_owner,
                            std::string const& a_name,
                            PEL_ModuleExplorer const* a_mod_exp,
                            GE_PointPoint_INT const* a_pt_pt_intersector=0,
                            GE_PointSegment_INT const* a_pt_seg_intersector=0 ) ;
      
      virtual GE_SegmentSegment_INT* create_clone( PEL_Object* a_owner ) const = 0 ;

   //-- Segments intersection

      // Intersection type :
      //    no_intersection : [P1P2] and [Q1Q2] do not intersect
      //    parallel : [P1P2] and [Q1Q2] are moreover strickly parallel
      //    colinear_disjoint : [P1P2] and [Q1Q2] are moreover colinear and disjoint
      //    colinear_one_in_the_other : [P1P2] and [Q1Q2] are moreover colinear and
      //                                one lies enterely in the other
      //    colinear : [P1P2] and [Q1Q2] are moreover colinear and each one contains
      //               one endpoint  of the other one
      //    one_intersection : [P1P2] and [Q1Q2] intersect in one point
      enum IntersectionType
      {
         no_intersection,
         parallel,
         colinear_disjoint,
         colinear_one_in_the_other,
         colinear,
         one_intersection
      } ;

      // Do [P1P2] and [Q1Q2] intersect ?
      virtual bool has_intersection( GE_Point const* P1,
                                     GE_Point const* P2,
                                     GE_Point const* Q1,
                                     GE_Point const* Q2 ) const = 0 ;
      
      // Compute the intersection of [P1P2] and [Q1Q2].
      virtual void compute_intersection( GE_Point const* P1,
                                         GE_Point const* P2,
                                         GE_Point const* Q1,
                                         GE_Point const* Q2 ) = 0 ;

      // intersection type
      virtual GE_SegmentSegment_INT::IntersectionType
                                      intersection_type( void ) const = 0 ;

      // alpha such that the intersection is expressed as
      // (1-alpha)*P1+(alpha)*P2
      virtual double alpha( void ) const = 0 ;

      // alpha such that the intersection is expressed as
      // (1-beta)*Q1+(beta)*Q2
      virtual double beta( void ) const = 0 ;
      
   protected: //----------------------------------------------------------------

      virtual ~GE_SegmentSegment_INT( void ) ;

      // Model registration constructor
      GE_SegmentSegment_INT( std::string const& a_name ) ;

      // Concrete constructor
      GE_SegmentSegment_INT( PEL_Object* a_owner ) ;

   //-- Internal status

      // is `self' a prototype
      bool is_a_prototype( void ) const ;      
      
   //-- Preconditions, Postconditions, Invariant
      
      virtual bool has_intersection_PRE( GE_Point const* P1,
                                         GE_Point const* P2,
                                         GE_Point const* Q1,
                                         GE_Point const* Q2 ) const ;
         
      virtual bool compute_intersection_PRE( GE_Point const* P1,
                                             GE_Point const* P2,
                                             GE_Point const* Q1,
                                             GE_Point const* Q2 ) const ;
      virtual bool compute_intersection_POST( GE_Point const* P1,
                                              GE_Point const* P2,
                                              GE_Point const* Q1,
                                              GE_Point const* Q2 ) const ;
      
      virtual bool alpha_PRE( void ) const ;
      virtual bool alpha_POST( double result ) const ;
      
      virtual bool beta_PRE( void ) const ;
      virtual bool beta_POST( double result ) const ;
      virtual bool invariant( void ) const ;
      
   private: //------------------------------------------------------------------

      GE_SegmentSegment_INT( void ) ;
      GE_SegmentSegment_INT( GE_SegmentSegment_INT const& other ) ;
      GE_SegmentSegment_INT const& operator=(
                             GE_SegmentSegment_INT const& other ) ;

      // Create replica of self from existing one :
      virtual GE_SegmentSegment_INT* create_replica( 
                     PEL_Object* a_owner,
                     PEL_ModuleExplorer const* a_mod_exp,
                     GE_PointPoint_INT const* a_pt_pt_intersector,
                     GE_PointSegment_INT const* a_pt_seg_intersector ) const = 0 ;

      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Attributes

      // Prototype :
      bool const PROTO ;      
      
} ;

#endif
