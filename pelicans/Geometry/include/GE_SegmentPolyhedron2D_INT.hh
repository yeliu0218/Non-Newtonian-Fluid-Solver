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

#ifndef GE_SEGMENT_POLYHEDRON_2D_INT_HH
#define GE_SEGMENT_POLYHEDRON_2D_INT_HH

#include <GE_SegmentPolyhedron_INT.hh>

class GE_Point ;

/*
Servers for segment to plane 2D polyhedron intersection checking.

PUBLISHED
*/

class PEL_EXPORT GE_SegmentPolyhedron2D_INT : public GE_SegmentPolyhedron_INT
{
   public: //------------------------------------------------------------------

   //-- Intersection

      virtual void check_intersection( GE_Point const* S0, 
                                       GE_Point const* S1, 
                                       GE_Mpolyhedron const* M ) ;

      virtual bool one_single_intersection( void ) const ;

      virtual void intersection_point( GE_Point* pt ) const ; 

   protected: //----------------------------------------------------------------

   private: //------------------------------------------------------------------

     ~GE_SegmentPolyhedron2D_INT( void ) ;
      GE_SegmentPolyhedron2D_INT( GE_SegmentPolyhedron2D_INT const& other ) ;
      GE_SegmentPolyhedron2D_INT& operator=(
                                  GE_SegmentPolyhedron2D_INT const& other ) ;
      
   //-- Plug in
      
      GE_SegmentPolyhedron2D_INT( void ) ;
      
      GE_SegmentPolyhedron2D_INT( PEL_Object* a_owner,
                                  PEL_ModuleExplorer const* a_mod_exp ) ;

      virtual GE_SegmentPolyhedron_INT* create_replica( 
                                  PEL_Object* a_owner,
                                  PEL_ModuleExplorer const* a_mod_exp ) const ;

  //-- Class Attributes

      void check_intersection( GE_Point const* S0,
                               GE_Point const* S1,
                               GE_Point const* V0,
                               GE_Point const* V1,
                               GE_Point const* V2 ) ;
      
      static GE_SegmentPolyhedron_INT const* PROTOTYPE ;

   //-- Attributes

      bool ONE_INTER ;
      GE_Point* const PI ;
} ;

#endif
