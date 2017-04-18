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

#ifndef GE_SEGMENT_POLYHEDRON_INT_HH
#define GE_SEGMENT_POLYHEDRON_INT_HH

#include <PEL_Object.hh>

class GE_Mpolyhedron ;
class GE_Point ;

class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;

/*
Servers for checking intersection between segment and  polyhedron
of dimension the space dimension minus one.

FRAMEWORK INSTANTIATION
   1. Derive a concrete subclass.
   2. Create a static class pointer which defines the model of the concrete
      class (prototype). It is build calling the prototype constructor :
          `GE_SegmentPolyhedron_INT( std::string const&, size_t )'
   3. Implement the function
          `create_replica( PEL_Object*, size_t, PEL_ModuleExplorer const* ) const'
      which create a new instance of the concrete class, calling the constructor :
          `GE_SegmentPolyhedron_INT( PEL_Object* a_owner, size_t )'
   4. Implement a destructor
   5. Implement other virtual functions

PUBLISHED
*/

class PEL_EXPORT GE_SegmentPolyhedron_INT : public PEL_Object
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      static GE_SegmentPolyhedron_INT* make( PEL_Object* a_owner,
                                             std::string const& a_name,
                                             size_t nb_space_dim,
                                             PEL_ModuleExplorer const* a_mod_exp ) ;

   //-- Status

      size_t nb_space_dimensions( void ) const ;

   //-- Intersection
      
      // Check intersection for segment S0S1 and polyhedron target M.
      virtual void check_intersection( GE_Point const* S0, 
                                       GE_Point const* S1, 
                                       GE_Mpolyhedron const* M ) = 0 ;

      // Is intersection performed ?
      virtual bool intersection_checked( void ) const ;

      // Is the segment cuts the polyhedron in one single point ?
      virtual bool one_single_intersection( void ) const = 0 ;

      // Set point `pt' as intersection point.
      virtual void intersection_point( GE_Point* pt ) const = 0 ;
      
      GE_Point const* segment_first_vertex( void ) const ;

      GE_Point const* segment_second_vertex( void ) const ;

      GE_Mpolyhedron const* target_polyhedron( void ) const ;

   //--  Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //---------------------------------------------------------------

   //-- Plug in
      
      virtual ~GE_SegmentPolyhedron_INT( void ) ;

      // Registration of an instance as `a_name'.
      GE_SegmentPolyhedron_INT( std::string const& a_name,
                                size_t nb_space_dim ) ;

      // Constructor called by `::create_replica'.
      GE_SegmentPolyhedron_INT( PEL_Object* a_owner,
                                size_t nb_space_dim ) ;

      // Create replica of self from existing one.
      virtual GE_SegmentPolyhedron_INT* create_replica(
                         PEL_Object* a_owner,
                         PEL_ModuleExplorer const* a_mod_exp ) const = 0 ;
      
  //-- Intersection
      
      virtual void reset( GE_Point const* S0, 
                          GE_Point const* S1, 
                          GE_Mpolyhedron const* M ) ;

      virtual void declare_intersection_checked( void ) ;

   //-- Internal status
      
      bool is_a_prototype( void ) const ;

  //-- Preconditions, Postconditions, Invariant
      
      bool check_intersection_PRE( GE_Point const* S0,
                                   GE_Point const* S1,
				   GE_Mpolyhedron const* M ) const ;
      bool check_intersection_POST( GE_Point const* S0,
                                    GE_Point const* S1,
				    GE_Mpolyhedron const* M ) const ;

      bool create_replica_PRE( PEL_Object const* a_owner,
                               PEL_ModuleExplorer const* a_mod_exp ) const ;
      bool create_replica_POST( GE_SegmentPolyhedron_INT const* result,
                                PEL_Object const* a_owner,
                                PEL_ModuleExplorer const* a_mod_exp ) const ;

      bool intersection_point_PRE( GE_Point const* pt ) const ; 

      virtual bool invariant( void ) const ;

   private: //-----------------------------------------------------------------

      GE_SegmentPolyhedron_INT( void ) ;
      GE_SegmentPolyhedron_INT( GE_SegmentPolyhedron_INT const& other ) ;
      GE_SegmentPolyhedron_INT& operator=(
                                GE_SegmentPolyhedron_INT const& other ) ;

      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Attributes
      
      bool const PROTO ;
      size_t const DIM ;
      
      bool INTER_CHECKED ;
      GE_Mpolyhedron const* M_SAVE ;
      GE_Point const* S0_SAVE ;
      GE_Point const* S1_SAVE ;
} ;

#endif
