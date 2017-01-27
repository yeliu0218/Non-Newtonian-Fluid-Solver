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

#ifndef PDE_POINT_IN_GRID_FE_HH
#define PDE_POINT_IN_GRID_FE_HH

#include <PEL_Object.hh>

class PEL_ModuleExplorer ;
class boolVector ;

class GE_Point ;
class GE_SegmentPolyhedron_INT ;

class PDE_CursorFEside ;
class PDE_DomainAndFields ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;

/*
servers devoted to the problem of point location in a grid

Each instance is attached to a `PDE_DomainAndFields::' object, say `dom', 
and provides the service that identifies which cell of the grid defined
in `dom', if any, contains a given query point.

PUBLISHED
*/

class PEL_EXPORT PDE_PointInGridFE : public PEL_Object
{

   public: //-----------------------------------------------------------
      
   //-- Instance delivery and initialization
      
      // Create and return an instance attached to `dom'.
      static PDE_PointInGridFE* create( PEL_Object* a_owner,
                                        PDE_DomainAndFields const* dom,
                                        PEL_ModuleExplorer const* exp = 0 ) ;

   //-- Grid intersection

      size_t nb_space_dimensions( void ) const ;
      
      // Is the query point `pt' interior to the grid of the attached 
      // `PDE_DomainAndFields::' object (the search being first
      // performed in the current cell of `cFE') ?
      bool is_in_grid( GE_Point const* pt, PDE_LocalFEcell* cFE ) const ; 

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_PointInGridFE( void ) ;
     ~PDE_PointInGridFE( void ) ;
      PDE_PointInGridFE( PDE_PointInGridFE const& other ) ;
      PDE_PointInGridFE& operator=( PDE_PointInGridFE const& other ) ;

      PDE_PointInGridFE( PEL_Object* a_owner,
                         PDE_DomainAndFields const* dom,
                         PEL_ModuleExplorer const* exp ) ;

      size_t cell_neighbour_id( GE_Point const* pt,
                                PDE_LocalFEcell const* cFE,
                                boolVector const& visited_cells ) const ;

      size_t cell_id_after_hole( GE_Point const* pt,
                                 PDE_LocalFEcell const* cFE,
                                 boolVector const& visited_cells ) const ;

      // Complete search : test all cells of the grid.
      void search_in_all_cells( GE_Point const* pt, PDE_LocalFEcell* cFE,
                                bool& found ) const ;
      
      static bool is_visited( boolVector const& visited_cells,
                              size_t id_cell ) ;
      static void set_visited( boolVector& visited_cells,
                               size_t id_cell ) ;
      
   //-- Attributes

      // Generalities 
      size_t const DIM ;

      // Segment to side polyhedron intersector:
      GE_SegmentPolyhedron_INT* INTERSECTOR ;

      // Cell face iterators:
      PDE_CursorFEside* const sFE ;
      PDE_LocalFEbound* const bFE ;
} ;

#endif
