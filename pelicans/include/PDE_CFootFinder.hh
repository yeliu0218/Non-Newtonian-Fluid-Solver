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

#ifndef PDE_C_FOOT_FINDER_HH
#define PDE_C_FOOT_FINDER_HH

#include <PEL_Object.hh>

#include <PDE_LocalFEcell.hh>

class PEL_ModuleExplorer ;

class GE_Color ;
class GE_Mpolyhedron ;
class GE_Point ;

class PDE_BoundFE ;
class PDE_CellFE ;
class PDE_DiscreteField ;
class PDE_GridFE ;
class PDE_MeshFE ;

/*
Servers which search for the foot of the characteristic issued from a given
point.

Those servers are deeply linked to the `PDE_LocalFEcell::' servers which
synchronize the search from its current mesh.

Searching the foot of the characteristic from `head' is finding a point `foot'
such as :     `head' = `foot' + `time_step'*`aa'
where - `aa' is the advective field of the problem ;
      - `time_step' the current time step.
The advective field `aa' my be expressed :
      - at `head' (backward-Euler method) ;
      - at `foot' (forward-Euler method) ;
      - ...

INSTANCE CREATION : `PDE_LocalFEcell::set_foot_finder'
INSTANCE DELIVERY : `PDE_LocalFEcell::foot_finder'
*/

class PEL_EXPORT PDE_CFootFinder : public PEL_Object
{
   public: //---------------------------------------------------------------

   //-- Search from characteristic head

      GE_Mpolyhedron const* polyhedron_of_head_cell( void ) const ;
      
      size_t mesh_id_of_head_cell( void ) const ;

      // Search the foot of the characteristic from `head'
      // (roughly:  `::foot' + `aa'*`time_step' = `head' ).
      void search_foot( PDE_DiscreteField const* aa,
                        size_t level,
                        double time_step,
                        GE_Point const* head ) ;
      
      bool search_has_been_performed( void ) const ;

   //-- Characteristic foot

      double value_at_foot( PDE_DiscreteField const* f,
                            size_t level,
                            size_t ic=0 ) const ;

      bool foot_has_been_found( void ) const ;
      
      bool foot_is_interior( void ) const ;
      
      bool foot_is_on_boundary( void ) const ;
      
      GE_Mpolyhedron const* polyhedron_of_foot_cell( void ) const ;
      
      GE_Color const* color_of_foot_cell( void ) const ;
      
      size_t mesh_id_of_foot_cell( void ) const ;
      
      virtual GE_Point const* foot( void ) const = 0 ;
      
      virtual GE_Point const* foot_ref( void ) const = 0 ;

      virtual double CFL_number( size_t dim ) const ;

   //-- Possible substitute if no foot has been found

      // Is a closest foot is investigate if the search has failed ?
      // IMPLEMENTATION : false
      virtual bool is_ersatz_foot_investigated( void ) const ;

      double value_at_ersatz_foot( PDE_DiscreteField const* f,
                                    size_t level,
                                    size_t ic=0 ) const ;      
      
      GE_Mpolyhedron const* polyhedron_of_ersatz_foot_cell( void ) const ;
      
      GE_Color const* color_of_ersatz_foot_cell( void ) const ;
      
      size_t mesh_id_of_ersatz_foot_cell( void ) const ;
      
      virtual GE_Point const* ersatz_foot( void ) const ;
      
      virtual GE_Point const* ersatz_foot_ref( void ) const ;
      
   protected: //------------------------------------------------------------

      virtual ~PDE_CFootFinder( void ) ;
      
      PDE_CFootFinder( PDE_LocalFEcell* fem,
                       PEL_ModuleExplorer const* exp ) ;

   //-- Perform search
      
      void raise_periodicity_error( void ) const ;

      virtual void search_foot( PDE_DiscreteField const* aa,
                                size_t level,
                                double time_step,
                                GE_Point const* head,
                                PDE_CellFE const* head_cell ) = 0 ;

   //-- Internal status

      void set_found_cell( PDE_CellFE const* a_cell ) ;

      void set_found_bound( PDE_BoundFE const* a_bound ) ;

      PDE_BoundFE const* foot_bound( void ) const ;

      PDE_CellFE const* foot_cell( void ) const ;
      
      PDE_MeshFE const* foot_mesh( void ) const ;

      virtual PDE_CellFE const* ersatz_foot_cell( void ) const ;
      
   //-- Preconditions, Postconditions, Invariant

      bool search_foot_PRE( PDE_DiscreteField const* aa,
                            size_t level,
                            double time_step,
                            GE_Point const* head,
                            PDE_CellFE const* head_cell ) const ;

      virtual bool search_foot_POST( PDE_DiscreteField const* aa,
                                     size_t level,
                                     double time_step,
                                     GE_Point const* head,
                                     PDE_CellFE const* head_cell ) const ;

      bool foot_PRE( void ) const ;

      bool foot_POST( GE_Point const* result ) const ;

      bool foot_ref_PRE( void ) const ;

      bool foot_ref_POST( GE_Point const* result ) const ;
      
      bool ersatz_foot_PRE( void ) const ;

      bool ersatz_foot_POST( GE_Point const* result ) const ;

      bool ersatz_foot_ref_PRE( void ) const ;

      bool ersatz_foot_ref_POST( GE_Point const* result ) const ;

      bool ersatz_foot_cell_PRE( void ) const ;

      bool ersatz_foot_cell_POST( PDE_CellFE const* result ) const ;
      
      virtual bool invariant( void ) const ;

   private: //--------------------------------------------------------------

      PDE_CFootFinder( void ) ;
      PDE_CFootFinder( PDE_CFootFinder const& other ) ;
      PDE_CFootFinder const& operator=( PDE_CFootFinder const& other ) ;
      
   //-- Instance delivery and initialization
      
      friend void PDE_LocalFEcell::set_foot_finder(
                                      PEL_ModuleExplorer const* exp ) ;
      
      static PDE_CFootFinder* create( PDE_LocalFEcell* fem,
                                      std::string const& concrete_name,
                                      PEL_ModuleExplorer const* exp,
                                      PDE_GridFE const* grid ) ;

      friend void PDE_LocalFEcell::synchronize_foot_finder( void ) ;
      
      void synchronize( PDE_LocalFEcell const* fem,
                        PDE_CellFE const* initial_cell ) ;

   //-- Attributes

      bool SEARCH_DONE ;
      PDE_CellFE const* HEAD_CELL ;
      GE_Point const* HEAD ;
      PDE_CellFE const* FOOT_CELL ;
      PDE_BoundFE const* FOOT_BOUND ;
} ;

#endif
