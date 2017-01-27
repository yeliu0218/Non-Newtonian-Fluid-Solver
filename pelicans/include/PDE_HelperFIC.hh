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

#ifndef PDE_HELPER_FIC_HH
#define PDE_HELPER_FIC_HH

#include <PEL_Object.hh>

#include <PDE_CursorFEside.hh>

#include <vector>

class PDE_CellFE ;
class PDE_FaceFE ;
class PDE_MeshFE ;

class PEL_EXPORT PDE_HelperFIC : public PEL_Object
{
   public: //-----------------------------------------------------------
      
   //-- Status
      
      size_t nb_space_dimensions( void ) const ;
      
      size_t current_side_id( void ) const ;
      
   //-- Coarser level
      
      size_t parent_cell_id( size_t i_adj ) const ;
      
   //-- Interpolation
      
      void prepare_for_interpolation( void ) ;
      
      bool ready_for_interpolation( void ) const ;
      
      size_t node( size_t i_pt, PDE_DiscreteField const* ff ) const ;
      
      GE_Point const* point_of_node( size_t i_pt ) const ;
            
      void set_calculation_point( GE_Point const* pt ) ;
      
      double interpolated_value( doubleVector const& value_at_nodes ) const ;
      
   protected: //--------------------------------------------------------
      
   private: //----------------------------------------------------------
      
      PDE_HelperFIC( void ) ;
     ~PDE_HelperFIC( void ) ;
      PDE_HelperFIC( PDE_HelperFIC const& other ) ;
      PDE_HelperFIC& operator=( PDE_HelperFIC const& other ) ;
      
   //-- Selective export for PDE_CursorFEside
      
      PDE_HelperFIC( PEL_Object* a_owner ) ;
      
      // to simulate selective export
      friend 
      PDE_HelperFIC* PDE_CursorFEside:: helper_FIC( void ) const ;
      
   //-- Internals
      
      void set_side( PDE_FaceFE const* a_side, 
                     size_t id_cell_0,
                     size_t id_cell_1 ) ;
            
   //-- Attributes
      
      size_t_vector ID_NUMBER ;
      size_t NB_DIMS ;
      PDE_FaceFE const* SIDE ;
      std::vector< PDE_MeshFE* > MESHES ;
      double X_REF ;
      double Y_REF ;
} ;

#endif
