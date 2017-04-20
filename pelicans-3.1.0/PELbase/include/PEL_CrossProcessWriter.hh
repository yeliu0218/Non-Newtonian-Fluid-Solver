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

#ifndef PEL_CROSS_PROCESS_WRITER_HH
#define PEL_CROSS_PROCESS_WRITER_HH

#include <PEL_DataOnMeshingWriter.hh>
#include <doubleArray2D.hh>
#include <size_t_vector.hh>
class PEL_List ;
class PEL_Module ;
class PEL_Communicator ;
class PEL_DoubleComparator ;
class intArray2D ;
class doubleArray2D ;

/*

  Writer which collects data on several meshings to build entire domain saving.

  PUBLISHED
*/

class PEL_EXPORT PEL_CrossProcessWriter : public PEL_DataOnMeshingWriter
{
   public: //-----------------------------------------------------------

   //-- Write

      virtual void write_cycle( PEL_ModuleExplorer const* exp ) ;
      
      // IMPLEMENTATION : true
      virtual bool is_parallel_writer( void ) const ;

   protected: //--------------------------------------------------------
      
   private: //----------------------------------------------------------

     ~PEL_CrossProcessWriter(void ) ;
      PEL_CrossProcessWriter( PEL_CrossProcessWriter const& other ) ;
      PEL_CrossProcessWriter& operator=( 
                              PEL_CrossProcessWriter const& other ) ;

      PEL_CrossProcessWriter( PEL_Object* a_owner, 
                              PEL_ModuleExplorer const* exp ) ;
      
   //-- Plug in
      
      PEL_CrossProcessWriter( void ) ;

      virtual PEL_CrossProcessWriter* create_replica( 
                                     PEL_Object* a_owner,
                      		     PEL_ModuleExplorer const* exp ) const ;

   //-- Internals

      PEL_Module* collect_variables( PEL_Object* a_owner,
                                     PEL_ModuleExplorer* exp,
                                     size_t& nb ) ;
      PEL_Module* collect_meshing( PEL_Object* a_owner,
                                   PEL_ModuleExplorer* exp ) ;
      PEL_Module* collect_field( PEL_Object* a_owner,
                                 PEL_ModuleExplorer* exp,
                                 size_t& nbf  ) ;
      PEL_Module* collect_data( PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) ;     
      
      void compute_global_vector_of_indices( 
                                 intVector& glob_vect_indices,
                                 intVector const& loc_vect_indices,
                                 size_t_vector const& loc2glob,
                                 int index_to_ignore ) const ;
      
      void compute_global_cell2face( 
                                 intVector& glob_cell_nb_faces,
                                 intArray2D& glob_cell2face,
                                 intVector const& loc_cell_nb_faces,
                                 intArray2D const& loc_cell2face,
                                 size_t_vector const& loc2glob_cells ) const ;
      
      void compute_global_face2cell( 
                                 intArray2D& glob_face2cell,
                                 intArray2D const& loc_face2cell,
                                 size_t_vector const& loc2glob_cells,
                                 size_t_vector const& loc2glob_faces ) const ;
      
   //-- Class attributes

      static PEL_CrossProcessWriter const* PROTOTYPE ;

   //-- Attributes
      
      PEL_List* const SUB_WRITERS ;
      PEL_Communicator const* const COM ;
      PEL_DoubleComparator const* const COORDS_CMP ;
      double const DBL_EPSI ;
      double const DBL_MINI ;
      size_t_vector loc2glob_VERTS ;
      size_t_vector loc2glob_FACES ;
      size_t_vector loc2glob_CELLS ;
      size_t NB_VERTS ;
      size_t NB_FACES ;
      size_t NB_CELLS ;
      int HALO_COLOR_INDEX ;
      
} ;

#endif
