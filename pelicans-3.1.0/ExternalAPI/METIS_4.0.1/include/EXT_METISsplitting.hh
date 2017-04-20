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

#ifndef EXT_METIS_SPLITTING_HH
#define EXT_METIS_SPLITTING_HH

#include <GE_SplittingStrategy.hh>

#include <intVector.hh>
class PEL_Timer ;

/*
Meshing splitting strategies using the METIS utility.

METIS is a family of programs for partitioning unstructured graphs and
hypergraphs and computing fill-reducing orderings of sparse matrices.
The underlying algorithms used by METIS are based on the state-of-the-art
multilevel paradigm that has been shown to produce high quality results
and scale to very large problems.
  
http://glaros.dtc.umn.edu/gkhome/views/metis

PUBLISHED  
*/

class EXT_METISsplitting : public GE_SplittingStrategy
{
   public: //------------------------------------------------------------

   //-- Cell balance

      virtual size_t cell_rank( size_t mesh_id ) const ;
            
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //---------------------------------------------------------

   private: //-----------------------------------------------------------

      EXT_METISsplitting( void ) ;
     ~EXT_METISsplitting( void ) ;
      EXT_METISsplitting( EXT_METISsplitting const& other ) ;
      EXT_METISsplitting& operator=( EXT_METISsplitting const& other ) ;

      EXT_METISsplitting( PEL_Object* a_owner,
			  PEL_ModuleExplorer const* exp,
			  GE_Meshing* meshing,
			  PEL_Communicator const* com  ) ;
      
      EXT_METISsplitting( PEL_Object* a_owner,
			  PEL_ModuleExplorer const* exp,
			  GE_Meshing* meshing,
			  size_t nb_rks, size_t rk ) ;

   //-- Plug in
      
      virtual EXT_METISsplitting* create_replica(
                          PEL_Object* a_owner,
                          PEL_ModuleExplorer const* exp,
                          GE_Meshing* meshing,
                          PEL_Communicator const* com ) const ;
      
      virtual EXT_METISsplitting* create_replica(
                          PEL_Object* a_owner,
                          PEL_ModuleExplorer const* exp,
                          GE_Meshing* meshing,
                          size_t nb_rks, size_t rk ) const ;

   //-- Cell balance

      void METIS_balancing( GE_Meshing* meshing ) ;

      int METIS_cell_type( std::string const& cell_poly_name ) const ;
 
   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Class attributes

      static EXT_METISsplitting const* PROTOTYPE ;

   //-- Attributes

      intVector CELL_RANK ;

      // Verbose:
      PEL_Timer* TIMER ;
      intVector RANK_CELLS ;
} ;

#endif
