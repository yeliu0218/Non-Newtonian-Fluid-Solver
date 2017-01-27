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

#ifndef GE_COORDINATE_SPLITTING_HH
#define GE_COORDINATE_SPLITTING_HH

#include <GE_SplittingStrategy.hh>

#include <intVector.hh>
class PEL_Data ;

/*
Meshing splitting strategies based on formulas that use the coordinate of the 
cell centers to determine the cells that will be assigned to the various
processes.

Example:
  
   MODULE splitting_strategy
       concrete_name = "GE_CoordinateSplitting"
       $DS_X = component( $DV_X, 0 )
       $DS_Y = component( $DV_X, 1 )
       coordinate_splitting_formula = 
            ( $DS_X < 1.0 && $DS_Y < 1.0 ? 2 :
              $DS_X < 1.0 && $DS_Y > 1.0 ? 0 :
              $DS_X > 1.0 && $DS_Y > 1.0 ? 3 :
              1 )
   END MODULE splitting_strategy

Remark:

   see `PEL_GroupExp::' for useful expressions for coordinate_splitting_formula
     
PUBLISHED
*/

class PEL_EXPORT GE_CoordinateSplitting : public GE_SplittingStrategy
{
   public: //------------------------------------------------------------

   //-- Cell balance
      
      virtual size_t cell_rank( size_t mesh_id ) const ;
      
   protected: //---------------------------------------------------------

   private: //-----------------------------------------------------------

      GE_CoordinateSplitting( void ) ;
     ~GE_CoordinateSplitting( void ) ;
      GE_CoordinateSplitting( GE_CoordinateSplitting const& other ) ;
      GE_CoordinateSplitting& operator=( 
                              GE_CoordinateSplitting const& other ) ;
      
      GE_CoordinateSplitting( PEL_Object* a_owner,
			      PEL_ModuleExplorer const* exp,
			      GE_Meshing* meshing,
                              PEL_Communicator const* com ) ;
      
      GE_CoordinateSplitting( PEL_Object* a_owner,
			      PEL_ModuleExplorer const* exp,
			      GE_Meshing* meshing,
                              size_t nb_rks, size_t rk ) ;

   //-- Plug in
      
      virtual GE_CoordinateSplitting* create_replica(
                              PEL_Object* a_owner,
                              PEL_ModuleExplorer const* exp,
                              GE_Meshing* meshing,
                              PEL_Communicator const* com ) const ;

      virtual GE_CoordinateSplitting* create_replica(
                              PEL_Object* a_owner,
                              PEL_ModuleExplorer const* exp,
                              GE_Meshing* meshing,
                              size_t nb_rks, size_t rk ) const ;
      
   //-- Cell balance
      
      void search_owners_from_coords( GE_Meshing* meshing,
                                      PEL_ModuleExplorer const* exp ) ;
      
      size_t point_owner( PEL_Data const* formula ) const ;
 
   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Class attributes

      static GE_CoordinateSplitting const* PROTOTYPE ;

   //-- Attributes
      
      intVector CELL_RANK ;
} ;

#endif
