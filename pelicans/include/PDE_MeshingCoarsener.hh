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

#ifndef PDE_MESHING_COARSENER_HH
#define PDE_MESHING_COARSENER_HH

#include <PEL_Object.hh>

#include <PDE_AdapterCHARMS.hh>

class PEL_Vector ;
class PDE_Activator ;

/*
HIGHLY UNSTABLE CLASS
*/

class PEL_EXPORT PDE_MeshingCoarsener : public PEL_Object
{
   public: //-----------------------------------------------------------
      
      void reset( void ) ;
      
      void prepare_for_coarsening( void ) ;
      
      size_t nb_levels( void ) const ;
      
      size_t current_fine_level( void ) const ;
      
      void do_one_coarsening( void ) ;

      void do_one_uncoarsening( void ) ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_MeshingCoarsener( void ) ;
     ~PDE_MeshingCoarsener( void ) ;
      PDE_MeshingCoarsener( PDE_MeshingCoarsener const& other ) ;
      PDE_MeshingCoarsener& operator=( 
                            PDE_MeshingCoarsener const& other ) ;

      friend PDE_MeshingCoarsener* 
             PDE_AdapterCHARMS:: meshing_coarsener( void ) const ;

      PDE_MeshingCoarsener( PEL_Object* a_owner, 
                            PDE_DomainAndFields const* dom,
                            PDE_GridFE* a_grid,
                            PDE_Activator* activator,
                            size_t verbose_level ) ;

    //-- Attributes

      PDE_DomainAndFields const* DOM ;
      PDE_GridFE* GRID ;
      PDE_Activator* ACTIVATOR ;
      PEL_Vector* BF_LISTS ;
      size_t LEVEL_MAX ;
      size_t FINE_LEVEL ;
      size_t VERB_LEVEL ;
} ;

#endif
