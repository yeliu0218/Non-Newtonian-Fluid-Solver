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

#ifndef PEL_TIC_READER_HH
#define PEL_TIC_READER_HH

#include <PEL_DataOnMeshingReader.hh>

class stringVector ;

/*
readers of savings performed by `::PEL_TICwriter' objects

PUBLISHED
*/

class PEL_EXPORT PEL_TICreader : public PEL_DataOnMeshingReader
{
   public: //-----------------------------------------------------------

   //-- Explorer on datas at current cycle

      virtual PEL_ModuleExplorer* meshing( void ) const ;

      virtual PEL_ModuleExplorer* fields( void ) const ;

      virtual PEL_ModuleExplorer* integration_domain( void ) const ;

      virtual PEL_ModuleExplorer* variables( void ) const ; 

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~PEL_TICreader(void ) ;
      PEL_TICreader( PEL_TICreader const& other ) ;
      PEL_TICreader& operator=( PEL_TICreader const& other ) ;

      PEL_TICreader( PEL_Object* a_owner, PEL_ModuleExplorer const* exp ) ;

      PEL_ModuleExplorer* restore_cycle( PEL_Module* m, size_t i_cycle ) ;
      PEL_ModuleExplorer* create_meshing( PEL_ModuleExplorer* m_exp ) ;
      PEL_ModuleExplorer* create_fields( PEL_ModuleExplorer* m_exp,
                                         size_t nb_vertices,
                                         size_t nb_meshes ) ;
      PEL_ModuleExplorer* create_idom( PEL_ModuleExplorer * m_exp ) ;

   //-- Plug in
      
      PEL_TICreader( void ) ;

      virtual PEL_TICreader* create_replica( 
                                     PEL_Object* a_owner,
                      		     PEL_ModuleExplorer const* exp ) const ;
      
   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;
      
   //-- Class attributes

      static PEL_TICreader const* PROTOTYPE ;

   //-- Attributes

      // Explorer at current cycle :
      PEL_ModuleExplorer* MESHING_EXP ;
      PEL_ModuleExplorer* FIELDS_EXP ;
      PEL_ModuleExplorer* I_DOM_EXP ;
      PEL_ModuleExplorer* VAR_EXP ;
} ;

#endif
