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

#ifndef EXT_SILO_WRITER_HH
#define EXT_SILO_WRITER_HH

#include <PEL_DataOnMeshingWriter.hh>

/*
writers for the data postprocessor MESHTV.
  
PUBLISHED
*/

class EXT_SiloWriter : public PEL_DataOnMeshingWriter
{
   public: //-----------------------------------------------------------

   //-- Write

      virtual void write_cycle( PEL_ModuleExplorer const* exp ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~EXT_SiloWriter(void ) ;
      EXT_SiloWriter( EXT_SiloWriter const& other ) ;
      EXT_SiloWriter& operator=( EXT_SiloWriter const& other ) ;
      
      EXT_SiloWriter( PEL_Object* a_owner,
		      PEL_ModuleExplorer const* exp  ) ;

   //-- Plug in

      EXT_SiloWriter( void ) ;

      virtual EXT_SiloWriter* create_replica( 
                                       PEL_Object* a_owner,
				       PEL_ModuleExplorer const* exp ) const ;

   //-- Write

      void write_grid( PEL_ModuleExplorer const* exp ) ;

      void write_field( int cycle_number, 
                        PEL_ModuleExplorer const* exp ) ;

   //-- Class attributes

      static EXT_SiloWriter const* PROTOTYPE ;

   //-- Attributes

      struct DBfile* dbfile ;
      std::string dir_name ;
      std::string meshing_name ;
} ;

#endif
