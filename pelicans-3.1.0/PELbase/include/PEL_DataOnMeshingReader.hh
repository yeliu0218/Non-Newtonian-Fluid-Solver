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

#ifndef PEL_DATA_ON_MESHING_READER_HH
#define PEL_DATA_ON_MESHING_READER_HH

#include <PEL_Object.hh>

class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;

/*
readers of savings performed by `::PEL_DataOnMeshingWriter' objects

FRAMEWORK INSTANTIATION
   1. Derive a concrete subclass.
   2. Create a static class pointer which defines the model of the concrete
      class (prototype). It is built calling the prototype constructor :
          `::PEL_DataOnMeshingReader( std::string const& )'
   3. Implement the function :
          `::create_replica( PEL_Object*, PEL_ModuleExplorer const* ) const'
      which create a new instance of the concrete class, calling the 
      constructor :
        `::PEL_DataOnMeshingReader( PEL_Object* )'
   4. Implement a destructor
   5. Implement all the virtual functions

PUBLISHED
*/

class PEL_EXPORT PEL_DataOnMeshingReader : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PEL_DataOnMeshingReader* make(
                                 PEL_Object* a_owner,
                                 std::string const& name,
                                 PEL_ModuleExplorer const* exp ) ;

   //-- Explorer on datas at current cycle

      virtual PEL_ModuleExplorer* meshing( void ) const = 0 ;

      virtual PEL_ModuleExplorer* fields( void ) const = 0 ;

      virtual PEL_ModuleExplorer* integration_domain( void ) const = 0 ;

      virtual PEL_ModuleExplorer* variables( void ) const = 0 ;
      
   protected: //--------------------------------------------------------

   //-- Plug in

      virtual ~PEL_DataOnMeshingReader( void ) ;

      PEL_DataOnMeshingReader( std::string const& name ) ;

      PEL_DataOnMeshingReader( PEL_Object* a_owner ) ;

      virtual PEL_DataOnMeshingReader* create_replica( 
                                  PEL_Object* a_owner,
				  PEL_ModuleExplorer const* exp ) const = 0 ;

      bool is_a_prototype( void ) const ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool meshing_POST( PEL_ModuleExplorer const* result ) const ;

      virtual bool fields_POST( PEL_ModuleExplorer const* result ) const ;

      virtual bool integration_domain_POST(
                                 PEL_ModuleExplorer const* result ) const ;

      virtual bool variables_POST( PEL_ModuleExplorer const* result ) const ;

      virtual bool create_replica_PRE( PEL_Object const* a_owner,
                                       PEL_ModuleExplorer const* exp ) const ;
      virtual bool create_replica_POST( PEL_DataOnMeshingReader const* result,
                                        PEL_Object const* a_owner,
                                        PEL_ModuleExplorer const* exp ) const ;

      virtual bool invariant( void ) const ;
      
   private: //----------------------------------------------------------

      PEL_DataOnMeshingReader( void ) ;
      PEL_DataOnMeshingReader( PEL_DataOnMeshingReader const& other ) ;
      PEL_DataOnMeshingReader& operator=( 
                               PEL_DataOnMeshingReader const& other ) ;

      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Attributes

      bool IS_PROTO ;
};

#endif
