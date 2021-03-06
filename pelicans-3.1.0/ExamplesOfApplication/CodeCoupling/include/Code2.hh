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

#ifndef Code2_HH
#define Code2_HH

#include <PEL_Application.hh>

/*
  See Code2 documentation for information.
*/

class Code2 : public PEL_Application
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static Code2* create( PEL_Object* a_owner,
                            PEL_ModuleExplorer const* exp ) ;
      
   //-- Program core execution

      virtual void run( void ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~Code2( void ) ; 
      Code2( Code2 const& other ) ;
      Code2& operator=( Code2 const& other ) ;

      Code2( PEL_Object* a_owner, PEL_ModuleExplorer const* exp ) ;
      Code2( PEL_Object* a_owner, stringVector& args ) ;
      
   //-- Plug in
      
      Code2( void ) ;
      
      virtual Code2* create_replica( PEL_Object* a_owner,
				     PEL_ModuleExplorer const* exp ) const ;

   //-- Class attributes
      
      static Code2 const* PROTOTYPE ;

   //-- Attributes
       
} ;

#endif



