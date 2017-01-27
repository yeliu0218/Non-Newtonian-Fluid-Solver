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

#ifndef PEL_ROOT_HH
#define PEL_ROOT_HH

#include <PEL_Object.hh>

#include <string>

class PEL_Vector ;

/*
Object devoted to be the upper-most owner in the ownership method
of lifetime management.

Implemented as a singleton.

The programm execution consists of five stages :
   1. Initial stage (Big-Bang time) : all statics are initialized and
      the only instance of PEL_Root is created.
   2. The data deck is read and stored in memory.
   3. An instance of a concrete subclass of PEL_Application is created.
   4. Program core execution : the program execution proceeds by performing
      its specific tasks. In particular, objects are created and organized
      into ownership trees whose root node is either the unique instance
      of PEL_Root or the NULL object.
   5. Final stage : termination of the only instance of PEL_Root, leading to 
      the termination of all objects belonging to a ownership tree whose
      root node is not the NULL object.
   The only instance of PEL_Root is mainly involved in the initial and the
   final stages, and at any time an object is created to exist until the end
   of the program execution
*/

class PEL_EXPORT PEL_Root : public PEL_Object
{
   public: //--------------------------------------------------------------

   //-- Instance delivery and initialization
      
      // the only instance 
      static PEL_Object* object( size_t i=0 ) ;

   //-- Termination
      
      // Terminate self: if i<j, object(i) is destroyed before object(j).
      static void cleanup( void ) ;

   protected: //-----------------------------------------------------------

   private: //-------------------------------------------------------------
  
      PEL_Root( void ) ;
     ~PEL_Root( void ) ;
      PEL_Root( PEL_Root const& other ) ;
      PEL_Root& operator=( PEL_Root const& other ) ;

      PEL_Root( PEL_Object* a_owner ) ;

      static PEL_Vector* objects( bool destroy = false ) ;

} ;

#endif
