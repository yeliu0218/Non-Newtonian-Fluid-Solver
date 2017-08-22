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

#ifndef AP_TUTORIAL_G_HH
#define AP_TUTORIAL_G_HH

#include <PEL_Application.hh>

class LA_Matrix ;
class LA_SeqVector ;
class LA_Solver ;
class LA_Vector ;

class PDE_DiscreteField ;
class PDE_LocalEquation ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_ResultSaver ;
class PDE_SetOfBCs ;
class PDE_SystemNumbering ;

/*
PUBLISHED
*/

class AP_TutorialG : public PEL_Application
{
   public: //-----------------------------------------------------------

   //-- Program core execution

      virtual void run( void ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~AP_TutorialG( void ) ;
      AP_TutorialG( AP_TutorialG const& other ) ;
      AP_TutorialG& operator=( AP_TutorialG const& other ) ;

      AP_TutorialG( PEL_Object* a_owner, 
                    PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      AP_TutorialG( void ) ;

      virtual AP_TutorialG* create_replica( 
                             PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp ) const ;

   //-- Discrete problem building

      void loop_on_cells( void ) const ;

      void loop_on_bounds( void ) const ;

   //-- Class attributes

      static AP_TutorialG const* PROTOTYPE ;
      
   //-- Attributes

      PDE_DiscreteField* TT ;
      double CONDUCTIVITY ;
      PDE_SetOfBCs const* BCs ;

      PDE_LocalEquation* ELEMENT_EQ ;
      PDE_LocalFEcell* cFE ;
      PDE_LocalFEbound* bFE ;

      PDE_ResultSaver* SAVER ;

      PDE_SystemNumbering* NMB ;
      LA_Matrix* A ;
      LA_Vector* F ;
      LA_Vector* X ;
      LA_SeqVector* X_LOC ;
      
      LA_Solver* SOLVER ;
} ;

#endif
