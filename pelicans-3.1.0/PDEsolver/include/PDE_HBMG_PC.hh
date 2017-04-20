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

#ifndef PDE_HBMG_PC_HH
#define PDE_HBMG_PC_HH

#include <PDE_GeometricMultilevel_PC.hh>

#include <vector>

class LA_Solver ;
class LA_SeqMatrix ;
class LA_SeqVector ;

class PEL_EXPORT PDE_HBMG_PC : public PDE_GeometricMultilevel_PC
{
   public: //----------------------------------------------------------

   //-- Instance delivery and initialization
      
      virtual PDE_HBMG_PC* create_clone( PEL_Object* a_owner ) const ;

   //-- Linear system solution
      
      virtual void solve( LA_Vector const* rhs,
                          LA_Vector* sol ) ;

      virtual bool successful_solve( void ) const ;

      virtual size_t nb_cycles_performed( void ) const ;

   protected: //-------------------------------------------------------

   private: //---------------------------------------------------------

     ~PDE_HBMG_PC( void ) ; 
      PDE_HBMG_PC( PDE_HBMG_PC const& other ) ;
      PDE_HBMG_PC& operator=( PDE_HBMG_PC const& other ) ;

      PDE_HBMG_PC( PEL_Object* a_owner,
                   PEL_ModuleExplorer const* exp ) ;

      PDE_HBMG_PC( PEL_Object* a_owner,
                   PDE_HBMG_PC const* other ) ;

   //-- Plug in

      PDE_HBMG_PC( void ) ;

      virtual PDE_HBMG_PC* create_replica( 
                                    PEL_Object* a_owner,
                                    PEL_ModuleExplorer const* exp ) const ;

   //-- Internals

      void do_one_multigrid_iteration( size_t level,
                                       LA_SeqMatrix const* mat,
                                       LA_SeqVector const* rhs,
                                       LA_SeqVector* sol ) const ;

   //-- Class attributes

      static PDE_HBMG_PC const* PROTOTYPE ;

   //-- Attributes

      double TOL ;
      bool CONVERGED ;
      size_t NB_CYCLES ;
      size_t NB_PRESMOO ;
      size_t NB_POSTSMOO ;
      size_t MC ;
      size_t I_CYCLE ;

      LA_Solver* SOLVER ;

      bool SOLVE_OK ;
      size_t VERBOSE ;
      std::string INDENT ;
} ;

#endif

        
