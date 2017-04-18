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

#ifndef PDE_ACTIVATOR_HH
#define PDE_ACTIVATOR_HH

#include <PEL_Object.hh>

#include <set>

class PDE_AdaptationRequest ;
class PDE_BasisFunctionCell ;
class PDE_CellFE ;
class PDE_DiscreteField ;

class PEL_EXPORT PDE_Activator : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization
      
      static PDE_Activator* create( PEL_Object* a_owner, 
                                    bool hierarchical_basis,
                                    size_t verbose_level ) ;
      
   //-- Configuration

      void add_excluded_field( PDE_DiscreteField const* ff ) ;

      bool is_excluded( PDE_DiscreteField const* ff ) const ;

   //-- Adaptation
      
      void split_meshes( PDE_AdaptationRequest* adap ) ;
            
      void refine( PDE_AdaptationRequest* adap,
                   bool make_refined_DOFs_bad_double ) ;
      
      void unrefine( PDE_AdaptationRequest* adap ) ;
      
      void unsplit_meshes_1( PDE_AdaptationRequest* adap ) ;
      
   //-- Utilities

      void activate_and_deactivate_child_cells( PDE_CellFE* ccell ) ;
      
      void deactivate_and_activate_child_cells( PDE_CellFE* ccell ) ;

   //-- Statistics
      
      void reset_counters( void ) ;
      
      size_t nb_activated_cells( void ) const ;
      
      size_t nb_deactivated_cells( void ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_Activator( void ) ;
     ~PDE_Activator( void ) ;
      PDE_Activator( PDE_Activator const& other ) ;
      PDE_Activator& operator=( PDE_Activator const& other ) ;

      PDE_Activator( PEL_Object* a_owner, 
                     bool hierarchical_basis,
                     size_t verbose_level )  ;

      void quasi_hierarchical_refine( 
                              PDE_AdaptationRequest* adap,
                              bool make_refined_DOFs_bad_double ) const ;

      void hierarchical_refine( 
                              PDE_AdaptationRequest* adap,
                              bool make_refined_DOFs_bad_double ) const ;

      void quasi_hierarchical_unrefine( PDE_AdaptationRequest* adap ) ;

      void hierarchical_unrefine( PDE_AdaptationRequest* adap ) ;

      static void print_bf( std::ostream& os, 
                            PDE_BasisFunctionCell* bf ) ;

    //-- Attributes

      bool HB ;
      bool SOMETHING_CHANGED ;
      size_t RLEVEL ;
      size_t VERB_LEVEL ;

      std::set< PDE_DiscreteField const* > EXCLUDED_FIELDS ;

      size_t NB_ACT_CELLS ;
      size_t NB_DEACT_CELLS ;
} ;

#endif
