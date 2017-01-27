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

#ifndef FE_GALERKIN_CHARACTERISTIC_HH
#define FE_GALERKIN_CHARACTERISTIC_HH

#include <FE_OneStepIteration.hh>

#include <doubleVector.hh>

class PEL_List ;
class PEL_ListIterator ;
class PEL_ModuleExplorer ;
class PEL_Timer ;
class PEL_Vector ;

class PDE_CFootFinder ;
class PDE_DomainAndFields ;
class PDE_DiscreteField ;
class PDE_LocalFEcell ;

class FE_Galerkin ;

/*
Objects decomposing the steps performed in `FE_StepByStepProgression::run' 
into substeps devoted to the solution of a system of partial differential 
equations (PDEs) discretized with a Characteristic-Galerkin finite element
method.

Instances of `FE_Galerkin::' are handled, one for each equation of
the considered system of PDEs. Those instances are responsible for
building and solving a problem written symbolically :
   ( u - u_0 ) / dt = creation(u)

As for the `FE_GalerkinCharacteristic::' instance, it is responsible
making u_0 represent the value of a given storage level of a discrete
field at the characteristic foot with respect to a given advection field
(through the use of `PDE_LocalFEcell::mask_value_at_IP')
*/

class PEL_EXPORT FE_GalerkinCharacteristic : public FE_OneStepIteration
{

   public: //-----------------------------------------------------------

   //-- Substeps of the step by step progression

      // IMPLEMENTATION : Call `FE_Galerkin::do_before_time_stepping'
      // for all handled `FE_Galerkin::' instances.
      virtual void do_before_time_stepping( FE_TimeIterator const* t_it ) ;

      // IMPLEMENTATION : Call `FE_Galerkin::do_before_inner_iterations_stage'
      // for all handled `FE_Galerkin::' instances.
      virtual void do_before_inner_iterations_stage( 
                                           FE_TimeIterator const* t_it ) ;

      /*
      IMPLEMENTATION :

      1. Call `FE_Galerkin::reset_discrete_problem' for each handled
         `FE_Galerkin::' instance

      2. For each cell :

         2.1 Find the characteristic foot of each integration point.
         2.2 Mask the values corresponding to the convected fields with
             their values at those feet.
         2.3 Call `FE_Galerkin::build_cell_contribution_to_material_derivative'
             for all handled `FE_Galerkin::' instance.

      3. Call `FE_Galerkin::terminate_discrete_problem' for each handled 
         `FE_Galerkin::' instance,
      */
      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

      // IMPLEMENTATION : 
      // Call `FE_OneStepIteration::do_after_inner_iterations_stage'
      // for all handled `FE_Galerkin::' instances.
      virtual void do_after_inner_iterations_stage(
                                     FE_TimeIterator const* t_it ) ;

      // IMPLEMENTATION : 
      // Call `FE_OneStepIteration::do_after_time_stepping'
      // for all handled `FE_Galerkin::' instances.
      virtual void do_after_time_stepping( void ) ;

   //-- Elapsed times

      virtual void print_additional_times( std::ostream& os,
                                           size_t indent_width ) const ;

   //-- Savings for post-processing

      // IMPLEMENTATION : Call `FE_Galerkin::do_additional_savings'
      // for all handled `FE_Galerkin::' instances.
      virtual void do_additional_savings( FE_TimeIterator const* t_it,
                                          PDE_ResultSaver* rs ) ;

  //-- Persistence
      
      // IMPLEMENTATION : Call `FE_Galerkin::add_storable_objects'
      // for all handled `FE_Galerkin::' instances.
      virtual void add_storable_objects( PEL_ListIdentity* list ) ;      
      
  //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~FE_GalerkinCharacteristic( void ) ;
      FE_GalerkinCharacteristic( FE_GalerkinCharacteristic const& other ) ;
      FE_GalerkinCharacteristic& operator=( 
                                 FE_GalerkinCharacteristic const& other ) ;

      FE_GalerkinCharacteristic( PEL_Object* a_owner,
                                 PDE_DomainAndFields const* dom,
                                 FE_SetOfParameters const* prms,
                                 PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      FE_GalerkinCharacteristic( void ) ;

      virtual FE_GalerkinCharacteristic* create_replica( 
                                  PEL_Object* a_owner,
                                  PDE_DomainAndFields const* dom,
                                  FE_SetOfParameters const* prms,
                                  PEL_ModuleExplorer* exp ) const ;

   //-- Internals


      void insert_Galerkin_problem( FE_Galerkin* gpb ) ;
      
      void reset_statistics( void ) ;
      void update_statistics( PDE_CFootFinder const* finder ) ;
      void print_transport_statistics( std::ostream& os,
                                       size_t indent_width ) const ;

   //-- Class attributes

      static FE_GalerkinCharacteristic const* PROTOTYPE ;

   //-- Attributes
     
      PDE_DomainAndFields const* DOM ;

      size_t L_MASKED ;

      // List of FE_Galerkin solvers :
      PEL_List* GLIST ;
      PEL_ListIterator* G_IT ;

      // List of FE_Galerkin solvers sorted for material derivative
      // assembling :
      PEL_List* SORTED_GLIST ;
      PEL_ListIterator* S_G_IT ;
      
      PEL_Timer* TR_TIMER ;

      PEL_Vector* CV_FIELDS ;
      PDE_DiscreteField const* ADV ;
      size_t L_ADV ;

      PDE_LocalFEcell* cFE ;

      // CFLs :
      doubleVector CFL_MEAN ;
      doubleVector CFL_MIN ;
      doubleVector CFL_MAX ;

      // IPs :
      size_t NB_ITG_PTS ;
      double ITG_WEIGHT ;
      size_t NB_NO_FOUND_ITG_PTS ;
      double W_NO_FOUND ;
      size_t NB_INFLOWS ;
      double W_INFLOW ;
} ;


#endif
