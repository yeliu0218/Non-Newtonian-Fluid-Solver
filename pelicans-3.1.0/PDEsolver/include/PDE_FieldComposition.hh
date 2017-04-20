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

#ifndef PDE_FIELD_COMPOSITION_HH
#define PDE_FIELD_COMPOSITION_HH

#include <PEL_Object.hh>

#include <string>

#include <doubleArray2D.hh>
#include <size_t_vector.hh>

class PDE_DiscreteField ;
class PDE_SetOfDiscreteFields ;
class PDE_SetOfFieldCompositions ;

class PEL_Iterator ;
class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;
class PEL_Vector ;

/*
Multicomponent valued functions of several variables where each variable
represents a PDE_DiscreteField object. 

Variables have as many components as the associated PDE_DiscreteField object

The dependance upon a PDE_DiscreteField object can be direct or through 
a composition with another PDE_FieldComposition object.

Creation of PDE_FieldComposition objects :
  `complete_internal_dependencies' should be called to avoid any difficulty
  due to  a possible composition with another PDE_FieldComposition object
  that was incompletely initialized when calling `add_one_composition'.

Use of PDE_FieldComposition objects :
  1. Iterate over all variables (either direct or implicit through 
     composition) and set their value.
  2. Call `compute'.
  3. Call `value'.

FRAMEWORK INSTANTIATION
  1. Derive a concrete subclass.
  2. Implement a static factory method and a private constructor.
     The PDE_FieldComposition subobject is initialized by calling
     `PDE_FieldComposition( PEL_Object*, std::string )'
     The variables are identified by calling the `add_one_variable' 
     or the `add_one_composition' method.
  3. Implement a destructor.
  4. Implement `nb_components()'
  5. Implement the computation of the value of the current instance
     from its specific dependance upon the variables in `compute_self()'.
     At this point, the required values of the PDE_DiscreteField object 
     variables are given by `variable_value()' and  the required values of 
     the PDE_FieldComposition object variables are given by calling the 
     `value()' method for these objects.
  6. Implement `value(size_t)', whose result has been computed in 
     `compute_self()'.

PUBLISHED
*/
   
class PEL_EXPORT PDE_FieldComposition : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_FieldComposition* make( PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp,
                                         size_t nb_sp_dimensions,
                                         PDE_SetOfDiscreteFields const* dfs ) ;

      // Perform the internal links with the `PDE_FieldComposition::'
      // objects that `self' depends on. These objects should be items of
      // `fcs' (a fatal error is raised if not).
      // IMPLEMENTATION : do nothing
      virtual void do_the_links( PDE_SetOfFieldCompositions const* fcs ) ;
      
      // Does all the variables of `self' and its compositions with internal
      // PDE_FieldComposition objects are also variables of `other' ?
      bool has_consistent_internal_dependencies( 
                                PDE_FieldComposition const* other ) const ;
      
      // Ensure that all variables of the compositions of `self' with other
      // PDE_FieldComposition objects are also variables of `self' .
      void complete_internal_dependencies( void ) ;

   //-- Description

      // name
      std::string const& name( void ) const ;

      // Concerning the result of calling `::value' with a given argument,
      // is it always the same ? 
      virtual bool is_constant( void ) const ;

      // number of components
      virtual size_t nb_components( void ) const = 0 ;

      // number of variables
      size_t nb_variables( void ) const ;

      // Is there a variable associated to `ff' ?
      bool is_a_variable( PDE_DiscreteField const* ff ) const ;

   //-- Iterator over variables (direct or implicit through composition)
      
      // Move iterator to first position.
      void start_variable_iterator( void ) ;

      // Move iterator one position.
      void go_next_variable( void ) ;

      // Is iterator position valid ?
      bool valid_variable( void ) const ;

      // the `PDE_DiscreteField::' object associated to current 
      // iterator position
      PDE_DiscreteField const* variable( void ) const ;

   //-- Computation and retrieval of self value

      // Has the `ic'-th component of the variable associated to `ff'
      // been set ?
      bool variable_value_is_set( PDE_DiscreteField const* ff,
                                  size_t ic ) const ;

      // Make `x' the value of the `ic'-th component of the variable
      // associated to `ff'.
      void set_variable_value( PDE_DiscreteField const* ff, 
                               size_t ic,
                               double x ) ;

      // value of `ic'-th component of the variable associated to `ff'
      double variable_value( PDE_DiscreteField const* ff,
                             size_t ic ) const ;
      /* 
      Compute the value of all components of self, in two steps :
      The algorithm is :
         1. Compute the value of all components of the PDE_FieldComposition
            objects composed with self ;
         2. Call `::compute_self'. */
      void compute( void ) ;

      // Has self been computed ?
      bool is_computed( void ) const ;
      
      // value of the `ic'-th component of self
      virtual double value( size_t ic ) const = 0 ;

   protected: //--------------------------------------------------------

   //-- Plug in

      virtual ~PDE_FieldComposition( void ) ;

      // Registration of the concrete subclass.
      PDE_FieldComposition( std::string const& a_concrete_name ) ;

      PDE_FieldComposition( PEL_Object* a_owner, std::string const& a_name ) ;
      
      virtual PDE_FieldComposition* create_replica(
                             PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp,
                             size_t nb_sp_dims,
                             PDE_SetOfDiscreteFields const* dfs ) const = 0 ;

      bool is_a_prototype( void ) const ;
      
   //-- Initialization

      // Add `fc' as a requirement of self.  
      void add_one_composition( PDE_FieldComposition* fc ) ;
      bool depends_on( PDE_FieldComposition const* fc ) const ;
      
      // Add `df' as a requirement of self.  
      void add_one_variable( PDE_DiscreteField const* df ) ;
      
   //-- Computation

      // Compute self.
      virtual void compute_self( void ) = 0 ;
      
   //-- Preconditions, Postconditions, Invariant

      virtual bool do_the_links_PRE(
                          PDE_SetOfFieldCompositions const* fcs ) const ;
      virtual bool do_the_links_POST(
                          PDE_SetOfFieldCompositions const* fcs ) const ;
      
      virtual bool value_PRE( size_t ic ) const ;

      virtual bool create_replica_PRE(
                             PEL_Object const* a_owner,
                             PEL_ModuleExplorer const* exp,
                             size_t nb_sp_dims,
                             PDE_SetOfDiscreteFields const* dfs ) const ;
      virtual bool create_replica_POST(
                             PDE_FieldComposition const* result,
                             PEL_Object const* a_owner,
                             PEL_ModuleExplorer const* exp,
                             size_t nb_sp_dims,
                             PDE_SetOfDiscreteFields const* dfs ) const ;

      virtual bool invariant( void ) const ;
      
   private: //----------------------------------------------------------
      
      PDE_FieldComposition( void ) ;
      PDE_FieldComposition( PDE_FieldComposition const& other ) ;
      PDE_FieldComposition& operator=( PDE_FieldComposition const& other ) ;
      
      void compute_from_provider( PDE_FieldComposition const* provider ) ;

      static PEL_ObjectRegister* plugins_map( void ) ;

   //-- Attributes

      std::string const NAME ;

      bool const PROTO ;
      
      // vector of PDE_DiscreteField*
      PEL_Vector* const VARS ;
      size_t NB_VARS ;
      size_t VAR_IDX ;
      size_t_vector GLOB_2_iVAR ;
      doubleArray2D VAR_VALUES ;

      // vector of PDE_FieldComposition*      
      PEL_Vector* const LAWS ; 

      bool COMPUTED ;
      
      PDE_FieldComposition const* VALUES_PROVIDER ;
      
} ;

#endif
