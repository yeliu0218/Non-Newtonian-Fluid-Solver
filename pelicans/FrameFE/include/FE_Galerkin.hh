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

#ifndef FE_GALERKIN_HH
#define FE_GALERKIN_HH

#include <FE_OneStepIteration.hh>

#include <intVector.hh>

class PEL_ObjectRegister ;
class PEL_Sequence ;
class PEL_Timer ;
class PEL_Vector ;

class GE_QuadratureRule ;
class GE_QRprovider ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LocalFE ;
class PDE_LocalFEcell ;
class PDE_LocalFEbound ;

/* 
Objects conceptually decomposing the steps performed in 
   `FE_StepByStepProgression::run' 
into substeps devoted to the solution of a system of partial differential 
equations discretized with a Galerkin finite element method. 
Servers of type `PDE_LocalFE::' provide the main tool involved in the 
associated computations.

The system of balance equations is written :

   material derivative of the unknowns = creation of the unknowns

Roughly :

   approximation of the material derivative = ( u - u_0 ) / dt

where dt is the current time step in `FE_StepByStepProgression::run'.
u and u_0 are computed at each integration point given by
a `PDE_LocalFEcell::' object. There are two possiblities :

1. u and u_0 simply represent two successive levels of storage of an 
   unknown field : the above approximation is a finite difference 
   Euler approximation of the time derivative of u.

2. u_0 represents the value of the unknown at the characteristic foot with
   respect to a given advection field : the above approximation represents
   the material derivative of u with respect to that advection field.
   Such a feature is provided by `FE_GalerkinCharacteristic::' objects.

Instances of `FE_Galerkin::' are mainly meant to be used in conjonction
with a `FE_GalerkinCharacteristic::' object.

FRAMEWORK INSTANTIATION

   CASE 1 : derivation of a concrete subclass

   1. Derive a concrete subclass, say MyGal.
   2. Choose a name for MyAppli, say "my_gal".
   3. Implement a private destructor.
   4. Declare all constructors private.
   5. Define the prototype to be registered :
      5.1 Implement a default constructor that initializes the
          `FE_Galerkin::' subobject by calling
               `FE_Galerkin( std::string const& )'
          with "my_gal" as argument.
          Example of pseudo-code :
          | MyGal:: MyGal( void ) : FE_Galerkin( "my_gal" ) {}
      5.2 Define and initialize a static instance by calling the default
          constructor.
             declaration (in the header file, eg MyGal.hh) :
             | static MyGal const* PROTOTYPE ;
             definition (in the implementation file, eg MyGal.cc) :
             | MyGal const* MyGal::PROTOTYPE = new MyGal() ;'
   6. Implement a private constructor that initializes the 
      `FE_Galerkin::' subobject by calling
            `FE_Galerkin( PEL_Object*, PDE_DomainAndFields const* )'.
      Example of pseudo-code :
      | MyGal:: MyGal( PEL_Object* a_owner,
      |                PDE_DomainAndFields const* dom,
      |                FE_SetOfParameters const* prms,
      |                PEL_ModuleExplorer const* exp )
      |    : FE_Galerkin( a_owner, dom, exp ), ...
      | { ... }
      Use
         `::add_one_convected_field'
      to register fields that occur through material derivatives.
   7. Implement the `::create_replica' method that allocates an object
      of type `MyGal' initialized using the private constructor described
      above, and subsequently return a pointer to that object.
      Example of pseudo-code :
      | MyGal* MyGal::create_replica( PEL_Object* a_owner,
      |                               PDE_DomainAndFields const* dom,
      |                               FE_SetOfParameters const* prms,
      |                               PEL_ModuleExplorer const* exp ) const
      | {
      |    PEL_LABEL( "MyGal::create_replica" ) ;
      |    PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;
      |    MyGal* result = new MyGal( a_owner, dom, prms, exp ) ;
      |    PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ;
      |    return result ;
      | }
   8. Implement the pure virtual methods, and possibly overrid the non
      pure virtual methods.


   CASE 2 : derivation of an abstract subclass

   1. Derive an abstract subclass, say MyGal.
   2. Implement a protected virtual destructor.
   3. Implement a protected constructor that initializes the
      `FE_Galerkin::' subobject by calling
               `FE_Galerkin( std::string const& )'
      Example of pseudo-code :
      | MyGal:: MyGal( std::string const& name ) 
      |    : FE_Galerkin( name ) {}
      This constructor is devoted to be used by the concrete subclasses 
      of MyGal for the registration of their prototype.
   4. Implement a protected constructor that initializes the 
      `PEL_Application::' subobject by calling
            `FE_Galerkin( PEL_Object*, PDE_DomainAndFields const* )'.
      Example of pseudo-code :
      | MyGal:: MyGal( PEL_Object* a_owner,
      |                PDE_DomainAndFields const* dom )
      |    : FE_Galerkin( a_owner, dom ), ...
      | { ... }
      This constructor is devoted to be used to initialize the MyGal
      base class subobject when creating objects of concrete subclasses
      of MyGal (such creations are performed in the `create_replica::'
      method whose implementation is deferred into those concrete subclasses).

PUBLISHED
*/

class PEL_EXPORT FE_Galerkin : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static FE_Galerkin* make( PEL_Object* a_owner,
                                PDE_DomainAndFields const* dom,
                                FE_SetOfParameters const* prms,
                                PEL_ModuleExplorer* exp ) ;

   //-- Features required by `FE_GalerkinCharacteristic::'

      // sequence of fields for which the material derivative discretization
      // is performed with a Galerkin Characteristic method
      PEL_Sequence const* convected_fields( void ) const ;

      // 1. storage level for which the values at IPs of the items of 
      //    `::convected_fields' are masked with the value at the 
      //    corresponding characteristic foot (computed with the same 
      //    storage level)
      // 2. storage level of the advection field used by 
      //    `FE_GalerkinCharacteristic::' to locate the characteristic foot
      size_t masked_level( void ) const ;

      // provider of quadrature rules used to perform the numerical
      // integration associated to material derivative terms
      virtual GE_QRprovider const* QRprovider_for_material_derivative( 
                                                             void ) const = 0 ;
   //-- Substeps of the step by step progression

      // Notify to `fe' of the calculation requirements that will be requested 
      // when calling `::build_cell_contribution_to_material_derivative'.
      virtual void transfer_calculation_requirements_for_material_derivative( 
                                        PDE_LocalFEcell* fe ) = 0 ;

      /*
      Perform the first step of the inner iteration conceptually
      decomposed as :
         1. reset_discrete_problem( t_it )

         2. for( fe->start() ; fe->is_valid() ; fe->go_next() )
              build_cell_contribution_to_material_derivative( t_it, fe )

         3. terminate_discrete_problem( t_it ).
      */
      virtual void reset_discrete_problem( FE_TimeIterator const* t_it ) = 0 ;


      // Compute and assemble the contribution of the material derivatives 
      // terms that may involve values at integration points that where
      // previouly masked, as defined by `PDE_LocalFE::mask_value_at_IP'
      // (this function must be called on behalf of `fe'
      // prior to `::build_cell_contribution_to_material_derivative').
      virtual void build_cell_contribution_to_material_derivative( 
                                              FE_TimeIterator const* t_it,
                                              PDE_LocalFEcell* fe ) = 0 ;

      /*
      Perform the first step of the inner iteration conceptually
      decomposed as :
         1. reset_discrete_problem( t_it )

         2. for( fe->start() ; fe->is_valid() ; fe->go_next() )
              build_cell_contribution_to_material_derivative( t_it, fe )

         3. terminate_discrete_problem( t_it ).
      */
      virtual void terminate_discrete_problem( 
                                      FE_TimeIterator const* t_it ) = 0 ;
      
   protected: //--------------------------------------------------------

      virtual ~FE_Galerkin( void ) ;

      FE_Galerkin( PEL_Object* a_owner, 
                   PDE_DomainAndFields const* dom,
                   PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      // for prototype registration only
      FE_Galerkin( std::string const& name ) ;

      virtual FE_Galerkin* create_replica( 
                                   PEL_Object* a_owner,
                                   PDE_DomainAndFields const* dom,
                                   FE_SetOfParameters const* prms,
                                   PEL_ModuleExplorer* exp ) const = 0 ;

   //-- Internals

      // Notify that the material derivative involving `ff' might be
      // discretized with a Characteristic Galerkin method. Successive calls
      // should always be done with the same `a_masked_level'.
      void add_one_convected_field( PDE_DiscreteField* ff,
                                    size_t a_masked_level ) ;

   //-- Preconditions, Postconditions, Invariant

      bool QRprovider_for_material_derivative_POST( 
                                         GE_QRprovider const* result ) const ;

      bool transfer_calculation_requirements_for_material_derivative_PRE( 
                                         PDE_LocalFEcell const* fe ) const ;

      bool reset_discrete_problem_PRE( FE_TimeIterator const* t_it ) const ;

      bool build_cell_contribution_to_material_derivative_PRE( 
  	                                    FE_TimeIterator const* t_it,
  	                                    PDE_LocalFEcell* fe ) const ;

      bool terminate_discrete_problem_PRE( 
                                       FE_TimeIterator const* t_it ) const ;

      virtual bool invariant( void ) const ;

   private: //----------------------------------------------------------

      FE_Galerkin( void ) ;
      FE_Galerkin( FE_Galerkin const& other ) ;
      FE_Galerkin& operator=( FE_Galerkin const& other ) ;

      static PEL_ObjectRegister* galerkin_map( void ) ;

   //-- Attributes

      PEL_Vector* CV_FIELDS ;
      size_t L_MASKED ;
} ;

#endif
