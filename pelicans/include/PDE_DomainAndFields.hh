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

#ifndef PDE_DOMAIN_AND_FIELDS_HH
#define PDE_DOMAIN_AND_FIELDS_HH

#include <PEL_Object.hh>

class PEL_Communicator ;
class PEL_Module ;
class PEL_ModuleExplorer ;
class PEL_Vector ;

class PDE_AdapterCHARMS ;
class PDE_AdapterHN ;
class PDE_CrossProcessBFNumbering ;
class PDE_CrossProcessNodeNumbering ;
class PDE_DomainBuilder ;
class PDE_IDomain ;
class PDE_ResultSaver ;
class PDE_SetOfBCs ;
class PDE_SetOfDiscreteFields ;
class PDE_SetOfFieldCompositions ;

class PDE_CursorFEside ;
class PDE_LocalFEcell ;
class PDE_LocalFEbound ;
class PDE_LocalFEinterface ;

class GE_Color ;
class GE_SetOfPoints ;

/*
Servers, associated to a Module Hierarchy, that deliver fundamental objects
required to approximate a system of PDEs (Partial Differential Equations)
involving fields defined on a geometrical domain triangulated with a
meshing.

The Module Hierarchy includes the specifications of the
delivered objects.

A meshing of a geometric domain is defined as a collection of finitely many
cells such that :
   - the union of all cells is the closure of the domain ;
   - two distinct cells have disjoint interiors ;
   - any side of any cell is either a side of another cell or a subset
     of the domain boundary.
*/

class PEL_EXPORT PDE_DomainAndFields : public PEL_Object
{
   public: //------------------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_DomainAndFields* create( PEL_Object* a_owner,
                                          PEL_ModuleExplorer const* exp ) ;

      // Create and return a new instance of the class and perform all
      // necessary tasks to reinterpret `self' as a subdomain of
      // a logical global domain in a processing distributed over
      // all processes in `com'.
      static PDE_DomainAndFields* create( PEL_Object* a_owner,
                                          PEL_ModuleExplorer const* exp,
                                          PEL_Communicator const* com ) ;

   //-- Characteristics

      // name, uniquely determining `self'
      std::string const& name( void ) const ;

   //-- Modifiers(5)

      // Create a new discrete field called `name_of_new_field', identical
      // to the already existing discrete field called `model_name',
      // and add it to `::set_of_discrete_fields()'.
      void duplicate_field( std::string const& model_name,
                            std::string const& name_of_new_field ) const ;

      // Create a new discrete field called `name_of_new_field', with the
      // same discretization than the already existing discrete field called
      // `model_name', and add it to `::set_of_discrete_fields()'.
      // The values of all the DOFs of the newly created discrete field
      // (of name `name_of_new_field') are set to zero.
      void duplicate_field( std::string const& model_name,
                            std::string const& name_of_new_field,
                            size_t nb_components_of_new_field,
                            size_t storage_depth_of_new_field ) const ;

      // Create new discrete field according to the Module Hierarchy
      // attached to `exp', and add them to `::set_of_discrete_fields()'.
      void append_fields( PEL_ModuleExplorer const* exp ) ;

      void apply_requests_of_DOFs_values_modules( bool verbose ) const ;

      void apply_requests_of_DOFs_imposed_value_modules( bool verbose ) const ;

   //-- Conformal Adjacent Domains(9.9)

      size_t nb_conformal_adjacent_domains( void ) const ;

      PDE_DomainAndFields const* conformal_adjacent_domain( size_t i ) const ;

   //-- Distributed processing(10)

      PDE_CrossProcessNodeNumbering* create_CrossProcessNodeNumbering(
                                                   PEL_Object* a_owner ) const ;

      PDE_CrossProcessBFNumbering* create_CrossProcessBFNumbering(
                                                   PEL_Object* a_owner ) const ;

      // Does `self' represent a subdomain of a logical global domain
      // involved in a processing distributed over several processes ?
      bool is_distributed( void ) const ;

   //-- Problem Definition(100)

      // number of space dimensions
      size_t nb_space_dimensions( void ) const ;

      // the `PDE_SetOfDiscreteFields::' object that conforms with the
      // associated Module Hierarchy
      PDE_SetOfDiscreteFields* set_of_discrete_fields( void ) const ;

      // Is the field called `fname' defined on the domain boundary ?
      bool is_defined_on_boundary( std::string const& fname ) const ;

      // the `PDE_SetOfFieldCompositions::' object that conforms with the
      // associated Module Hierarchy
      PDE_SetOfFieldCompositions const* set_of_field_compositions( void ) const ;

      // the `PDE_SetOfBCs::' object that conforms with the
      // associated Module Hierarchy
      PDE_SetOfBCs const* set_of_boundary_conditions( void ) const ;

      // the set of mesh vertices of the grid generated according to the
      // associated Module Hierarchy
      GE_SetOfPoints const* set_of_vertices( void ) const ;

      // explorer on the (sub)Module called "decohesion" in the
      // associated Module Hierarchy (if any)
      PEL_ModuleExplorer const* decohesion_explorer( void ) const ;

   //-- Discretization on a Finite Element Grid(120)

      // Create and return a `PDE_LocalFEcell::' object that conforms with the
      // associated Module Hierarchy.
      PDE_LocalFEcell* create_LocalFEcell( PEL_Object* a_owner ) const ;

      // Create and return a `PDE_CursorFEside::' object that conforms with
      // the associated Module Hierarchy.
      PDE_CursorFEside* create_CursorFEside( PEL_Object* a_owner ) const ;

      // Create and return a `PDE_LocalFEbound::' object that conforms with
      // the associated Module Hierarchy.
      PDE_LocalFEbound* create_LocalFEbound( PEL_Object* a_owner ) const ;

      // Create and return a `PDE_LocalFEinterface::' object that conforms
      // with the associated Module Hierarchy, for which all inward meshes
      // have the color `inward_domain' and all outward meshes have the color
      // `outward_domain'.
      PDE_LocalFEinterface* create_LocalFEinterface(
                               PEL_Object* a_owner,
                               GE_Color const* inward_domain,
                               GE_Color const* outward_domain ) const ;

   //-- Save for postprocessing(140)

      // the instance of  `PDE_ResultSaver::' that conforms with the
      // associated Module Hierarchy (any call on behalf of `self' will
      // provide the same instance, two calls on behalf of different
      // `PDE_DomainAndFields::' object will provides two different instances :
      // each `PDE_DomainAndFields::' object is in a one-to-one correspondance
      // with the `PDE_ResultSaver::' object returned by this method)
      PDE_ResultSaver* result_saver( void ) const ;

      // novel and unused instance of `PDE_ResultSaver::' built
      // according to the the Module Hierarchy attached to `exp'
      PDE_ResultSaver* novel_result_saver(
                                    PEL_ModuleExplorer const* exp ) const ;

   //-- Persistence

      virtual void save_state( PEL_ObjectWriter* writer ) const ;

      virtual void restore_state( PEL_ObjectReader* reader ) ;

   //-- Adaptation

      PDE_AdapterCHARMS* adapter_CHARMS( void ) const ;
      
      PDE_AdapterHN* adapter_HN( void ) const ;

   //-- Input - Output

      void print_grid( std::ostream& os, size_t indent_width ) const ;

   protected: //-------------------------------------------------------------

   private: //---------------------------------------------------------------

      PDE_DomainAndFields( void ) ;
     ~PDE_DomainAndFields( void ) ;
      PDE_DomainAndFields( PDE_DomainAndFields const& other ) ;
      PDE_DomainAndFields& operator=( PDE_DomainAndFields const& other ) ;

      PDE_DomainAndFields( PEL_Object* a_owner,
                           PEL_ModuleExplorer const* exp,
                           PEL_Communicator const* com ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Attributes

      PEL_ModuleExplorer const* EXP ;
      std::string DISC ;
      PDE_DomainBuilder* BUILDER_FE ;
      mutable PDE_ResultSaver* RSAVER ;
      mutable PEL_ModuleExplorer const* DECOHESION_EXP ;
      mutable PDE_AdapterCHARMS* DA ;
      mutable PDE_AdapterHN* TA ;
      PEL_Communicator const* COM ;
} ;

#endif
