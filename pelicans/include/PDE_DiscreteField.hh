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

#ifndef PDE_DISCRETE_FIELD_HH
#define PDE_DISCRETE_FIELD_HH

#include <PEL_Object.hh>

#include <string>

#include <boolVector.hh>
#include <size_t_array2D.hh>
#include <doubleArray3D.hh>
#include <doubleVector.hh>

class PEL_ModuleExplorer ;
class PEL_Vector ;
class stringVector ;

class LA_SeqVector ;

class PDE_CrossProcessNodeNumbering ;
class PDE_DOFconstraintsIterator ;
class PDE_DOFconstraints ;
class PDE_LinkDOF2Unknown ;

/*
Discrete fields.

A discrete field is a set of DOFs.

A DOF (standing for Degree Of Freedom) is defined by two notions :
     1. a Node ;
     2. a component.
  A Node has nothing to do with a geometric location. It is merely a shortcut
     for an index.
  Each Node might be tagged active or not.

For each DOF :
   - an amount of storage_depth() values of type double are stored
   - a supplementary value, called imposed value, might possibly be stored.
  Full access to these values is provided.

*/

class PEL_EXPORT PDE_DiscreteField : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_DiscreteField* create( PEL_Object* a_owner,
                                        std::string const& a_name,
                                        size_t a_nb_components,
                                        size_t a_depth ) ;

      virtual PDE_DiscreteField* create_duplication( 
                                        PEL_Object* a_owner,
                                        std::string const& a_name ) const ;

      virtual PDE_DiscreteField* create_duplication( 
                                        PEL_Object* a_owner,
                                        std::string const& a_name,
                                        size_t a_nb_components,
                                        size_t a_depth ) const ;

      static size_t nb_objects( void ) ;

      void set_nb_nodes( size_t a_nb_nodes ) ;
      
      void add_nodes( size_t nb_supp_nodes ) ;

   //-- Identification

      // number, uniquely determining `self'
      size_t id_number( void ) const ;

      // name (two instances might have the same name)
      std::string const& name( void ) const ;
      
   //-- Distributed processing

      // Does `self' represent a part of a logical global discrete field
      // involved in a processing distributed over several processes ?
      bool is_distributed( void ) const ;

      // cross-process numbering
      PDE_CrossProcessNodeNumbering* cross_process_numbering( void ) const ;
      
      // Perform all necessary tasks to reinterpret `self' as a part of 
      // a logical global discrete field in a processing distributed over
      // `numbering'->communicator().
      void build_cross_process_globalization( 
                               PDE_CrossProcessNodeNumbering* numbering ) ;

   //-- DOFs status(130)

      // number of nodes
      size_t nb_nodes( void ) const ;

      // number of DOFs at each node 
      size_t nb_components( void ) const ;
      
      // Is node `n' tagged active ?
      bool node_is_active( size_t n ) const ;

      // Tag node `n' active.
      void set_node_active( size_t n ) ;

      // Tag node `n' not active.
      void set_node_inactive( size_t n ) ;

      // Is the DOF defined by node `n' and component `ic' fixed ?
      bool DOF_is_fixed( size_t n, size_t ic ) const ;

      // Is the DOF defined by node `n' and component `ic' free ?
      bool DOF_is_free( size_t n, size_t ic ) const ;

      // number of values of type double stored for each DOF 
      // (other than imposed values)
      size_t storage_depth( void ) const ;

      // Does the DOF defined by node `n' and component `ic' 
      // has an imposed value ?
      bool DOF_has_imposed_value( size_t n, size_t ic=0 ) const ;

   //-- DOF modification(140)

      // Modify `self' by copying the DOF characteristics of `other'.
      void set( PDE_DiscreteField const* other ) ;

      // Modify the stored values of the DOF defined by node `n' 
      // and component `ic'.
      void modify_DOF( bool imposed, double value, size_t n, size_t ic=0 ) ;

      void equip_DOF_with_imposed_value( size_t n, double x, size_t ic=0 ) ;
      
   //-- Access to the DOF values(150)

      // imposed value of the DOF defined by node `n' and component `ic'
      double DOF_imposed_value( size_t n, size_t ic=0 ) const ;
      
      // `level'th storage value of the DOF defined by 
      // node `n' and component `ic'
      double DOF_value( size_t level, size_t n, size_t ic=0 ) const ;

      // For the `level'-th storage and the `ic'-th component, extract the 
      // value of all DOFs and copy it in `vec' at the index of the 
      // associated node.
      void extract_DOFs_value( size_t level, 
                               LA_SeqVector* vec, 
                               size_t ic = 0 ) const ;

      // For the `level'-th storage, extract the value of all DOFs that are 
      // linked to an unknown and copy it in `vec' at the index of the 
      // associated unknown as specified by `link'.
      void extract_unknown_DOFs_value( size_t level,
                                       LA_SeqVector* vec,
                                       PDE_LinkDOF2Unknown const* link ) const ;

   //-- Setting DOF values(150)

      // Set the `level'-th storage value of the DOF defined by node `n' and 
      // component `ic'.
      void set_DOF_value( size_t level, size_t n, double x, size_t ic=0 ) ;

      // Set the imposed value of the DOF defined by node `n' and 
      // component `ic'.
      void set_DOF_imposed_value( size_t n, double x, size_t ic=0 ) ;

      // For the `level'-th storage, set the value of all DOFs that are 
      // linked to an unknown equal to the value of the `vec' element at the  
      // index of the associated unknown as specified by link.
      void update_DOFs_value( size_t level,
                              LA_SeqVector const* vec,
                              PDE_LinkDOF2Unknown const* link ) ;

      // For the `level'-th storage, set the value of all free DOFs that are 
      // linked to an unknown equal to the value of the `vec' element at the 
      // index of the associated unknown as specified by link.
      void update_free_DOFs_value( size_t level,
                                   LA_SeqVector const* vec,
                                   PDE_LinkDOF2Unknown const* link ) ;

      // Do the same as `::update_free_DOFs_value' but instead of using `vec'
      // to set the value of all free DOFs, add `alpha'*`vec' to
      // the already existing values of all free DOFs (`level'-th storage).
      void add_to_free_DOFs_value( size_t level,
                                   LA_SeqVector const* vec,
                                   PDE_LinkDOF2Unknown const* link,
                                   double alpha ) ;

   //-- Copying DOF values between different storage levels(160)

      // Set the value of `level'-th storage of all DOFs equal to the
      // imposed value if any.
      void enforce_imposed_values_to_DOFs( size_t level ) ;

      // Set the value of the `target_level'-th storage of all DOFs equal to
      // that of the `source_level'-th storage.
      void copy_DOFs_value( size_t source_level, size_t target_level ) ;

      // Set the value of the `target_level'-th storage of all DOFs that
      // are linked to an unknown (as specified by `link') equal to
      // that of the `source_level'-th storage.
      void copy_unknown_DOFs_value( size_t source_level,
                                    size_t target_level,
                                    PDE_LinkDOF2Unknown const* link ) ;

      // Set the value of the `target_level'-th storage of all free DOFs that
      // are linked to an unknown (as specified by `link') equal to
      // that of the `source_level'-th storage.
      void copy_unknown_free_DOFs_value( size_t source_level,
                                         size_t target_level,
                                         PDE_LinkDOF2Unknown const* link ) ;
      
   //-- Constraining
      
      void remove_constraint_for_DOF( size_t slave_n, size_t slave_ic ) ;
      
      void add_constraint_for_DOF( size_t slave_n, size_t slave_ic,
                                   size_t master_n, size_t master_ic, 
                                   double coef ) ;
      
      bool DOF_is_constrained( size_t n, size_t ic ) const ;
      
      PDE_DOFconstraintsIterator* create_constraints_iterator( 
                                         PEL_Object* a_owner ) const ;
            
      void enforce_constraints_for_DOFs( size_t level ) ;
            
   //-- Persistence   

      virtual void save_state( PEL_ObjectWriter* writer ) const ;

      virtual void restore_state( PEL_ObjectReader* reader ) ;

   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_DiscreteField( void ) ;
     ~PDE_DiscreteField( void ) ;
      PDE_DiscreteField( PDE_DiscreteField const& other ) ;
      PDE_DiscreteField& operator=( PDE_DiscreteField const& other ) ;

      PDE_DiscreteField( PEL_Object* a_owner,
                         std::string const& a_name,
                         size_t a_nb_components,
                         size_t a_depth ) ;

      PDE_DiscreteField( PEL_Object* a_owner,
                         PDE_DiscreteField const* other,
                         std::string const& a_name ) ;
      
      PDE_DiscreteField( PEL_Object* a_owner,
                         PDE_DiscreteField const* other,
                         std::string const& a_name,
                         size_t a_nb_components,
                         size_t a_depth ) ;

   //-- Class attributes

      static size_t NB_INSTANCES ;
      
   //-- Attributes

      std::string const FNAME ;
      size_t const ID ;

      size_t NB_COMPS ;
      size_t NB_NODES ;
      size_t STO_DEPTH ;
    
      doubleArray3D VALUES ;
 
      size_t_array2D DBC_IDX ;
      boolVector     DBC_FLAG ;
      doubleVector   DBC_VALUES ;

      boolVector ACTIVE_NODES ;

      PDE_CrossProcessNodeNumbering* DF ;
      
      PDE_DOFconstraints* CSTR ;
      PDE_DOFconstraintsIterator* CSTR_IT ;
} ;


#endif

