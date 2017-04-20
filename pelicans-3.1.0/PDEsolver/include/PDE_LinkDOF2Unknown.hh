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

#ifndef PDE_LINK_DOF_2_UNKNOWN_HH
#define PDE_LINK_DOF_2_UNKNOWN_HH

#include <PEL_Object.hh>

#include <size_t_vector.hh>
#include <boolVector.hh>

#include <PDE_DiscreteField.hh>

class LA_Vector ;

class PDE_CrossProcessUnknownNumbering ;

/*
Arrangements of some of the DOFs of `PDE_DiscreteField::' objects in a linear
sequence (called unknown vector) where each DOF is identified by a 
sequential (ordinal) number (called unknown).

PUBLISHED
*/

class PEL_EXPORT PDE_LinkDOF2Unknown : public PEL_Object
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization
      
      // Create and return an instance.
      static PDE_LinkDOF2Unknown* create( PEL_Object* a_owner,
                                          PDE_DiscreteField const* ff,  
                                          std::string const& ordering,
                                          bool imposed_out ) ;

      // Create and return an instance.
      static PDE_LinkDOF2Unknown* create( PEL_Object* a_owner,
                                          PDE_DiscreteField const* ff,  
                                          size_t ic,
                                          bool imposed_out ) ;

      // Create and return an instance.
      static PDE_LinkDOF2Unknown* create( PEL_Object* a_owner,
                                          PDE_DiscreteField const* ff,  
                                          bool imposed_out ) ;

      virtual PDE_LinkDOF2Unknown* create_clone( PEL_Object* a_owner ) const ;

      // Reset the internal state, as if the create method used to instantiate
      // `self' was just completed (the field whose DOFs are handled might 
      // have changed since the creation of `self'). 
      void reset( void ) ;

      void reset( boolVector const& observed_nodes ) ;

   //-- Characteristics

      // field whose DOFs are handled
      PDE_DiscreteField const* field( void ) const ;

      // list of the DOF components that are handled
      size_t_vector const& components_table( void ) const ;

      // for the field whose DOFs are handled, number of nodes when `self'
      // was created, or just after calling `::reset'
      size_t nb_field_nodes( void ) const ;

      // Are the DOFs having an imposed value handled ?
      bool DOFs_with_imposed_value_are_dropped( void ) const ;

      // option determining the dependance of the linear sequence of indices
      // representing the handled DOFs with respect to the nodes and components
      // "sequence_of_the_nodes" means : the node is varied first, 
      //                                 then the component and so on
      // "sequence_of_the_components" means : the component is varied first, 
      //                                      then the node and so on
      std::string DOFs_ordering_in_unknown( void ) const ;

   //-- Distributed processing

      PDE_CrossProcessUnknownNumbering const* 
                                    cross_process_numbering( void ) const ;

   //-- Linear sequence of indices

      // number of indices in the linear sequence representing the handled DOFs
      size_t unknown_vector_size( void ) const ;

      // Does the DOF defined by node `n' and component `ic' is represented in 
      // the unknown vector ?
      bool DOF_is_unknown( size_t n, size_t ic=0 ) const ;

      // the index of the DOF defined by node `n' and component `ic' in the 
      // unknown vector
      size_t unknown_linked_to_DOF( size_t n, size_t ic=0 ) const ;
       
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_LinkDOF2Unknown( void ) ;
     ~PDE_LinkDOF2Unknown( void ) ;
      PDE_LinkDOF2Unknown( PDE_LinkDOF2Unknown const& other ) ;
      PDE_LinkDOF2Unknown& operator=( PDE_LinkDOF2Unknown const& other ) ;

      PDE_LinkDOF2Unknown( PEL_Object* a_owner,
                           PDE_DiscreteField const* ff,
                           std::string const& ordering,
                           bool imposed_out ) ;

      PDE_LinkDOF2Unknown( PEL_Object* a_owner,
                           PDE_DiscreteField const* ff,
                           size_t ic,
                           bool imposed_out ) ;

      PDE_LinkDOF2Unknown( PEL_Object* a_owner, 
                           PDE_LinkDOF2Unknown const* other ) ;

   //-- Internals

      enum OrderingType{ sequenceOfTheComponents, sequenceOfTheNodes } ;

      size_t local_index_of_DOF( size_t n, size_t idx_in_COMPS ) const ;

      void set_ordering_option( std::string const& option ) ;
      
   //-- Preconditions, Postconditions, Invariant      

      virtual bool invariant( void ) const ;      

   //-- Attributes

      PDE_DiscreteField const* FIELD ;
      size_t NB_NODES ;

      size_t NB_COMPS ;
      size_t_vector COMPS ;
      size_t_vector INDEX_in_COMPS ;

      bool DROP_IMPOSED ;
      OrderingType ORDERING ;

      size_t NB_DOFs ;
      size_t_vector DOF_2_UNKNOWN ;
      boolVector DOF_IN_UNKNOWNS ;

      PDE_CrossProcessUnknownNumbering* DIS_LINK ;
} ;

#ifndef OUTLINE
   #include <PDE_LinkDOF2Unknown.icc>
#endif

#endif
