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

#ifndef PDE_DOF_CONTRAINTS_HH
#define PDE_DOF_CONTRAINTS_HH

#include <PEL_Object.hh>

#include <vector>

#include <PDE_DiscreteField.hh> 

#include <size_t_array2D.hh>

class PDE_DOFconstraintsIterator ;

class PEL_EXPORT PDE_DOFconstraints : public PEL_Object
{
   public: //-----------------------------------------------------------
      
   //-- Indices
      
      void raise_nb_nodes( size_t a_nb_nodes ) ;
      
      size_t nb_nodes( void ) const ;
      
      size_t nb_components( void ) const ;
      
   //-- Constraints
      
      void remove( size_t slave_n, size_t slave_ic ) ;
      
      void add( size_t slave_n, size_t slave_ic,
                size_t master_n, size_t master_ic, 
                double coef ) ;
      
      bool has( size_t n, size_t ic ) const ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_DOFconstraints( void ) ;
     ~PDE_DOFconstraints( void ) ;
      PDE_DOFconstraints( PDE_DOFconstraints const& other ) ;
      PDE_DOFconstraints& operator=( PDE_DOFconstraints const& other ) ;
      
      PDE_DOFconstraints( PEL_Object* a_owner,
                          size_t a_nb_nodes,
                          size_t a_nb_comps ) ;
   //-- Friends
      
      friend class PDE_DOFconstraintsIterator ;
      
      friend void PDE_DiscreteField::add_constraint_for_DOF(
                                         size_t slave_n, size_t slave_ic,
                                         size_t master_n, size_t master_ic, 
                                         double coef ) ;
                  
   //-- Internals
      
      struct ConstraintElement
      {
          ConstraintElement( size_t master_n, size_t master_ic, double xx )
             : n( master_n ), ic( master_ic ), coef( xx ) {}
          size_t n ;
          size_t ic ;
          double coef ;
      } ;
      
   //-- Attributes

      size_t NB_NODES ;
      size_t NB_COMPS ;
      std::vector< std::vector< ConstraintElement > > CONSTRAINTS ;
      size_t_array2D CONSTRAINT_IDX ;
} ;


#endif

