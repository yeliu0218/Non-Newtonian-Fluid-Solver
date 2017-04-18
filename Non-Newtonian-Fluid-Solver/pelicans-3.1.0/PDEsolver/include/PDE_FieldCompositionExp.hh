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

#ifndef PDE_FIELD_COMPOSITION_EXP_HH
#define PDE_FIELD_COMPOSITION_EXP_HH

#include <PDE_FieldComposition.hh>

#include <doubleVector.hh>
#include <stringVector.hh>

#include <string>

class PEL_ContextSimple ;
class PEL_Data ;
class PEL_Vector ;

class PDE_SetOfDiscreteFields ;
class PEL_Sequence ;

/*
PDE_FieldComposition objects built from expressions.

PUBLISHED
*/
   
class PEL_EXPORT PDE_FieldCompositionExp : public PDE_FieldComposition
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization
      
      // Create and return an instance defined by `expression'. The variables
      // of `expression' should be of type `PEL_Data::'DoubleVector. The
      // right part of their name should be the name of an item of `dfs'
      // or alternatively the name of another `PDE_FieldCompositionExp::' (for
      // composition).
      static PDE_FieldCompositionExp* create(
                                      PEL_Object* a_owner,
                                      std::string const& a_name,
                                      PEL_Data* expression,
                                      PDE_SetOfDiscreteFields const* dfs ) ;
 
      virtual void do_the_links( PDE_SetOfFieldCompositions const* fcs ) ;
      
   //-- Description
      
      virtual size_t nb_components( void ) const ;
      
   //-- Computation and retrieval of self value

      virtual double value( size_t ic ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------
      
      PDE_FieldCompositionExp( void ) ;
     ~PDE_FieldCompositionExp( void ) ;
      PDE_FieldCompositionExp( PDE_FieldCompositionExp const& other ) ;
      PDE_FieldCompositionExp& operator=( 
                               PDE_FieldCompositionExp const& other ) ;

   //-- Plug in

      PDE_FieldCompositionExp( PEL_Object* a_owner,
                               std::string const& a_name,
                               PEL_Data* expression,
                               PDE_SetOfDiscreteFields const* dfs ) ;
      
      virtual PDE_FieldComposition* create_replica(
                             PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp,
                             size_t nb_sp_dims,
                             PDE_SetOfDiscreteFields const* dfs ) const ;
      
   //-- Computation

      virtual void compute_self( void ) ;

   //-- Attributes

      bool IS_CSTE ;
      PEL_Data const* EXPR ;
      PEL_ContextSimple* CONTEXT ;
      PEL_Vector* FIELDS ;
      doubleVector VAL ;
      stringVector INNER_LAW_NAMES ;
      PEL_Vector* INNER_LAWS ;
} ;

#endif
