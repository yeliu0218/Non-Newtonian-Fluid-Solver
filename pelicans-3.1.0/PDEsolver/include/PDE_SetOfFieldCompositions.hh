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

#ifndef PDE_SET_OF_FIELD_COMPOSITIONS_HH
#define PDE_SET_OF_FIELD_COMPOSITIONS_HH

#include <PEL_Object.hh>

#include <string>

class PDE_FieldComposition ;
class PDE_SetOfDiscreteFields ;

class PEL_ModuleExplorer ;
class PEL_Vector ;

/*
Sets of `PDE_CompositionOfFields::' instances, uniquely determined by their name. 
*/

class PEL_EXPORT PDE_SetOfFieldCompositions : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_SetOfFieldCompositions* create(  PEL_Object* a_owner ) ;

      // Build and append to `self' the `PDE_FieldComposition::' defined
      // in `exp'.
      void build_compositions( PEL_ModuleExplorer* exp,
                               size_t nb_sp_dimensions,
                               PDE_SetOfDiscreteFields const* fields ) ;

   //-- Access(50)

      // number of instances of PDE_FieldCompositions included
      size_t nb_compositions( void ) const ;

      // Does self include an instance of PDE_FieldComposition of 
      // name `name' ?
      bool has( std::string const& name ) const ;

      // the included instance of `PDE_FieldComposition::' of name `name'
      // ( raise an error if `name' is unknown )
      PDE_FieldComposition* item( std::string const& name ) const ;

   //-- Iteration over the included instances of `PDE_FieldComposition::'(70)

      // Move iterator to the first position.
      void start( void ) const ;

      // Is iterator position valid ? 
      bool is_valid( void ) const ;

      // Move iterator one position.
      void go_next( void ) const ;

      // instance of `PDE_FieldComposition::' at current iterator position
      PDE_FieldComposition* item( void ) const ;
      
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------
      
      PDE_SetOfFieldCompositions( void ) ;
     ~PDE_SetOfFieldCompositions( void ) ;
      PDE_SetOfFieldCompositions(
                   PDE_SetOfFieldCompositions const& other ) ;
      PDE_SetOfFieldCompositions& operator=(
                   PDE_SetOfFieldCompositions const& other ) ;

      PDE_SetOfFieldCompositions( PEL_Object* a_owner ) ;
      
      void append( PDE_FieldComposition* composition ) ;
      size_t find_index( std::string const& name ) const ;

   //-- Attributes

      PEL_Vector* const COMPO_TABLE ;
      mutable size_t INDEX ;
} ;

#endif // PDE_SET_OF_FIELD_COMPOSITIONS_HH
