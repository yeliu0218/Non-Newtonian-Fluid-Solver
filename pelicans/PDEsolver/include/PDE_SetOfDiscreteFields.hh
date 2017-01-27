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

#ifndef PDE_SET_OF_DISCRETE_FIELDS_HH
#define PDE_SET_OF_DISCRETE_FIELDS_HH

#include <PEL_Object.hh>
#include <boolVector.hh>

#include <string>

class PDE_DiscreteField ;
class PEL_Vector ;

/*
Sets of `PDE_DiscreteField::' instances, uniquely determined by their name. 
*/

class PEL_EXPORT PDE_SetOfDiscreteFields : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance.
      static PDE_SetOfDiscreteFields* create( PEL_Object* a_owner ) ;

   //-- Access(50)

      // number of instances of PDE_DiscreteField included
      size_t nb_fields( void ) const ;
      
      // the included instance of PDE_DiscreteField of index `i'
      PDE_DiscreteField* item( size_t i ) const ;
      
      bool boundary_located( size_t i ) const ;
      
      PEL_Vector const* reference_elements( size_t i ) const ;

      // Does self include an instance of PDE_DiscreteField of 
      // name `field_name' ?
      bool has( std::string const& field_name ) const ;

      // the included instance of PDE_DiscreteField of name `field_name'
      // ( raise an error if `field_name' is unknown )
      PDE_DiscreteField* item( std::string const& field_name ) const ;
      
      size_t index_of( std::string const& field_name ) const ;

   //-- Element change(60)

      // Add `field'.
      void append( PDE_DiscreteField* field,
                   bool located_on_bounds,
                   PEL_Vector* elements ) ;

   //-- Iteration over the included instances of PDE_DiscreteField(70)

      // Move iterator to the first position.
      void start( void ) const ;

      // Is iterator position valid ? 
      bool is_valid( void ) const ;

      // Move iterator one position.
      void go_next( void ) const ;

      // instance of PDE_DiscreteField at current iterator position
      PDE_DiscreteField* item( void ) const ;

   //-- Persistence   

      virtual void save_state( PEL_ObjectWriter* writer ) const ;

      virtual void restore_state( PEL_ObjectReader* reader ) ;
      
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_SetOfDiscreteFields( void ) ;
     ~PDE_SetOfDiscreteFields( void ) ;
      PDE_SetOfDiscreteFields( PDE_SetOfDiscreteFields const& other ) ;   
      PDE_SetOfDiscreteFields& operator=(
                               PDE_SetOfDiscreteFields const& other ) ;

      PDE_SetOfDiscreteFields( PEL_Object* a_owner ) ;

      size_t find_index( std::string const& field_name ) const ;

   //-- Attributes

      PEL_Vector* FIELDS ;
      PEL_Vector* ELMS ;
      boolVector ON_BOUNDS ;
      mutable size_t INDEX ;
} ;

#endif
