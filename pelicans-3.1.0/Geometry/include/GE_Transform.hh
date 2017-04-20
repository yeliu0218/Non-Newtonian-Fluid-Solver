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

#ifndef GE_TRANSFORM_HH
#define GE_TRANSFORM_HH

#include <PEL_Object.hh>

class PEL_ModuleExplorer ;

class GE_Color ;
class GE_Point ;
class GE_Vector ;

/*
PUBLISHED
*/

class PEL_EXPORT GE_Transform : public PEL_Object
{
   public: //----------------------------------------------------------

   //-- Instance delivery and initialization

      static GE_Transform* create( PEL_Object* a_owner,
                                   size_t nb_sp_dims,
                                   PEL_ModuleExplorer const* exp ) ;
      
      virtual GE_Transform const* inverse( void ) const = 0 ;
      
   //-- Status

      size_t nb_space_dimensions( void ) const ;
      
      GE_Color const* source_color( void ) const ;
      
      GE_Color const* target_color( void ) const ;
      
   //-- Transformation

      virtual void apply( GE_Point* pt ) const = 0 ;

   protected: //--------------------------------------------------------   

      virtual ~GE_Transform( void ) ;
      
      GE_Transform( PEL_Object* a_owner,
                    size_t nb_sp_dims,
                    PEL_ModuleExplorer const* exp ) ;
      
      GE_Transform( PEL_Object* a_owner ) ;
      
      void initialize_as_inverse( GE_Transform const* other ) ;
            
   //-- Preconditions, Postconditions, Invariant
      
      virtual bool inverse_POST( GE_Transform const* result ) const ;
      
      virtual bool apply_PRE( GE_Point* pt ) const ;

   private: //----------------------------------------------------------

      GE_Transform( void ) ;
      GE_Transform( GE_Transform const& other ) ;
      GE_Transform& operator=( GE_Transform const& other ) ;
      
   //-- Attributes
      
      size_t NB_DIMS ;
      GE_Color const* S_COLOR ;
      GE_Color const* T_COLOR ;
} ;

#endif
