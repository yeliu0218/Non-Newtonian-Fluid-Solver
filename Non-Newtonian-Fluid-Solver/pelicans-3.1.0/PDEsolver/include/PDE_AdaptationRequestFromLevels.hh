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

#ifndef PDE_ADAPTATION_REQUEST_FROM_LEVELS_HH
#define PDE_ADAPTATION_REQUEST_FROM_LEVELS_HH

#include <PDE_AdaptationRequest.hh>

class PEL_EXPORT PDE_AdaptationRequestFromLevels : public PDE_AdaptationRequest
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_AdaptationRequestFromLevels* create_for_refinement( 
                                                     PEL_Object* a_owner,
                                                     size_t target_level,
                                                     size_t verbose_level ) ;
 
      static PDE_AdaptationRequestFromLevels* create_for_unrefinement( 
                                                     PEL_Object* a_owner,
                                                     size_t target_level,
                                                     size_t verbose_level ) ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_AdaptationRequestFromLevels( void ) ;
     ~PDE_AdaptationRequestFromLevels( void ) ;
      PDE_AdaptationRequestFromLevels( 
                    PDE_AdaptationRequestFromLevels const& other ) ;
      PDE_AdaptationRequestFromLevels& operator=( 
                    PDE_AdaptationRequestFromLevels const& other ) ;

      PDE_AdaptationRequestFromLevels( PEL_Object* a_owner,
                                       size_t target_level,
                                       bool refine,
                                       size_t verbose_level ) ;
      
   //-- Internal configuration by applying adaptation rules
      
      virtual bool to_be_refined( PDE_CellFE const* cell,
                                  PDE_BasisFunctionCell const* bf,
                                  PDE_ReferenceElement const* elm,
                                  size_t local_node ) const ;

      virtual bool to_be_unrefined( PDE_CellFE const* cell,
                                    PDE_BasisFunctionCell const* bf,
                                    PDE_ReferenceElement const* elm,
                                    size_t local_node ) const ;
      
   //-- Attributes
      
      size_t TLEVEL ;
      bool REF ;
} ;

#endif
