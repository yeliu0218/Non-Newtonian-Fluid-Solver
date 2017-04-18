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

#ifndef PDE_REFINEMENT_PATTERN_PROVIDER_HH
#define PDE_REFINEMENT_PATTERN_PROVIDER_HH

#include <PEL_Object.hh>

#include <map>

class GE_ReferencePolyhedron ;
class GE_ReferencePolyhedronRefiner ;

class PDE_CellFE ;
class PDE_ReferenceElement ;
class PDE_ReferenceElementRefiner ;

class PEL_EXPORT PDE_RefinementPatternProvider : public PEL_Object
{
   public: //------------------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_RefinementPatternProvider const* object( 
                                                  std::string const& a_name ) ;

   //-- Access to refiners

      GE_ReferencePolyhedronRefiner const* cell_refiner( 
                                      PDE_CellFE const* cell ) const ;

      PDE_ReferenceElementRefiner const* reference_element_refiner(
                                      PDE_ReferenceElement const* elm ) const ;

   protected: //---------------------------------------------------------------

   private: //-----------------------------------------------------------------

      PDE_RefinementPatternProvider( void ) ;
     ~PDE_RefinementPatternProvider( void ) ;
      PDE_RefinementPatternProvider( 
                          PDE_RefinementPatternProvider const& other ) ;
      PDE_RefinementPatternProvider& operator=( 
                          PDE_RefinementPatternProvider const& other ) ;

      PDE_RefinementPatternProvider( PEL_Object* a_owner ) ;

   //-- Internals

      void add_pattern( PDE_ReferenceElement const* elm,
                        PDE_ReferenceElementRefiner const* refiner ) ;

   //-- Attributes

      std::map< PDE_ReferenceElement const*,
                PDE_ReferenceElementRefiner const* > ELEM_RFS ;
      mutable std::map< PDE_ReferenceElement const*,
                PDE_ReferenceElementRefiner const* >::const_iterator IT_E ;

      std::map< GE_ReferencePolyhedron const*,
                GE_ReferencePolyhedronRefiner const* > CELL_RFS ;
      mutable std::map< GE_ReferencePolyhedron const*,
                GE_ReferencePolyhedronRefiner const* >::const_iterator IT_C ;
} ;

#endif
