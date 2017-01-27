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

#ifndef GE_REFERENCE_SQUARE_WITH_TRIANGLES_HH
#define GE_REFERENCE_SQUARE_WITH_TRIANGLES_HH

#include <GE_ReferencePolyhedronRefiner.hh>

/*
PUBLISHED
*/ 

class PEL_EXPORT GE_ReferenceSquareWithTriangles : public GE_ReferencePolyhedronRefiner
{
   public: //-----------------------------------------------------------

   //-- Fine polyhedra

      virtual void compute_location_in_subcell( 
                                       size_t ic, 
                                       GE_Point const* pt_cell,
                                       GE_Point* pt_subcell ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~GE_ReferenceSquareWithTriangles( void ) ;
      GE_ReferenceSquareWithTriangles( 
                            GE_ReferenceSquareWithTriangles const& other ) ;
      GE_ReferenceSquareWithTriangles& operator=( 
                            GE_ReferenceSquareWithTriangles const& other ) ;

      GE_ReferenceSquareWithTriangles( PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp ) ;
      
   //-- Plug in
      
      GE_ReferenceSquareWithTriangles( void ) ;
      
      virtual GE_ReferenceSquareWithTriangles const* create_replica(
                                    PEL_Object* a_owner,
                                    PEL_ModuleExplorer const* exp ) const ;      

   //-- Class attributes

      static GE_ReferenceSquareWithTriangles const* PROTOTYPE ;

   //-- Internals

      void initialize_right( void ) ;

      void initialize_left( void ) ;

      void initialize_cross( void ) ;

      void initialize_26_acute_triangles( void ) ;
} ;

#endif
