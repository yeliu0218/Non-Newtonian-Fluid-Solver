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

#ifndef PDE_FACES_OF_CELL_FE_HH
#define PDE_FACES_OF_CELL_FE_HH

#include <PDE_CellFE.hh>
#include <vector>

class PEL_VectorIterator ;

class PEL_EXPORT PDE_FacesOfCellFE : public PEL_Object
{
   public: //---------------------------------------------------------------

      void start( void ) ;
      
      bool is_valid( void ) const ;

      void go_next( void ) ;

      PDE_FaceFE* item( void ) const ;
             
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      PDE_FacesOfCellFE( void ) ;
     ~PDE_FacesOfCellFE( void ) ;
      PDE_FacesOfCellFE( PDE_FacesOfCellFE const& other ) ;
      PDE_FacesOfCellFE& operator=( PDE_FacesOfCellFE const& other ) ;
      
      PDE_FacesOfCellFE( PEL_Object* a_owner, PEL_Vector const* faces ) ;
      
      friend PDE_FacesOfCellFE* PDE_CellFE::create_faces_iterator( 
                                    PEL_Object* a_owner ) const ;
      
      static void extend_faces( std::vector< PDE_FaceFE* >& vec_of_faces,
                                PDE_FaceFE* a_face ) ;

   //-- Attributes

      std::vector< PDE_FaceFE* > ALL_FACES ;
      std::vector< PDE_FaceFE* >::const_iterator IT ;
      bool OK ;
} ;

#endif
