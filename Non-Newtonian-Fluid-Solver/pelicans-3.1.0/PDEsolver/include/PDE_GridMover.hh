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

#ifndef PDE_GRID_MOVER_HH
#define PDE_GRID_MOVER_HH

#include <PEL_Object.hh>

class PEL_Vector ;

class GE_Point ;
class GE_SetOfPoints ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LocalFEcell ;

/*
Servers that strain the underlying grid of `PDE_DomainAndFields::' objects.

The vertices of the grid are moved, and then all geometrical status of the 
polyhedrons are updated.
  
PUBLISHED
*/

class PEL_EXPORT PDE_GridMover : public PEL_Object
{
   public: //----------------------------------------------------------  

   //-- Instance delivery and initialization

      static PDE_GridMover* create( PEL_Object* a_owner,
                                    PDE_DomainAndFields const* dom ) ;

   //-- Instance characteristics

      size_t nb_space_dimensions( void ) const ;
      
   //-- Move grid :

      // Move the computationnal grid with the `level'-th storage of
      // `strain' premultiplied by `coef'.
      void move_grid( PDE_DiscreteField const* strain, size_t level,
                      double coef = 1.0 ) ;
      
   protected: //-------------------------------------------------------

   private: //----------------------------------------------------------
      
      PDE_GridMover( void ) ;
     ~PDE_GridMover( void ) ;
      PDE_GridMover( PDE_GridMover const& other ) ;
      PDE_GridMover& operator=( PDE_GridMover const& other ) ;

      PDE_GridMover( PEL_Object* a_owner,
                     PDE_DomainAndFields const* dom ) ;

   //-- Attributes
   
      size_t NB_DIMS ;   
      GE_SetOfPoints const* const VERTICES ;
      PDE_LocalFEcell* cFE ;
      PEL_Vector* const DEFO ;
} ;


#endif
