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

#ifndef GE_PERTURBATED_MESHING_EXP_HH
#define GE_PERTURBATED_MESHING_EXP_HH

#include <PEL_Expression.hh>

#include <doubleVector.hh>

class PEL_Randomizer ;

/*
Expression for pseudo-randomly move points (eg vertices of a meshing)

name       : perturbated_coordinates
arguments  : Double, DoubleVector, Double, Bool
type       : DoubleVector

The return value `new_coords' of

   perturbated_coordinates( `epsilon', `coords', `h_min', `is_boundary' )

is such that:

   * `new_coords' has the same size as `coords'

   * if `is_boundary' is true: `new_coords' is equal to `coords' 

   * otherwise: `new_coords' is obtained by adding to `coords' a displacement
         of length `epsilon'*`h_min' in a pseudo-random direction (which is
         computed in a sequence of pseudo-random numbers)

Possible Usage
--------------

   `perturbated_coordinates' may be used as the transformation expression in 
   `GE_TransformedMeshing::' to generated a pseudo-randomly perturbated 
   meshing.

   Example:
    
      MODULE GE_Meshing
         concrete_name = "GE_TransformedMeshing"
         transformation = perturbated_coordinates(
                          0.2, // Perturbation coefficient
                          $DV_X,
                          $DS_VERT_HMIN, $BS_VERT_ON_BOUND )
         ...
      END MODULE GE_Meshing

   in this case 
      `$DV_X' is the coordinates of the current vertex
      `$DS_VERT_HMIN' is the minimal distance between the current vertex and
         its neighbor vertices
      `$BS_VERT_ON_BOUND' is `true' iff the current vertex is located on the
         domain boundary.
      `epsilon' should be such that (0.<`epsilon'<0.5)

PUBLISHED
*/

class PEL_EXPORT GE_PerturbatedMeshingExp : public PEL_Expression
{
   public: //---------------------------------------------------------

   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual doubleVector const& to_double_vector( PEL_Context const* ct ) const ;

   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------
      
      GE_PerturbatedMeshingExp( void ) ;
     ~GE_PerturbatedMeshingExp( void ) ;
      GE_PerturbatedMeshingExp( GE_PerturbatedMeshingExp const& other ) ;
      GE_PerturbatedMeshingExp& operator=(
                                GE_PerturbatedMeshingExp const& other ) ;

      enum PerturbatedMeshingExp{ pert_coords } ;
         
      GE_PerturbatedMeshingExp( PEL_Object* a_owner,
                                PerturbatedMeshingExp id,
                                std::string const& a_name,
                                PEL_Sequence const* argument_list ) ;

      double random_value( void ) const ;
      
   //-- Plug in

      GE_PerturbatedMeshingExp( PerturbatedMeshingExp id,
                                std::string const& a_name ) ;
      
      virtual GE_PerturbatedMeshingExp* create_replica( 
                                PEL_Object * a_owner,
                                PEL_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( PEL_Sequence const* some_arguments ) const ;
      
   //-- Class attributes

      static GE_PerturbatedMeshingExp const* PROTOTYPE_PERT_COORDS ;

   //-- Attributes

      PerturbatedMeshingExp ID ;
      PEL_Randomizer* const RAND ;
} ;

#endif
