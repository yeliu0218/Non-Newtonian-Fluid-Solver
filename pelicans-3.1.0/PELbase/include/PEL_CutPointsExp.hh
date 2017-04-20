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

#ifndef PEL_CUT_POINTS_EXP_HH
#define PEL_CUT_POINTS_EXP_HH

#include <PEL_Expression.hh>

/*
---
name     : middle_point
argument : Double, DoubleVector
type     : Double

middle_point(x,v) return the middle point of the interval defined by two
successive elements of v containing x.

Example:
   middle_point( 2.1, (1.,2.,3.,4.) ): value is 2.5

---
name     : middle_points
argument : DoubleVector[,DoubleVector]
type     : DoubleVector

middle_points(v) builds a vector containaing the middle points of
the intervals defined by two successive elements of v.

Example:
   middle_points( (1.,2.,3.,4.) ): value is (1.5,2.5,3.5)

middle_points(x_table,v) builds a vector containaing the middle points of
the intervals defined by two successive elements of v containing the elements
of x_table.

Example:
   middle_points( (1.1, 3.9), (1.,2.,3.,4.) ): value is (1.5,3.5)
   
---
name     : x_cut_points
argument : < x_table, y, [z] >
type     : doubleArray2D

x_cut_points(x_table,y) builds size(x_table)-1 points of 2 dimensions,
the i-th of which is (0.5*(x_table(i)+xtable(i+1)),y).
x_cut_points(x_table,y,z) builds size(x_table)-1 points of 3 dimensions,
the i-th of which is (0.5*(x_table(i)+xtable(i+1)),y,z).

example:
   x_cut_points((1.,2.,3.),3.9):
                 value is ( (1.5,3.9),(2.5,3.9) )
   x_cut_points((1.,2.,3.,5.),3.9,6.):
                 value is ( (1.5,3.9,6.),(2.5,3.9,6.),(4.,3.9,6.) )

---
name     : y_cut_points
argument : < x, y_table, [z] >
type     : doubleArray2D

y_cut_points(x,y_table) builds a set of size(y_table)-1 points of 2 dimensions,
the i-th of which is (x,0.5*(y_table(i)+ytable(i+1))).
y_cut_points(x,y_table,z) builds a set of size(y_table)-1 points of 3 dimensions,
the i-th of which is (x,0.5*(y_table(i)+ytable(i+1)),z).

example:
   y_cut_points(3.9,(1.,2.,3.),):
                 value is ( (3.9,1.5),(3.9,2.5) )
   y_cut_points(3.9,(1.,2.,3.,5.),6.):
                 value is ( (3.9,1.5,6.),(3.9,2.5,6.),(3.9,4.,6.) )

---
name     : z_cut_points
argument : < x, y, z_table >
type     : doubleArray2D

z_cut_points(x,y,z_table) builds a set of size(z_table)-1 points of 3 dimensions,
the i-th of which is (x,y,0.5*(z_table(i)+ztable(i+1))).

example:
   z_cut_points(3.9,6.,(1.,2.,3.,5.)):
                 value is ( (3.9,6.,1.5),(3.9,6.,2.5),(3.9,6.,4.) )
                 
PUBLISHED
*/

class PEL_EXPORT PEL_CutPointsExp : public PEL_Expression
{
   public: //----------------------------------------------------------

   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
      
   //-- Value

      virtual double to_double( PEL_Context const* ct ) const ;
      
      virtual doubleVector const& to_double_vector(
                                        PEL_Context const* ct ) const ;
      
      virtual doubleArray2D const& to_double_array2D(
                                        PEL_Context const* ct ) const ;
      
   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

      PEL_CutPointsExp( void ) ;
     ~PEL_CutPointsExp( void ) ;
      PEL_CutPointsExp( PEL_CutPointsExp const& other ) ;
      PEL_CutPointsExp& operator=( PEL_CutPointsExp const& other ) ;

      enum IS_CutOp { xcut, ycut, zcut, mid_point, mid_points } ;
      
      PEL_CutPointsExp( PEL_Object* a_owner,
                        std::string const& a_name,
                        PEL_Sequence const* argument_list,
                        IS_CutOp a_op ) ;

      void check_table( doubleVector const& verts_table ) const ;
      
      double m_pt( double x, doubleVector const& verts_table ) const ;
      
      void build_coords( doubleVector const& verts_table,
                         doubleVector& coords_table ) const ;
      
   //-- Plug in
      
      PEL_CutPointsExp( std::string const& a_name, IS_CutOp a_op ) ;

      virtual PEL_CutPointsExp* create_replica( 
                       PEL_Object* a_owner,
                       PEL_Sequence const* argument_list ) const ;
      
  //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments(
                      PEL_Sequence const* some_arguments ) const ;
      
   //-- Class attributes

      static PEL_CutPointsExp const* PROTOTYPE_MP ;
      static PEL_CutPointsExp const* PROTOTYPE_MPS ;
      static PEL_CutPointsExp const* PROTOTYPE_X ;
      static PEL_CutPointsExp const* PROTOTYPE_Y ;
      static PEL_CutPointsExp const* PROTOTYPE_Z ;
          
   //-- Attributes
      
      IS_CutOp const OP ;
} ;

#endif
