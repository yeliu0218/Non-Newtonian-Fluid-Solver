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

#ifndef PEL_INTERPOL_EXP_HH
#define PEL_INTERPOL_EXP_HH

#include <PEL_Expression.hh>
#include <doubleVector.hh>

/*
Expressions interpolating between given values.

---
name      : interpol
arguments : DoubleVector, DoubleVector, Double 
type      : double

interpol( x_values, fx_values, x ) computes at x the linear interpolation 
between the data points (x_values, fx_values).
More precisely:
   * x_values and fx_values are two vectors of the same size, say N
   * the elements of x_values are ordered by increasing values
        x_values(0) < x_values(1) < ... < x_values(N-1)
   * if x is smaller than x_values(0) then return fx_values(0)
     else if x is greater than x_values(N-1) then return fx_values(N-1)
     else find j such that x is between x_values(j) and x_values(j+1)
          return the linear interpolation between 
                 fx_values(j) and fx_values(j+1)

example :
   $DV_Xval = < 1. 2. 5. >
   $DV_Yval = < 1. 4. 2. >
   interpol( $DV_Xval, $DV_Yval, 3. ) : value is y given by 
                                        (y-4.0)/(3.0-2.0)=(2.0-4.0)/(5.0-2.0)

---
name      : interpol
arguments : String, Double 
type      : double

interpol( filename, x ) computes at x the linear interpolation between the
data points stored in file of name filename.
the value is identical as that of interpol( x_values, fx_values, x )
(see above) where filename stores in sequence x_values(0), fx_values(0), 
x_values(1), fx_values(1), ... , x_values(N-1), fx_values(N-1)

example :
   file "values.txt":
      1. 1.
      2. 4.
      5. 2.
   interpol( join( this_file_dir(), "values.txt" ), 3. ) : same value as above
      ie y given by (y-4.0)/(3.0-2.0)=(2.0-4.0)/(5.0-2.0)
   
PUBLISHED
*/

class PEL_EXPORT PEL_InterpolExp : public PEL_Expression
{
   public: //-------------------------------------------------------
      
   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual double to_double( PEL_Context const* ct = 0 ) const ;
      
   protected: //-------------------------------------------------------
      
   private: //-------------------------------------------------------

      PEL_InterpolExp( void ) ;
     ~PEL_InterpolExp( void ) ;
      PEL_InterpolExp( PEL_InterpolExp const& other ) ;
      PEL_InterpolExp& operator=( PEL_InterpolExp const& other ) ;

      enum PEL_InterpolOp { lin_inter_1D } ;

      PEL_InterpolExp( PEL_Object* a_owner,
                       std::string const& a_name,
                       PEL_Sequence const* argument_list,
                       PEL_InterpolOp a_op ) ;

      void read_tables_1D( std::string const& filename,
                           doubleVector& X_table,
                           doubleVector& FX_table ) const ;
      
      double linear_interpol_1D( doubleVector const& X_table,
                                 doubleVector const& FX_table,
                                 double x ) const ;

      void check_tables_1D( doubleVector const& X_table,
                            doubleVector const& FX_table ) const ;

   //-- Plug in

      PEL_InterpolExp( std::string const& a_name,
                       PEL_InterpolOp a_op ) ;

      virtual PEL_Expression* create_replica( 
                                   PEL_Object * a_owner,
                                   PEL_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( PEL_Sequence const* some_arguments ) const ;
      
   //-- Class attributes

      static PEL_InterpolExp const* PROTOTYPE_LIN_1D ;
      
   //-- Attributes

      PEL_InterpolOp const OP ;
      bool FROM_FILE ;
      mutable bool X1_IS_SET ;
      mutable doubleVector X1 ;
      mutable bool F_IS_SET ;
      mutable doubleVector FX1 ;
      mutable bool CHECK ;
} ;

#endif
