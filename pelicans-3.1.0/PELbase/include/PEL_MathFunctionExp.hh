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

#ifndef PEL_MATH_FUNCTION_HH
#define PEL_MATH_FUNCTION_HH

#include <PEL_Expression.hh>

/* 
Expressions defined by the usual mathematical functions :
      sqr, sqrt, pow, exp, log, log10,
      sin, cos, tan, asin, acos, atan, atan2,
      sinh, cosh, tanh, asinh, acosh, atanh
      j0, j1, jn, y0, y1, yn,
      gamma, lgamma, erf, erfc, incomplete_gamma, En, Ei,
      abs, floor, ceil
and by the comparison between double values :
      double_equality

PUBLISHED
*/

class PEL_EXPORT PEL_MathFunctionExp : public PEL_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
      
   //-- Value

      virtual bool to_bool( PEL_Context const* ct = 0 ) const ;
      
      virtual double to_double( PEL_Context const* ct = 0 ) const ;

      virtual int to_int( PEL_Context const* ct = 0 ) const ;
      
   //-- Formal calculus

      virtual PEL_Data* create_derivative( PEL_Object* a_owner,
                                           PEL_Variable const* var,
                                           PEL_Context const* ct ) const ;

   protected: //-------------------------------------------------------
      
   private: //-------------------------------------------------------

      PEL_MathFunctionExp( void ) ;
     ~PEL_MathFunctionExp( void ) ;
      PEL_MathFunctionExp( PEL_MathFunctionExp const& other ) ;
      PEL_MathFunctionExp& operator=( PEL_MathFunctionExp const& other ) ;

      enum MathOp { Sqr, Sqrt, Pow, Exp, Log, Log10,
                    Sin, Cos, Tan, ASin, ACos, ATan, ATan2,
                    Sinh, Cosh, Tanh, ASinh, ACosh, ATanh,
                    J0, J1, Jn, Y0, Y1, Yn,
                    Gamma, LGamma, Erf, Erfc,
                    IGamma, En, Ei, 
                    Abs, Floor, Ceil,
                    DblEq,
                    RandDbl, RandInt } ;
            
      PEL_MathFunctionExp( PEL_Object* a_owner,
                           std::string const& a_name,
                           PEL_Sequence const* argument_list,
                           MathOp a_op ) ;

      PEL_Data const* alternative_result( void ) const ;
      
   //-- Plug in

      PEL_MathFunctionExp( std::string const& a_name, MathOp a_op ) ;

      virtual PEL_MathFunctionExp* create_replica( 
                                   PEL_Object * a_owner,
                                   PEL_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( PEL_Sequence const* some_arguments ) const ;

      
   //-- Class attributes

      static PEL_MathFunctionExp const* PROTOTYPE_Sqr ;
      static PEL_MathFunctionExp const* PROTOTYPE_Sqrt ;
      static PEL_MathFunctionExp const* PROTOTYPE_Pow ;
      static PEL_MathFunctionExp const* PROTOTYPE_Exp ;
      static PEL_MathFunctionExp const* PROTOTYPE_Log ;
      static PEL_MathFunctionExp const* PROTOTYPE_Log10 ;
      static PEL_MathFunctionExp const* PROTOTYPE_Sin ;
      static PEL_MathFunctionExp const* PROTOTYPE_Cos ;
      static PEL_MathFunctionExp const* PROTOTYPE_Tan ;
      static PEL_MathFunctionExp const* PROTOTYPE_ASin ;
      static PEL_MathFunctionExp const* PROTOTYPE_ACos ;
      static PEL_MathFunctionExp const* PROTOTYPE_ATan ;
      static PEL_MathFunctionExp const* PROTOTYPE_ATan2 ;
      static PEL_MathFunctionExp const* PROTOTYPE_Sinh ;
      static PEL_MathFunctionExp const* PROTOTYPE_Cosh ;
      static PEL_MathFunctionExp const* PROTOTYPE_Tanh ;
      static PEL_MathFunctionExp const* PROTOTYPE_ASinh ;
      static PEL_MathFunctionExp const* PROTOTYPE_ACosh ;
      static PEL_MathFunctionExp const* PROTOTYPE_ATanh ;
      static PEL_MathFunctionExp const* PROTOTYPE_J0 ;
      static PEL_MathFunctionExp const* PROTOTYPE_J1 ;
      static PEL_MathFunctionExp const* PROTOTYPE_Jn ;
      static PEL_MathFunctionExp const* PROTOTYPE_Y0 ;
      static PEL_MathFunctionExp const* PROTOTYPE_Y1 ;
      static PEL_MathFunctionExp const* PROTOTYPE_Yn ;
      static PEL_MathFunctionExp const* PROTOTYPE_Gamma ;
      static PEL_MathFunctionExp const* PROTOTYPE_LGamma ;
      static PEL_MathFunctionExp const* PROTOTYPE_Erf ;
      static PEL_MathFunctionExp const* PROTOTYPE_Erfc ;
      static PEL_MathFunctionExp const* PROTOTYPE_IGamma ;
      static PEL_MathFunctionExp const* PROTOTYPE_En ;
      static PEL_MathFunctionExp const* PROTOTYPE_Ei ;
      static PEL_MathFunctionExp const* PROTOTYPE_ABS ;
      static PEL_MathFunctionExp const* PROTOTYPE_Floor ;
      static PEL_MathFunctionExp const* PROTOTYPE_Ceil ;
      static PEL_MathFunctionExp const* PROTOTYPE_DblEq ;
      static PEL_MathFunctionExp const* PROTOTYPE_RandDbl ;
      static PEL_MathFunctionExp const* PROTOTYPE_RandInt ;
      
   //-- Attributes
      
      MathOp const OP ;
      PEL_Data const* const ARG1 ;
      PEL_Data const* const ARG2 ;
      PEL_Data const* const ARG3 ;
      PEL_Data const* const ARG4 ;
      PEL_Data const* const ARG5 ;
} ;

#endif
