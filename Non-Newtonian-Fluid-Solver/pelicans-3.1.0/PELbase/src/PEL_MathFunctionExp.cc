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

#include <PEL_MathFunctionExp.hh>

#include <PEL_assertions.hh>
#include <PEL_Double.hh>
#include <PEL_Error.hh>
#include <PEL_Int.hh>
#include <PEL_List.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <iostream>

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Sqr = new PEL_MathFunctionExp( "sqr", Sqr ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Sqrt = new PEL_MathFunctionExp( "sqrt", Sqrt ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Pow = new PEL_MathFunctionExp( "pow", Pow ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Exp = new PEL_MathFunctionExp( "exp", Exp ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Log = new PEL_MathFunctionExp( "log", Log ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Log10 = new PEL_MathFunctionExp( "log10", Log10 ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Sin = new PEL_MathFunctionExp( "sin", Sin ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Cos = new PEL_MathFunctionExp( "cos", Cos ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Tan = new PEL_MathFunctionExp( "tan", Tan ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_ASin = new PEL_MathFunctionExp( "asin", ASin ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_ACos = new PEL_MathFunctionExp( "acos", ACos ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_ATan = new PEL_MathFunctionExp( "atan", ATan ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_ATan2 = new PEL_MathFunctionExp( "atan2", ATan2 ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Sinh = new PEL_MathFunctionExp( "sinh", Sinh ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Cosh = new PEL_MathFunctionExp( "cosh", Cosh ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Tanh = new PEL_MathFunctionExp( "tanh", Tanh ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_ASinh = new PEL_MathFunctionExp( "asinh", ASinh ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_ACosh = new PEL_MathFunctionExp( "acosh", ACosh ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_ATanh = new PEL_MathFunctionExp( "atanh", ATanh ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_J0 = new PEL_MathFunctionExp( "j0", J0 ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_J1 = new PEL_MathFunctionExp( "j1", J1 ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Jn = new PEL_MathFunctionExp( "jn", Jn ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Y0 = new PEL_MathFunctionExp( "y0", Y0 ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Y1 = new PEL_MathFunctionExp( "y1", Y1 ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Yn = new PEL_MathFunctionExp( "yn", Yn ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Gamma = new PEL_MathFunctionExp( "gamma", Gamma ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_LGamma = new PEL_MathFunctionExp( "lgamma", LGamma ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Erf = new PEL_MathFunctionExp( "erf", Erf ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Erfc = new PEL_MathFunctionExp( "erfc", Erfc ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_IGamma =
                        new PEL_MathFunctionExp( "incomplete_gamma", IGamma ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_En = new PEL_MathFunctionExp( "En", En ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Ei = new PEL_MathFunctionExp( "Ei", Ei ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_ABS = new PEL_MathFunctionExp( "abs", Abs ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Floor = new PEL_MathFunctionExp( "floor", Floor ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_Ceil = new PEL_MathFunctionExp( "ceil", Ceil ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_DblEq = new PEL_MathFunctionExp( "double_equality", DblEq ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_RandDbl = new PEL_MathFunctionExp( "random_double", RandDbl ) ;

PEL_MathFunctionExp const* 
PEL_MathFunctionExp::PROTOTYPE_RandInt = new PEL_MathFunctionExp( "rand", RandInt ) ;

struct PEL_MathFunctionExp_ERROR
{
   static void n0( std::string const& f_name, std::string const& op_name ) ;
} ;

//----------------------------------------------------------------------
PEL_MathFunctionExp:: PEL_MathFunctionExp( std::string const& a_name,
                                           MathOp a_op ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , OP( a_op )
   , ARG1( 0 )
   , ARG2( 0 )
   , ARG3( 0 )
   , ARG4( 0 )
   , ARG5( 0 )
{
   PEL_LABEL( "PEL_MathFunctionExp:: PEL_MathFunctionExp" ) ;
}

//----------------------------------------------------------------------
PEL_MathFunctionExp*
PEL_MathFunctionExp:: create_replica(
                             PEL_Object* a_owner,
                             PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MathFunctionExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_MathFunctionExp* result = new PEL_MathFunctionExp( a_owner, 
                                                          name(), 
                                                          argument_list, 
                                                          OP ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_MathFunctionExp:: PEL_MathFunctionExp( PEL_Object* a_owner,
                                           std::string const& a_name,
                                           PEL_Sequence const* argument_list,
                                           MathOp a_op ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , OP( a_op )
   , ARG1( nb_arguments()>0 ? arg(0) : 0 )
   , ARG2( nb_arguments()>1 ? arg(1) : 0 )
   , ARG3( nb_arguments()>2 ? arg(2) : 0 )
   , ARG4( nb_arguments()>3 ? arg(3) : 0 )
   , ARG5( nb_arguments()>4 ? arg(4) : 0 )
{
   PEL_LABEL( "PEL_MathFunctionExp:: PEL_MathFunctionExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_MathFunctionExp:: ~PEL_MathFunctionExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MathFunctionExp:: ~PEL_MathFunctionExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_MathFunctionExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MathFunctionExp:: data_type" ) ;

   PEL_Data::Type result = PEL_Data::Undefined ;

   switch( OP )
   {
      case Sqr :
      case Sqrt :
      case Pow :
      case Exp :
      case Log :
      case Log10 :
      case Sin :
      case Cos :
      case Tan :
      case ASin :
      case ACos :
      case ATan :
      case ATan2 :
      case Sinh :
      case Cosh :
      case Tanh :
      case ASinh :
      case ACosh :
      case ATanh :
      case J0 :
      case J1 :
      case Jn :
      case Y0 :
      case Y1 :
      case Yn :
      case Gamma :
      case LGamma :
      case IGamma :
      case Erf :
      case Erfc :
      case En :
      case Ei :
      case Floor :
      case Ceil :
      case RandDbl :
         result = PEL_Data::Double ;
         break ;
      case Abs :
         result = ARG1->data_type() ;
         break ;
      case DblEq :
         result = PEL_Data::Bool ;
         break ;
      case RandInt :
         result = PEL_Data::Int ;
         break ;
      default :
         PEL_MathFunctionExp_ERROR::n0( "data_type", name() ) ;
         break ;
   } 
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
PEL_MathFunctionExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result = "undefined" ;

   switch( OP )
   {
      case Sqr :
      case Sqrt :
      case Exp :
      case Log :
      case Log10 :
      case Sin :
      case Cos :
      case Tan :
      case ASin :
      case ACos :
      case ATan :
      case Sinh :
      case Cosh :
      case Tanh :
      case ASinh :
      case ACosh :
      case ATanh :
      case J0 :
      case J1 :
      case Y0 :
      case Y1 :
      case Gamma :
      case LGamma :
      case Erf :
      case Erfc :
      case Floor :
      case Ceil :
         result = name() + "(DS)" ;
         break ;
      case Jn :
      case Yn :
         result = name() + "(IS,DS)" ;
         break ;
      case Pow :
      case ATan2 :
         result =  name() + "(DS,DS)" ;
         break ;
      case IGamma :
         result =  name() + "(DS,DS[,DS,DS,IS])" ;
         break ;
      case En :
         result =  name() + "(IS,DS[,DS,DS,IS])" ;
         break ;
      case Ei :
         result =  name() + "(DS[,DS,DS,IS])" ;
         break ;
      case Abs :
         result =  name() + "(DS|IS)" ;
         break ;
      case DblEq :
         result = name() + "(DS,DS,DS,DS)" ;
         break ;
      case RandDbl :
      case RandInt :
         result = name() + "()" ;
         break ;
      default :
         PEL_MathFunctionExp_ERROR::n0( "usage", name() ) ;
         break ;
   }    
   
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_MathFunctionExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MathFunctionExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = false ;
   
   switch( OP )
   {
      case Sqr :
      case Sqrt :
      case Exp :
      case Log :
      case Log10 :
      case Sin :
      case Cos :
      case Tan :
      case ASin :
      case ACos :
      case ATan :
      case Sinh :
      case Cosh :
      case Tanh :
      case ASinh :
      case ACosh :
      case ATanh :
      case J0 :
      case J1 :
      case Y0 :
      case Y1 :
      case Gamma :
      case LGamma :
      case Erf :
      case Erfc :
      case Floor :
      case Ceil :
         result = ( some_arguments->count()==1 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            result = t0==Double ;
         }
         break ;
      case Abs :
         result = ( some_arguments->count()==1 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            result = t0==Double || t0==Int ;
         }
         break ;
      case Jn :
      case Yn :
         result = ( some_arguments->count()==2 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
            result = ( t0==Int && t1==Double ) ;
         }
         break ;
      case Pow :
      case ATan2 :  
         result = ( some_arguments->count()==2 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
            result = ( t0==Double && t1==Double ) ;
         }
         break ;
      case IGamma :
         result = ( some_arguments->count()==2 ||
                    some_arguments->count()==5 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
            result = ( t0==Double && t1==Double ) ;
            if( result && some_arguments->count()==5 )
            {
               Type t2 = extract_arg( some_arguments, 2 )->data_type() ;
               Type t3 = extract_arg( some_arguments, 3 )->data_type() ;
               Type k4 = extract_arg( some_arguments, 4 )->data_type() ;
               result = ( t2==Double && t3==Double && k4==Int ) ;
            }
         }
         break ;
      case En :
         result = ( some_arguments->count()==2 ||
                    some_arguments->count()==5 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
            result = ( t0==Int && t1==Double ) ;
            if( result && some_arguments->count()==5 )
            {
               Type t2 = extract_arg( some_arguments, 2 )->data_type() ;
               Type t3 = extract_arg( some_arguments, 3 )->data_type() ;
               Type k4 = extract_arg( some_arguments, 4 )->data_type() ;
               result = ( t2==Double && t3==Double && k4==Int ) ;
            }
         }
         break ;
      case Ei :
         result = ( some_arguments->count()==1 ||
                    some_arguments->count()==4 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            result = ( t0==Double ) ;
            if( result && some_arguments->count()==4 )
            {
               Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
               Type t2 = extract_arg( some_arguments, 2 )->data_type() ;
               Type t3 = extract_arg( some_arguments, 3 )->data_type() ;
               result = ( t1==Double && t2==Double && t3==Int ) ;
            }
         }
         break ;
      case DblEq :
         result = ( some_arguments->count()==4 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
            Type t2 = extract_arg( some_arguments, 2 )->data_type() ;
            Type t3 = extract_arg( some_arguments, 3 )->data_type() ;
            result = ( t0==Double && t1==Double && t2==Double && t3==Double ) ;
         }
         break ;
      case RandDbl :
      case RandInt :
         result = ( some_arguments->count()==0 ) ;
         break ;
      default :
         PEL_MathFunctionExp_ERROR::n0( "matches_args", name() ) ;
         break ;
   }
   
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_MathFunctionExp:: to_bool( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MathFunctionExp:: to_bool" ) ;
   PEL_CHECK_PRE( to_bool_PRE(ct) ) ;
   
   bool result = false ;

   switch( OP )
   {
      case DblEq :
         result = PEL::double_equality( ARG1->to_double( ct ),
                                        ARG2->to_double( ct ),
                                        ARG3->to_double( ct ),
                                        ARG4->to_double( ct ) ) ;
         break ;
      default :
         PEL_MathFunctionExp_ERROR::n0( "to_bool", name() ) ;
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
double
PEL_MathFunctionExp:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MathFunctionExp:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE(ct) ) ;
   
   double result = PEL::max_double() ;

   switch( OP )
   {
      case Sqr :
         result = PEL::sqr( ARG1->to_double(ct) ) ;
         break ;
      case Sqrt :
         result = PEL::sqrt( ARG1->to_double(ct) ) ;
         break ;
      case Pow :
         result = PEL::pow( ARG1->to_double(ct), ARG2->to_double(ct) ) ;
         break ;
      case Exp :
         result = PEL::exp( ARG1->to_double(ct) ) ;
         break ;
      case Log :
         result = PEL::log( ARG1->to_double(ct) ) ;
         break ;
      case Log10 :
         result = PEL::log10( ARG1->to_double(ct) ) ;
         break ;
      case Sin :
         result = PEL::sin( ARG1->to_double(ct) ) ;
         break ;
      case Cos :
         result = PEL::cos( ARG1->to_double(ct) ) ;
         break ;
      case Tan :
         result = PEL::tan( ARG1->to_double(ct) ) ;
         break ;
      case ASin :
         result = PEL::asin( ARG1->to_double(ct) ) ;
         break ;
      case ACos :
         result = PEL::acos( ARG1->to_double(ct) ) ;
         break ;
      case ATan :
         result = PEL::atan( ARG1->to_double(ct) ) ;
         break ;
      case ATan2 :
         result = PEL::atan2( ARG1->to_double(ct), ARG2->to_double(ct) ) ;
         break ;
      case Sinh :
         result = PEL::sinh( ARG1->to_double(ct) ) ;
         break ;
      case Cosh :
         result = PEL::cosh( ARG1->to_double(ct) ) ;
         break ;
      case Tanh :
         result = PEL::tanh( ARG1->to_double(ct) ) ;
         break ;
      case ASinh :
         result = PEL::asinh( ARG1->to_double(ct) ) ;
         break ;
      case ACosh :
         result = PEL::acosh( ARG1->to_double(ct) ) ;
         break ;
      case ATanh :
         result = PEL::atanh( ARG1->to_double(ct) ) ;
         break ;
      case J0 :
         result = PEL::j0( ARG1->to_double(ct) ) ;
         break ;
      case J1 :
         result = PEL::j1( ARG1->to_double(ct) ) ;
         break ;
      case Jn :
         result = PEL::jn( ARG1->to_int(ct), ARG2->to_double(ct) ) ;
         break ;
      case Y0 :
         result = PEL::y0( ARG1->to_double(ct) ) ;
         break ;
      case Y1 :
         result = PEL::y1( ARG1->to_double(ct) ) ;
         break ;
      case Yn :
         result = PEL::yn( ARG1->to_int(ct), ARG2->to_double(ct) ) ;
         break ;
      case Gamma :
         result = PEL::gamma( ARG1->to_double(ct) ) ;
         break ;
      case LGamma :
         result = PEL::lgamma( ARG1->to_double(ct) ) ;
         break ;
      case Erf :
         result = PEL::erf( ARG1->to_double(ct) ) ;
         break ;
      case Erfc :
         result = PEL::erfc( ARG1->to_double(ct) ) ;
         break ;
      case IGamma :
         {
            double const a = ARG1->to_double(ct) ;
            if( a <= 0 )
            {
               raise_error( "first argument should be stricktly positive" ) ;
            }
            double const x = ARG2->to_double(ct) ;
            if( x <= 0 )
            {
               raise_error( "second argument should be stricktly positive" ) ;
            }
            if( ARG3 == 0 )
            {
               result = PEL::incomplete_gamma( a, x ) ;
            }
            else
            {
               double const eps = ARG3->to_double(ct) ;
               if( eps <= 0 )
               {
                  raise_error( "third argument should be stricktly positive" ) ;
               }
               double const err = ARG4->to_double(ct) ;
               if( err <= 0 )
               {
                  raise_error( "fourth argument should be stricktly positive" ) ;
               }
               int const iter_max = ARG5->to_int(ct) ;
               if( iter_max <= 0 )
               {
                  raise_error( "fifth argument should be stricktly positive" ) ;
               }
               result = PEL::incomplete_gamma( a, x,
                                               eps, err, (size_t) iter_max ) ;
            }
         }
         break ;
      case En :
         {
            int const n = ARG1->to_int(ct) ;
            if( n < 0 )
            {
               raise_error( "first argument should be positive" ) ;
            }
            double const x = ARG2->to_double(ct) ;
            if( x < 0 )
            {
               raise_error( "second argument should be positive" ) ;
            }
            if( ARG3 == 0 )
            {
               if( x<1.E-30 && n<2 )
               {
                  raise_error( "first argument less than 2 and null second argument" ) ;
               }
               result = PEL::En( (size_t) n, x ) ;
            }
            else
            {
               double const eps = ARG3->to_double(ct) ;
               if( eps <= 0 )
               {
                  raise_error( "third argument should be stricktly positive" ) ;
               }
               double const err = ARG4->to_double(ct) ;
               if( err <= 0 )
               {
                  raise_error( "fourth argument should be stricktly positive" ) ;
               }
               int const iter_max = ARG5->to_int(ct) ;
               if( iter_max <= 0 )
               {
                  raise_error( "fifth argument should be stricktly positive" ) ;
               }
               if( x<eps && n<2 )
               {
                  raise_error( "first argument less than 2 and null second argument" ) ;
               }
               result = PEL::En( (size_t) n, x, eps, err, (size_t) iter_max ) ;
            }
         }
         break ;
      case Ei :
         {
            double const x = ARG1->to_double(ct) ;
            if( x < 0 )
            {
               raise_error( "first argument should be positive" ) ;
            }
            if( ARG2 == 0 )
            {
               result = PEL::Ei( x ) ;
            }
            else
            {
               double const eps = ARG2->to_double(ct) ;
               if( eps < 0 )
               {
                  raise_error( "second argument should be positive" ) ;
               }
               double const err = ARG3->to_double(ct) ;
               if( err <= 0 )
               {
                  raise_error( "third argument should be stricktly positive" ) ;
               }
               int const iter_max = ARG4->to_int(ct) ;
               if( iter_max <= 0 )
               {
                  raise_error( "fourth argument should be stricktly positive" ) ;
               }
               result = PEL::Ei( x, eps, err, (size_t) iter_max ) ;
            }
         }
         break ;
      case Abs :
         result = PEL::abs( ARG1->to_double(ct) ) ;
         break ;
      case Floor :
         result = PEL::floor( ARG1->to_double(ct) ) ;
         break ;
      case Ceil :
         result = PEL::ceil( ARG1->to_double(ct) ) ;
         break ;
      case  RandDbl :
         result = PEL::random_double() ;
         break ;
      default :
         PEL_MathFunctionExp_ERROR::n0( "to_double", name() ) ;
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
int
PEL_MathFunctionExp:: to_int( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MathFunctionExp:: to_int" ) ;
   PEL_CHECK_PRE( to_int_PRE(ct) ) ;
   
   int result = PEL::max_int() ;

   switch( OP )
   {
      case Abs :
         result = PEL::abs( ARG1->to_int(ct) ) ;
         break ;
      case RandInt :
         result = PEL::rand() ;
         break ;
      default :
         PEL_MathFunctionExp_ERROR::n0( "to_int", name() ) ;
         break ;
   }
   return result ;
}


//----------------------------------------------------------------------
PEL_Data*
PEL_MathFunctionExp:: create_derivative( PEL_Object* a_owner,
                                         PEL_Variable const* var,
                                         PEL_Context const* ct ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MathFunctionExp:: create_derivative" ) ;
   PEL_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;
   
   PEL_Data* result = 0 ;
   if( OP==Pow )
   {
      PEL_List * lst = PEL_List::create( 0 ) ;
      lst->append( ARG1->create_clone( lst ) ) ;
      PEL_Data* log = PEL_Expression::create( 0, "log", lst ) ;
      lst->set_owner( log ) ;
      
      lst = PEL_List::create( 0 ) ;
      lst->append( ARG2->create_derivative( lst, var, ct ) ) ;
      lst->append( log ) ; log->set_owner( lst ) ;
      PEL_Data* yprime_log = PEL_Expression::create( 0, "*", lst ) ;
      lst->set_owner( yprime_log ) ;
      
      lst = PEL_List::create( 0 ) ;
      lst->append( yprime_log ) ; yprime_log->set_owner( lst ) ;
      lst->append( create_clone( lst ) ) ;
      PEL_Data* yprime_log_x_pow_y = PEL_Expression::create( 0, "*", lst ) ;
      lst->set_owner( yprime_log_x_pow_y ) ;

      lst = PEL_List::create( 0 ) ;
      lst->append( ARG2->create_clone( lst ) ) ; 
      lst->append( ARG1->create_derivative( lst, var, ct ) ) ;
      PEL_Data* y_xprime = PEL_Expression::create( 0, "*", lst ) ;
      lst->set_owner( y_xprime ) ;
      
      lst = PEL_List::create( 0 ) ;
      lst->append( ARG2->create_clone( lst ) ) ; 
      lst->append( PEL_Double::create( lst, 1.0 ) ) ; 
      PEL_Data* ym1 = PEL_Expression::create( 0, "-", lst ) ;
      lst->set_owner( ym1 ) ;

      lst = PEL_List::create( 0 ) ;
      lst->append( ARG1->create_clone( lst ) ) ; 
      lst->append( ym1 ) ;  ym1->set_owner( lst ) ;
      PEL_Data* x_pow_ym1 = PEL_Expression::create( 0, "pow", lst ) ;
      lst->set_owner( x_pow_ym1 ) ;
      
      lst = PEL_List::create( 0 ) ;
      lst->append( y_xprime ) ; y_xprime->set_owner( lst ) ;
      lst->append( x_pow_ym1 ) ; x_pow_ym1->set_owner( lst ) ;
      PEL_Data* y_xprime_x_pow_ym1 = PEL_Expression::create( 0, "*", lst ) ;
      lst->set_owner( y_xprime_x_pow_ym1 ) ;
      
      lst = PEL_List::create( 0 ) ;
      lst->append( yprime_log_x_pow_y ) ;yprime_log_x_pow_y->set_owner( lst ) ;
      lst->append( y_xprime_x_pow_ym1 ) ;y_xprime_x_pow_ym1->set_owner( lst ) ;
      result = PEL_Expression::create( a_owner, "+", lst ) ;
      lst->set_owner( result ) ;
      
   }
   else
   {
      PEL_Data* composed_derivative = 0 ;
      PEL_Data const* der_variable = ARG1 ;
      PEL_List * lst = PEL_List::create( 0 ) ;
      switch( OP )
      {
         case Sqr :
            lst->append( PEL_Double::create( lst, 2.0 ) ) ;
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = PEL_Expression::create( 0, "*", lst ) ;
            lst->set_owner( composed_derivative ) ;
         
            break ;
         case Sqrt :
            lst->append( ARG1->create_clone( lst ) ) ;
            lst->append( PEL_Double::create( lst, -0.5 ) ) ;
            composed_derivative = PEL_Expression::create( 0, "pow", lst ) ;
            lst->set_owner( composed_derivative ) ;
            lst = PEL_List::create( 0 ) ;
            lst->append( PEL_Double::create( lst, 0.5 ) ) ;
            lst->append( composed_derivative ) ; composed_derivative->set_owner( lst ) ;
            composed_derivative = PEL_Expression::create( 0, "*", lst ) ;
            lst->set_owner( composed_derivative ) ;
         
            break ;
         case Exp :
            composed_derivative = create_clone( 0 ) ;
            lst->destroy() ;
            break ;
         case Log :
            lst->append( PEL_Double::create( lst, 1.0 ) ) ;
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = PEL_Expression::create( 0, "/", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case Sin :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = PEL_Expression::create( 0, "cos", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case Cos :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = PEL_Expression::create( 0, "sin", lst ) ;
            lst->set_owner( composed_derivative ) ;
            lst = PEL_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = PEL_Expression::create( 0, "unary_minus", lst ) ;
            lst->set_owner( composed_derivative ) ;         
            break ;
         case Tan :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = PEL_Expression::create( 0, "cos", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = PEL_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            lst->append( PEL_Double::create( lst, -2.0 ) ) ;
            composed_derivative = PEL_Expression::create( 0, "pow", lst ) ;
            lst->set_owner( composed_derivative ) ;         
            break ;
         case ASin :
         case ACos :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = PEL_Expression::create( 0, "sqr", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = PEL_List::create( 0 ) ;
            lst->append( PEL_Double::create( lst, 1.0 ) ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = PEL_Expression::create( 0, "-", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = PEL_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            lst->append( PEL_Double::create( lst, -0.5 ) ) ;
            composed_derivative = PEL_Expression::create( 0, "pow", lst ) ;
            lst->set_owner( composed_derivative ) ;
            if( OP==ASin ) break ;
            
            lst = PEL_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = PEL_Expression::create( 0, "unary_minus", lst ) ;
            lst->set_owner( composed_derivative ) ;         
            break ;
         case ATan :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = PEL_Expression::create( 0, "sqr", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = PEL_List::create( 0 ) ;
            lst->append( PEL_Double::create( lst, 1.0 ) ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = PEL_Expression::create( 0, "+", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = PEL_List::create( 0 ) ;
            lst->append( PEL_Double::create( lst, 1.0 ) ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = PEL_Expression::create( 0, "/", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case Sinh :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = PEL_Expression::create( 0, "cosh", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case Cosh :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = PEL_Expression::create( 0, "sinh", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case Tanh :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = PEL_Expression::create( 0, "cosh", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = PEL_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            lst->append( PEL_Double::create( lst, -2.0 ) ) ;
            composed_derivative = PEL_Expression::create( 0, "pow", lst ) ;
            lst->set_owner( composed_derivative ) ;         
            break ;
         case Erf :
         case Erfc :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = PEL_Expression::create( 0, "sqr", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = PEL_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = PEL_Expression::create( 0, "unary_minus", lst ) ;
            lst->set_owner( composed_derivative ) ;

            lst = PEL_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = PEL_Expression::create( 0, "exp", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = PEL_List::create( 0 ) ;
            lst->append(
               PEL_Double::create( lst,
                                   ( OP==Erf ? 1.0 : -1.0 )*2.0/PEL::sqrt( PEL::pi() ) ) ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = PEL_Expression::create( 0, "*", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case ASinh :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = PEL_Expression::create( 0, "sqr", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = PEL_List::create( 0 ) ;
            lst->append( PEL_Double::create( lst, 1.0 ) ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = PEL_Expression::create( 0, "+", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = PEL_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            lst->append( PEL_Double::create( lst, -0.5 ) ) ;
            composed_derivative = PEL_Expression::create( 0, "pow", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case ACosh :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = PEL_Expression::create( 0, "sqr", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = PEL_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            lst->append( PEL_Double::create( lst, 1.0 ) ) ;
            composed_derivative = PEL_Expression::create( 0, "-", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = PEL_List::create( 0 ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            lst->append( PEL_Double::create( lst, -0.5 ) ) ;
            composed_derivative = PEL_Expression::create( 0, "pow", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case ATanh :
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = PEL_Expression::create( 0, "sqr", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = PEL_List::create( 0 ) ;
            lst->append( PEL_Double::create( lst, 1.0 ) ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = PEL_Expression::create( 0, "-", lst ) ;
            lst->set_owner( composed_derivative ) ;
            
            lst = PEL_List::create( 0 ) ;
            lst->append( PEL_Double::create( lst, 1.0 ) ) ;
            lst->append( composed_derivative ) ;composed_derivative->set_owner(lst) ;
            composed_derivative = PEL_Expression::create( 0, "/", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case Ei :
            // e(x)/x
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = PEL_Expression::create( 0, "exp", lst ) ;
            lst->set_owner( composed_derivative ) ;

            lst = PEL_List::create( 0 ) ;
            lst->append( composed_derivative ) ;
            composed_derivative->set_owner( lst ) ;
            lst->append( ARG1->create_clone( lst ) ) ;
            composed_derivative = PEL_Expression::create( 0, "/", lst ) ;
            lst->set_owner( composed_derivative ) ;
            break ;
         case En :
            // n==0 ? -(x+1)*e(-x)/x^2 : -E(n-1,x)
            {
               der_variable = ARG2 ;
               int const n = ARG1->to_int( ct ) ;
               if( n == 0 )
               {
                  lst->append( ARG2->create_clone( lst ) ) ;
                  composed_derivative =
                     PEL_Expression::create( 0, "unary_minus", lst ) ;
                  lst->set_owner( composed_derivative ) ;

                  lst = PEL_List::create( 0 ) ;
                  lst->append( composed_derivative ) ;
                  composed_derivative->set_owner( lst ) ;
                  composed_derivative =
                     PEL_Expression::create( 0, "exp", lst ) ;
                  lst->set_owner( composed_derivative ) ;

                  lst = PEL_List::create( 0 ) ;
                  lst->append( PEL_Double::create( lst, -1. ) ) ;
                  lst->append( ARG2->create_clone( lst ) ) ;
                  PEL_Expression* exp = 
                     PEL_Expression::create( 0, "-", lst ) ;
                  lst->set_owner( exp ) ;

                  lst = PEL_List::create( 0 ) ;
                  lst->append( exp ) ;
                  exp->set_owner( lst ) ;
                  lst->append( composed_derivative ) ;
                  composed_derivative->set_owner( lst ) ;
                  composed_derivative =
                     PEL_Expression::create( 0, "*", lst ) ;
                  lst->set_owner( composed_derivative ) ;

                  lst = PEL_List::create( 0 ) ;
                  lst->append( composed_derivative ) ;
                  composed_derivative->set_owner( lst ) ;
                  lst->append( ARG2->create_clone( lst ) ) ;
                  composed_derivative =
                     PEL_Expression::create( 0, "/", lst ) ;
                  lst->set_owner( composed_derivative ) ;
                  
                  lst = PEL_List::create( 0 ) ;
                  lst->append( composed_derivative ) ;
                  composed_derivative->set_owner( lst ) ;
                  lst->append( ARG2->create_clone( lst ) ) ;
                  composed_derivative =
                     PEL_Expression::create( 0, "/", lst ) ;
                  lst->set_owner( composed_derivative ) ;
               }
               else
               {
                  lst->append( PEL_Int::create( lst, n-1 ) ) ;
                  lst->append( ARG2->create_clone( lst ) ) ;
                  if( ARG3 != 0 )
                  {
                     lst->append( ARG3->create_clone( lst ) ) ;
                     lst->append( ARG4->create_clone( lst ) ) ;
                     lst->append( ARG5->create_clone( lst ) ) ;
                  }
                  composed_derivative = PEL_Expression::create( 0, "En", lst ) ;
                  lst->set_owner( composed_derivative ) ;

                  lst = PEL_List::create( 0 ) ;
                  lst->append( composed_derivative ) ;
                  composed_derivative->set_owner( lst ) ;
                  composed_derivative =
                     PEL_Expression::create( 0, "unary_minus", lst ) ;
                  lst->set_owner( composed_derivative ) ;
               }
            }
            break ;
         case J0 :
         case J1 :
         case Jn :
         case Y0 :
         case Y1 :
         case Yn :
         case Gamma :
         case LGamma :
         case ATan2 :
         case Log10 :
         case IGamma :
         default :
            PEL_MathFunctionExp_ERROR::n0( "create_derivative", name() ) ;
            break ;
      }
      PEL_Data* der = der_variable->create_derivative( 0, var, ct ) ;
      lst = PEL_List::create( 0 ) ;
      lst->append( der ) ; der->set_owner( lst ) ;
      lst->append( composed_derivative ) ; composed_derivative->set_owner( lst ) ;
      result = PEL_Expression::create( a_owner, "*", lst ) ;
      lst->set_owner( result ) ;
   }   
   PEL_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//internal--------------------------------------------------------------
void 
PEL_MathFunctionExp_ERROR:: n0( std::string const& f_name,
                                std::string const& op_name )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** PEL_MathFunctionExp::" + f_name +"\n" ;
   mesg += "    operation " + op_name + " not implemented." ;
   PEL_Error::object()->raise_internal( mesg ) ;
}
