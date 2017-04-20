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

#ifndef RS_AVULA_EXP_HH
#define RS_AVULA_EXP_HH

#include <PEL_Expression.hh>

#include <doubleVector.hh>

/*
Expressions related to the "Starting Flow in a Long Circular Tube I" benchmark

Reference:
   document: Implementation of Solvers for Partial Differential Equations 
             with PELICANS
   chapter:  Reference Solution of Some Systems of 
             Partial Differential Equations
   section:  Starting Flow in a Long Circular Tube I
   
PUBLISHED
*/

class PEL_EXPORT RS_AvulaEXP : public PEL_Expression
{
   public: //---------------------------------------------------------------
   
   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual doubleVector const& to_double_vector( 
                                          PEL_Context const* ct = 0 ) const ;

   protected: //-----------------------------------------------------------
      
   private: //-------------------------------------------------------------
    
      RS_AvulaEXP( void ) ;
     ~RS_AvulaEXP( void ) ;
      RS_AvulaEXP( RS_AvulaEXP const& other ) ;
      RS_AvulaEXP& operator=( RS_AvulaEXP const& other ) ;

      enum Func { U_AXI, P_AXI, U_3D, P_3D } ;
            
      RS_AvulaEXP( PEL_Object* a_owner,
                   std::string const& a_name,
                   PEL_Sequence const* argument_list,
                   Func an_expr ) ;

   //-- Plug in

      RS_AvulaEXP( std::string const& a_name, Func an_expr ) ;

      virtual RS_AvulaEXP* create_replica( 
                                   PEL_Object * a_owner,
                                   PEL_Sequence const* argument_list ) const ;
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( PEL_Sequence const* some_arguments ) const ;

   //-- Internals

      double uz_adim( double eta, double tau, double tau_cut, 
                      double a, double c, double reynolds, double tol ) const ;

      double hconvol( double tau, double a, double c, double gamma_n ) const ;

      double root_of_J0( size_t n ) const  ;

   //-- Class attributes

      static RS_AvulaEXP const* PROTOTYPE_U_AXI ; 
      static RS_AvulaEXP const* PROTOTYPE_P_AXI ;
      static RS_AvulaEXP const* PROTOTYPE_U_3D ;
      static RS_AvulaEXP const* PROTOTYPE_P_3D ;
     
   //-- Attributes

      Func EXPR ;
      mutable doubleVector DV_result_1 ;
      mutable doubleVector DV_result_2 ;
      mutable doubleVector DV_result_3 ;
} ;

#endif
