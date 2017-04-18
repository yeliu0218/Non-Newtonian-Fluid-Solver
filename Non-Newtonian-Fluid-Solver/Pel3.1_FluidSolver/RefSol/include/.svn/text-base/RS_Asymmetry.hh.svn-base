/*
 *  Copyright : 
 *    "Institut de Radioprotection et de S�ret� Nucl�aire - IRSN" (1995-2008)
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

#ifndef RS_Asymmetry_HH
#define RS_Asymmetry_HH

#include <PEL_Expression.hh>

#include <doubleVector.hh> 
#include <doubleArray2D.hh>

/*
Expressions related to the "Variable Density Incompressible Flows I"
benchmark

Reference:
   document: Implementation of Solvers for Partial Differential Equations 
             with PELICANS
   chapter:  Reference Solution of Some Systems of 
             Partial Differential Equations
   section:  Variable Density Incompressible Flows I

PUBLISHED
*/

class RS_Asymmetry : public PEL_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual doubleVector const& to_double_vector( 
	                                  PEL_Context const* ct = 0 ) const ; 
      
      
      
   protected: //-------------------------------------------------------
      
   private: //-------------------------------------------------------

      RS_Asymmetry( void ) ;
     ~RS_Asymmetry( void ) ;
      RS_Asymmetry( RS_Asymmetry const& other ) ;
      RS_Asymmetry& operator=( 
                                  RS_Asymmetry const& other ) ;

      enum Func { U, CC} ;
            
      RS_Asymmetry( PEL_Object* a_owner,
		                  std::string const& a_name,
    		                  PEL_Sequence const* argument_list,
		                  Func an_exp ) ;

      PEL_Data const* alternative_result( void ) const ;
      
   //-- Plug in

      RS_Asymmetry( std::string const& a_name, Func an_exp ) ;

      virtual RS_Asymmetry* create_replica( 
                                   PEL_Object * a_owner,
                                   PEL_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( PEL_Sequence const* some_arguments ) const ;
      
   //-- Class attributes

      static RS_Asymmetry const* PROTOTYPE_U ;
      static RS_Asymmetry const* PROTOTYPE_CC ;
    
      
   //-- Attributes
      
      Func const EXPR ;
      mutable doubleVector DV_result_1 ;
      mutable doubleVector DV_result_2 ;
    
} ;

#endif
