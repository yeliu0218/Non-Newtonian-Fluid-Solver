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

#ifndef PEL_CONVERT_TYPE_EXP_HH
#define PEL_CONVERT_TYPE_EXP_HH

#include <PEL_Expression.hh>

/* 
Expressions for type conversion

--
name     : double
argument : Int
type     : Double

double( i ) converts the integer i into a double using the standard
C conversion

example :
   double( 1 ) : value is 1.0

--
name     : int
argument : Double
type     : Int

int( x ) converts the double x into an int

example :
   int( 1.0 )   : value is 1
   int( 1.01 )  : value is 1
   int( 1.99 )  : value is 1
   int( -2.01 ) : value is -2
   int( -2.9 )  : value is -2
--

PUBLISHED
*/

class PEL_EXPORT PEL_ConvertTypeExp : public PEL_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual double to_double( PEL_Context const* ct = 0 ) const ;

      virtual int to_int( PEL_Context const* ct = 0 ) const ;
      
   protected: //-------------------------------------------------------
      
   private: //-------------------------------------------------------

      PEL_ConvertTypeExp( void ) ;
     ~PEL_ConvertTypeExp( void ) ;
      PEL_ConvertTypeExp( PEL_ConvertTypeExp const& other ) ;
      PEL_ConvertTypeExp& operator=( PEL_ConvertTypeExp const& other ) ;

      PEL_ConvertTypeExp( PEL_Object* a_owner,
                          std::string const& a_name,
                          PEL_Sequence const* argument_list ) ;

   //-- Plug in

      PEL_ConvertTypeExp( std::string const& a_name ) ;

      virtual PEL_ConvertTypeExp* create_replica( 
                                   PEL_Object * a_owner,
                                   PEL_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( PEL_Sequence const* some_arguments ) const ;
      
   //-- Class attributes

      
      static PEL_ConvertTypeExp const* PROTOTYPE_DOUBLE ;
      static PEL_ConvertTypeExp const* PROTOTYPE_INT ;
      
   //-- Attributes

      PEL_Data const* const P0 ;
} ;

#endif
