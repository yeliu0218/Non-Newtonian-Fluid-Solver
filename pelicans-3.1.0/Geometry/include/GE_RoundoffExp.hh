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

#ifndef GE_ROUNDOFF_EXP_HH
#define GE_ROUNDOFF_EXP_HH

#include <PEL_Expression.hh>

/*
Expression for rounding off a double, used by GE_Meshing derivatives

name      : default_roundoff
arguments : Double, Int, Double
type      : Double

default_roundoff( x, n, xmin ) is the round off of x computed by keeping
n significant digits if x>xmin, 0 otherwise

PUBLISHED
*/

class PEL_EXPORT GE_RoundoffExp : public PEL_Expression
{
   public: //--------------------------------------------------------
      
   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual double to_double( PEL_Context const* ct = 0 ) const ;
      
   protected: //-----------------------------------------------------
      
   private: //-------------------------------------------------------

      GE_RoundoffExp( void ) ;
     ~GE_RoundoffExp( void ) ;
      GE_RoundoffExp( GE_RoundoffExp const& other ) ;
      GE_RoundoffExp& operator=( GE_RoundoffExp const& other ) ;

      GE_RoundoffExp( PEL_Object* a_owner,
                      std::string const& a_name,
                      PEL_Sequence const* argument_list ) ;

   //-- Plug in

      GE_RoundoffExp( std::string const& a_name ) ;

      virtual PEL_Expression* create_replica( 
                                   PEL_Object * a_owner,
                                   PEL_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( PEL_Sequence const* some_arguments ) const ;
      
   //-- Class attributes

      static GE_RoundoffExp const* PROTOTYPE ;
      
   //-- Attributes

} ;

#endif
