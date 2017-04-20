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

#ifndef PEL_Randomizer_HH
#define PEL_Randomizer_HH

#include <PEL_Object.hh>
#include <doubleVector.hh>

class size_t_vector ;

// Random number generator.
//     LONG PERIOD (>2.D18) RANDOM NUMBER GENERATOR OF L'ECUYER
//     WITH BAYS-DURHAM SHUFFLE AND ADDED SAFEGUARDS. RETURNS A
//     UNIFORM RANDOM DEVIATE BETWEEN 0.0 AND 1.0. CALL WITH IDUM
//     A NEGATIVE INTEGER TO INITIALIZE. THERAFTER, DO NOT ALTER
//     IDUM BETWEEN SUCCESSIVE DEVIATES IN A SEQUENCE. RNMX SHOULD
//     APPROXIMATE THE LARGEST FLOATING VALUE THAT IS LESS THAN 1.

class PEL_EXPORT PEL_Randomizer : public PEL_Object
{
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static PEL_Randomizer* create( PEL_Object* a_owner,
                                    int series ) ;

      virtual PEL_Randomizer* create_clone( PEL_Object* a_owner ) const ;

  //-- Cursor movement

      // random value
      virtual double item( void ) const ;

      // Start random iterator.
      virtual void start( void ) ;
  
      // Go to next value.
      virtual void go_next( void ) ;

      void build_permutation( size_t_vector& vec ) ;
      
   protected: //------------------------------------------------------------

      virtual ~PEL_Randomizer( void ) ;

      PEL_Randomizer( PEL_Object* a_owner, int series ) ;
      
   private: //--------------------------------------------------------------

      PEL_Randomizer( void ) ;
      PEL_Randomizer( PEL_Randomizer const& other ) ;
      PEL_Randomizer const& operator=( PEL_Randomizer const& other ) ;

      int IDUM ;
      int IDUM2 ;
      int IY ;
      int my_series ;
      doubleVector IV ;
      double value ;
      
      
} ;

#endif
