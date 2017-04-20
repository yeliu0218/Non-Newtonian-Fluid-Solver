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


#include <iostream>
#include <map>

#include <PEL_export.hh>

class MI_OneBCbuilder ;

// Sand http adress is http://sand .
// Maintenance pelicans mail adress is mailto:maintenance.pelican@irsn.fr .

/*
  ...

  // error Fault : 00253
  Example :
     MODULE toto
       titi = 3. // titi is a dimensionless parameter
     END MODULE toto

*/

class PEL_EXPORT DOC_SyntaxTest
{
   public://-----------------------------------------------
      
   //-- Input - Output

       friend std::ostream& operator<<( std::ostream& out, DOC_SyntaxTest const& ex ) ;
      
   //-- Stuct declaration

       // Intern struct declaration.
       struct index { size_t i ; size_t l ; } ;

       // Reference to intern struct.
       static void fool( DOC_SyntaxTest::index idx ) ;
      
   //-- Polymorphism on pre and post conditions

      // @SIGNET fint
      void f( int i ) ;
      // @SIGNET fdouble
      void f( double d ) ;

      // First way to solve problem.
      bool fint_PRE( int i ) const ;
      bool fdouble_PRE( double d ) const ;
      // Second way to solve problem.
      bool f_POST( int i ) const ;
      bool f_POST( double d ) const ;
      
   private://-----------------------------------------------

      std::map< std::string, MI_OneBCbuilder* >::const_iterator IT_BC_to_build ;
};
