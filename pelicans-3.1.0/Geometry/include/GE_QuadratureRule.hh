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

#ifndef GE_QUADRATURE_RULE_HH
#define GE_QUADRATURE_RULE_HH

#include <PEL_Object.hh>

#include <doubleVector.hh>

class PEL_ObjectRegister ;
class PEL_Vector ;

class GE_Point ;
class GE_ReferencePolyhedron ;


/*
Sets of points and sets of weights that define a method of numerical
integration, in the Gaussian sense, over a polyhedron.

Each quadrature rule can be characterized by its order. 
Integrals of any function over the polyhedron are approximated by the sum of 
function values at the quadrature points multiplied by the quadrature weights.
This method is exact for polynomials of degree up to the rule order.

Most rules are "instancied" once and for all by static method object ( unique 
instance ) but blank rules of type GE_Customized_QR can also been 
instancied.

IMPLEMENTATION :
   A register is maintained that maps the name of each concrete subclass
   (apart from `GE_Customized_QR::') with the unique instance of 
   this subclass.

FRAMEWORK INSTANTIATION :
   1. Derive a concrete subclass, say `MyRule'.
   2. Choose a name for `MyRule', say "MyRule".
   3. Declare all constructors private.
   4. Define the unique instance :
      4.1 Implement a default constuctor that initializes the 
          `GE_QuadratureRule::' subobject by calling
     `GE_QuadratureRule( std::string, GE_ReferencePolyhedron const*, size_t )'
          with "MyRule" as a first argument.
          The default constructor should call in sequence
             `append_point( GE_Point*, double )'
          for all integration points, and finally
             `set_sum_of_weights( double )'
      4.2 Define and initialize a private static instance of `MyRule' by
          calling the default private constructor :
            declaration :
               static MyRule* REGISTRATOR ;
            definition :
               MyRule* MyRule::REGISTRATOR = new MyRule() ;
   5. Implement a private destructor, that sets to zero REGISTRATOR
      defined above.
*/

class PEL_EXPORT GE_QuadratureRule : public PEL_Object
{
   public: //------------------------------------------------------------------

   //-- Instance delivery and initialization

      // a descendant instance
      static GE_QuadratureRule const* object( std::string a_name ) ;
      
   //-- Identification

      // registration name of `self'
      std::string const& name( void ) const ;

   //-- Rule characteristics

      // polyhedron associated to `self'
      GE_ReferencePolyhedron const* reference_polyhedron( void ) const  ;

      // if non-zero, maximum degree of the polynomial for which the method
      // of numerical integration defined by `self' is exact
      size_t order( void ) const ;

      // number of integration points
      size_t nb_points( void ) const ;

      // geometrical location of the `i'-th integration point
      GE_Point const* point( size_t i ) const ;

      // weight of the `i'-th integration point
      double weight( size_t i ) const ;

      // sum of the weight of all integration points
      double sum_of_weights( void ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;


   protected: //---------------------------------------------------------------

      // Initialize a `GE_QuadratureRule::' subobject of a complete
      // object whose generating class is called `a_name',
      // whose order is `a_order', associated to the polyhedron `poly'.
      GE_QuadratureRule( std::string a_name,
                         GE_ReferencePolyhedron const* poly,
                         size_t a_order )  ;

      // exclusively for `GE_Customized_QR::'
      GE_QuadratureRule( PEL_Object* a_owner,
                         GE_ReferencePolyhedron const* poly,
                         size_t a_order ) ;

      virtual ~GE_QuadratureRule( void ) ;

      // Make `pt' an element of the set of integration points with an
      // associated quadrature weight `pt_weight' (to be called in the
      // default constructor of concrete descendants).
      // On exit, `self' still refers to the object identified by `pt'.
      void append_point( GE_Point* pt, double pt_weight ) ;

      // Set the value of the sum of the weight of all integration points 
      // at `sum' (that value will be subsequently returned 
      // by `::sum_of_weights', which avoids useless calculations).
      void set_sum_of_weights( double sum ) ;

      // Clear all points.
      void clear_points( void ) ;
      
   //--  Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   private: //-----------------------------------------------------------------

      GE_QuadratureRule( void ) ;
      GE_QuadratureRule( GE_QuadratureRule const& other ) ;
      GE_QuadratureRule const& operator=( GE_QuadratureRule const& other ) ;

      static PEL_ObjectRegister* plugins_map( void ) ;
 
   //--  Preconditions, Postconditions, Invariant

      bool set_sum_of_weights_PRE( double sum ) const ;
      bool sum_of_weights_POST( double result ) const ;
         
   //-- Attributes
   
      std::string const RULE_NAME ;
      GE_ReferencePolyhedron const* const POLY ;
      size_t const ORDER ;
      PEL_Vector* const POINTS ;
      doubleVector WEIGHTS ;
      double TOT_WEIGHT ;
} ;

#endif
