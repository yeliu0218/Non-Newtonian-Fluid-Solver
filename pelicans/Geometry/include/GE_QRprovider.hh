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

#ifndef GE_QR_PROVIDER_HH
#define GE_QR_PROVIDER_HH

#include <PEL_Object.hh>

class PEL_Sequence ;

class GE_Customized_QR ;
class GE_QuadratureRule ;
class GE_Mpolyhedron ;
class GE_ReferencePolyhedron ;
class GE_SimplePolygon2D ;

class PEL_ObjectRegister ;
class PEL_Vector ;

/*
Providers of instances of `GE_QuadratureRule::' subclasses.

A one-to-one relationship between `GE_QuadratureRule::' and
`GE_ReferencePolyhedron::' objects.

IMPLEMENTATION :
   - Each concrete subclass has a unique instance.
   - A register is maintained that maps the name of each concrete subclass
     with the unique instance of this subclass.

FRAMEWORK INSTANTIATION :
   1. Derive a concrete subclass, say MyPro (derivation of abstract
      subclasses is improbable although possible).
   2. Choose a name for MyPro, say "my_pro".
   3. Declare all constructors private.
   4. Define the unique instance :
      4.1 Implement a default constructor that initializes the 
          `GE_QRprovider::' subobject by calling
                `GE_QRprovider( std::string )'
          with "my_pro" as argument value.
          Example of pseudo-code :
          | MyPro:: MyPro( void ) : GE_QRprovider( "my_pro" ) {}
      4.2 Define and initialize a private static instance
          by calling the default private constructor :
             declaration (in the header file, eg MyPro.hh) :
             |  static MyPro* REGISTRATOR ;
             definition (in the implementation file, eg MyPro.cc) :
             |  MyPro* MyPro::REGISTRATOR = new MyPro() ;
   5. Implement a private destructor, that sets to zero REGISTRATOR
      defined above.
   6. Implement the `build' method by successively calling  `add_QR' 
      for each desired association between a `GE_Mpolyhedron::' object and a 
      `GE_QuadratureRule::' object.
*/

class PEL_EXPORT GE_QRprovider : public PEL_Object
{
   public: //------------------------------------------------------------------

   //-- Instance delivery and initialization

      // the unique instance of the concrete subclass called `a_name'
      static GE_QRprovider const* object( std::string a_name ) ;

   //-- Internal status

      // indicator that changes only if some instance returned by
      // `::quadrature_rule' has changed
      size_t status( void ) const ;

   //-- Access to GE_QuadratureRule objects

      // instance of a subclass of `GE_QuadratureRule::' associated
      // to `poly'
      GE_QuadratureRule const* quadrature_rule(
                              GE_ReferencePolyhedron const* poly ) const ;

      // Create and return a quadrature rule defined on the intersection 
      // of `pp' and the polyhedra of `poly_list' for which integration
      // points are given by `self'.
      virtual GE_Customized_QR const* create_subset(
                              PEL_Object* a_owner,
                              GE_Mpolyhedron const* pp,
                              PEL_Sequence const* poly_list ) const ;

   protected: //---------------------------------------------------------------

   //-- Plug in

      virtual ~GE_QRprovider( void ) ;

      // only for the initialization of the unique instance
      // of the concrete class name of name `a_name'
      GE_QRprovider( std::string const& a_name ) ;

      // exclusively used by `GE_CustomizedQR_provider::'
      GE_QRprovider( PEL_Object* a_owner, std::string const& a_name ) ;

      // Build the one-to-one relationship between `GE_QuadratureRule::' and
      // `GE_ReferencePolyhedron::' objects, using one call to `::add_QR'
      // for each correspondance.
      virtual void build( void ) = 0 ;
      
      // Notify that `qr' is associated to `poly'.
      void add_QR( GE_ReferencePolyhedron const* poly,
                   GE_QuadratureRule const* qr ) ;

      // Change inner status.
      void notify_status_change( void ) const ;
      
   //--  Preconditions, Postconditions, Invariant

      virtual bool create_subset_PRE( PEL_Object const* a_owner,
                                      GE_Mpolyhedron const* pp,
                                      PEL_Sequence const* poly_list ) const ;
      
      virtual bool create_subset_POST( GE_Customized_QR const* result,
                                       PEL_Object const* a_owner,
                                       GE_Mpolyhedron const* pp,
                                       PEL_Sequence const* poly_list ) const ;

      virtual bool invariant( void ) const ;

   private: //-----------------------------------------------------------------

      GE_QRprovider( void ) ;
      GE_QRprovider( GE_QRprovider const& other ) ;
      GE_QRprovider& operator=( GE_QRprovider const& other ) ;

      static PEL_ObjectRegister* plugins_map( void ) ;
      
   //-- Attributes

      std::string const NAME ;
      PEL_Vector* const QRS ;
      mutable size_t STAT ;

} ;

#endif
