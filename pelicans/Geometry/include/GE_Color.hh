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

#ifndef GE_COLOR_HH
#define GE_COLOR_HH

#include <PEL_Object.hh>

#include <string>
#include <stringVector.hh>
#include <intArray2D.hh>

class PEL_ListIdentity ;
class PEL_ObjectRegister ;

/*
Markers used to identify clusters of geometrical entities.

PUBLISHED*/

class PEL_EXPORT GE_Color : public PEL_Object
{
   public: //------------------------------------------------------------

   //-- Instance delivery and initialization

      // the instance called `a_name' if any (a fatal error is raised if none)
      static GE_Color const* object( std::string const& a_name ) ;

      // the "null color" instance.
      static GE_Color const* null_color( void ) ;

      // the "halo color" instance.
      static GE_Color const* halo_color( void ) ;

   //-- Management of the set of instances(0.2)

      // Does a color of name `a_name' exist ? 
      static bool exist( std::string const& a_name ) ;

      // Create leaf instance called `a_name' if it does not exist.
      static void extend( std::string const& a_name ) ;

      // Create a composite instance called `a_name' whose parts 
      // are identified by the name container `a_name_list' if it does not 
      // exist. If any of the items of `a_name_list' is not the name of
      // an existing color, a fatal error is raised.
      static void extend( std::string const& a_name,
                          stringVector const& a_name_list ) ;

      // vector of the names of all instances
      static stringVector const& color_table( void ) ;

      // connectivity between macro colors and single ones
      static intArray2D const& color_table_connectivity( void ) ;
      
   //-- Identification

      // identifier uniquely identifying `self'
      int identifier( void) const ;

      // name of `self' 
      std::string const& name( void ) const ;

   //-- Characteristics

      // Is `self' a composition of other instances ? 
      bool is_composite( void ) const ;

      // Does `self' contains a color of name `a_name' in its composition ? 
      bool has( std::string const& a_name ) const ;

   //-- Comparison

      // Are `self' and `other' sharing some color ? 
      bool is_overlapping( GE_Color const* other ) const ;

      // Are `self' and `other' matching up ? 
      bool is_matching( GE_Color const* other ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //---------------------------------------------------------

   private: //-----------------------------------------------------------

      GE_Color( void ) ;
     ~GE_Color( void ) ;
      GE_Color( GE_Color const& other ) ;
      GE_Color& operator=( GE_Color const& other ) ;

      GE_Color( PEL_Object* a_owner, std::string const& a_name ) ;

      // map whose items are the instances and whose keys are their name
      static PEL_ObjectRegister* colors( void ) ;

   //-- Class attributes

      static int NEXT_ID ;
      static stringVector NAMES ;
      static intArray2D CONNECTIVITY ;
      
   //-- Attributes

      std::string MY_NAME ;
      int ID ;
      bool IS_COMPOSITE ;
      PEL_ListIdentity* COMPOSING_COLORS ;
      
} ;

#endif
