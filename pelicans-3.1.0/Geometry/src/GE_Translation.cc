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

#include <GE_Translation.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>
#include <PEL.hh>

#include <GE_Color.hh>
#include <GE_Point.hh>
#include <GE_Vector.hh>

#include <iostream>
#include <sstream>

using std::ostringstream ;
using std::endl ;

//----------------------------------------------------------------------
GE_Translation*
GE_Translation:: create( PEL_Object* a_owner,
                         size_t nb_sp_dims,
                         PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Translation:: create" ) ;
   
   GE_Translation* result = new GE_Translation( a_owner, nb_sp_dims, exp ) ;
   
   result->build_inverse( a_owner ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Translation:: GE_Translation( PEL_Object* a_owner,
                                 size_t nb_sp_dims,
                                 PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : GE_Transform( a_owner, nb_sp_dims, exp )
   , VEC( GE_Vector::create( this, 
                             exp->doubleVector_data( "translation_vector" ) ) )
   , INVERSE( 0 )
{
   if( VEC->nb_components() != nb_sp_dims )
   {
      ostringstream mesg ;
      mesg << "the entry of keyword:" << endl ;
      mesg << "   \"translation\"" << endl ;
      mesg << "should be of length: " << nb_sp_dims ;
      PEL_Error::object()->raise_module_error( exp, mesg.str() ) ;
   }
}

//----------------------------------------------------------------------
GE_Translation:: GE_Translation( PEL_Object* a_owner )
//----------------------------------------------------------------------
   : GE_Transform( a_owner )
   , VEC( 0 )
{   
}

//----------------------------------------------------------------------
GE_Translation:: ~GE_Translation( void  )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
GE_Translation const*
GE_Translation:: inverse( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Translation:: inverse" ) ;
   
   GE_Translation const* result = INVERSE ;
   
   PEL_CHECK_POST( inverse_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
GE_Translation:: apply( GE_Point* pt ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Translation:: apply" ) ;
   PEL_CHECK_PRE( apply_PRE( pt ) ) ;

   pt->translate( 1.0, VEC ) ;
}

//----------------------------------------------------------------------
void
GE_Translation:: print( std::ostream& os, size_t indent_width ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Translation:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   os << space << "translation: " << endl ;
   VEC->print( os, indent_width+3 ) ; os << std::endl ;
   os << space << "   source color: \"" << source_color()->name() 
               << "\"" << endl ;
   os << space << "   target color: \"" << target_color()->name() 
               << "\"" << endl ;
}

//----------------------------------------------------------------------
void
GE_Translation:: build_inverse( PEL_Object* a_owner )
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_Translation:: build_inverse" ) ;
   
   GE_Translation* tr = new GE_Translation( a_owner ) ;
   tr->initialize_as_inverse( this ) ;
   tr->VEC = VEC->create_clone( tr ) ;
   tr->VEC->scale( -1.0 ) ;
   
   INVERSE = tr ;
   tr->INVERSE = this ; 
}
