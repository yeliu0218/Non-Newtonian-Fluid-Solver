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

#include <PDE_SetOfDiscreteFields.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ObjectReader.hh>
#include <PEL_ObjectWriter.hh>
#include <PEL_Vector.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_ReferenceElement.hh>

#include <iostream>
#include <sstream>

using std::ostringstream ; using std::endl ;

//----------------------------------------------------------------------
PDE_SetOfDiscreteFields*
PDE_SetOfDiscreteFields:: create( PEL_Object* a_owner )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDiscreteFields:: create" ) ;

   PDE_SetOfDiscreteFields* result = new PDE_SetOfDiscreteFields( a_owner ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_fields() == 0 ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_SetOfDiscreteFields:: PDE_SetOfDiscreteFields( PEL_Object* a_owner )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , FIELDS( PEL_Vector::create( this, 0 ) )
   , ELMS( PEL_Vector::create( this, 0 ) )
   , ON_BOUNDS( 0 )
   , INDEX( PEL::bad_index() )
{
}

//----------------------------------------------------------------------
PDE_SetOfDiscreteFields:: ~PDE_SetOfDiscreteFields( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
size_t
PDE_SetOfDiscreteFields:: nb_fields( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDiscreteFields:: nb_fields" ) ;
   return( FIELDS->count() ) ;
}

//----------------------------------------------------------------------
PDE_DiscreteField*
PDE_SetOfDiscreteFields:: item( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDiscreteFields:: item" ) ;
   PEL_CHECK_PRE( i < nb_fields() ) ;

   PDE_DiscreteField* result = 
                      static_cast<PDE_DiscreteField*>( FIELDS->at( i ) ) ;
                      
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_SetOfDiscreteFields:: boundary_located( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDiscreteFields:: boundary_located" ) ;
   PEL_CHECK_PRE( i < nb_fields() ) ;
   
   return( ON_BOUNDS( i ) ) ;
}

//----------------------------------------------------------------------
PEL_Vector const*
PDE_SetOfDiscreteFields:: reference_elements( size_t i ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDiscreteFields:: reference_elements" ) ;
   PEL_CHECK_PRE( i < nb_fields() ) ;
   
   PEL_Vector const* result = static_cast< PEL_Vector* >( ELMS->at( i ) ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->index_limit() >= 1 ) ;
   PEL_CHECK_POST( result->count() == result->index_limit() ) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t j=0 ; j<result->index_limit() ; ++j ),
              dynamic_cast< PDE_ReferenceElement const* >( result->at( j ) ) 
              != 0  ) ) ;
   return( result ) ; 
}

//----------------------------------------------------------------------
bool
PDE_SetOfDiscreteFields:: has( std::string const& field_name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDiscreteFields:: has" ) ;
   PEL_CHECK_PRE( !field_name.empty() ) ;

   return( find_index( field_name )<nb_fields() ) ;
}

//----------------------------------------------------------------------
PDE_DiscreteField*
PDE_SetOfDiscreteFields:: item( std::string const& field_name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDiscreteFields:: item( std::string const& )" ) ;
   PEL_CHECK_PRE( !field_name.empty() ) ;

   size_t const idx = find_index( field_name ) ;
   if( idx == PEL::bad_index() )
   {
      ostringstream mesg ;
      mesg << "PDE_SetOfDiscreteFields :" << endl ;
      mesg << "   there is no recorded" << endl ;
      mesg << "   instance of PDE_DiscreteField" << endl ;
      mesg << "   called \"" << field_name << "\"" ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }
 
   PDE_DiscreteField* result =
               static_cast<PDE_DiscreteField*>( FIELDS->at( idx ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->name() == field_name ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_SetOfDiscreteFields:: index_of( std::string const& field_name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDiscreteFields:: index_of" ) ;
   PEL_CHECK_PRE( !field_name.empty() ) ;

   size_t result = find_index( field_name ) ;

   PEL_CHECK_POST( ( result < nb_fields() ) || 
                   ( result == PEL::bad_index() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_SetOfDiscreteFields:: append( PDE_DiscreteField* field,
                                  bool located_on_bounds,
                                  PEL_Vector* elements )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDiscreteFields:: append" ) ;
   PEL_CHECK_PRE( field != 0 ) ;
   PEL_CHECK_PRE( field->owner()==0 ) ;
   PEL_CHECK_PRE( elements->owner()==0 ) ;
   PEL_CHECK_PRE( elements->index_limit() >= 1 ) ;
   PEL_CHECK_PRE( elements->count() == elements->index_limit() ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t i=0 ; i<elements->index_limit() ; ++i ),
         dynamic_cast< PDE_ReferenceElement const* >( elements->at( i ) ) 
         != 0  ) ) ;
   PEL_CHECK_PRE( !has( field->name() ) ) ;
   PEL_SAVEOLD( size_t, nb_fields, nb_fields() ) ;

   FIELDS->append( field ) ;
   field->set_owner( this ) ;
   
   ON_BOUNDS.append( located_on_bounds ) ;
   
   ELMS->append( elements ) ;
   elements->set_owner( this ) ;
   
   PEL_CHECK_POST( field->owner() == this ) ;
   PEL_CHECK_POST( elements->owner() == this ) ;
   PEL_CHECK_POST( nb_fields() == OLD(nb_fields)+1 ) ;
   PEL_CHECK_POST( has( field->name() ) ) ;
   PEL_CHECK_POST( item( field->name() ) == field ) ;
   PEL_CHECK_POST( item( index_of( field->name() ) ) == field ) ;
   PEL_CHECK_POST( reference_elements( index_of( field->name() ) ) == 
                   elements ) ;
   PEL_CHECK_POST( boundary_located( index_of( field->name() ) ) == 
                   located_on_bounds ) ;
}

//----------------------------------------------------------------------
void
PDE_SetOfDiscreteFields:: start( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDiscreteFields:: start" ) ;
   INDEX = 0 ;
}

//----------------------------------------------------------------------
void
PDE_SetOfDiscreteFields:: go_next( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDiscreteFields:: go_next" ) ;
   INDEX++ ;
}

//----------------------------------------------------------------------
bool
PDE_SetOfDiscreteFields:: is_valid( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDiscreteFields:: is_valid" ) ;
   return( INDEX<nb_fields() ) ;
}

//----------------------------------------------------------------------
PDE_DiscreteField*
PDE_SetOfDiscreteFields:: item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDiscreteFields:: item()" ) ;
   PEL_CHECK_PRE( is_valid() ) ;

   PDE_DiscreteField* result = 
          static_cast<PDE_DiscreteField*>( FIELDS->at( INDEX ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_SetOfDiscreteFields:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string space( indent_width, ' ' ) ;
   os << space << nb_fields() << " discrete fields : "
      << std::endl ;
   for( start() ; is_valid() ; go_next() )
   {
      os << space << "   - \"" << item()->name() << "\"" << std::endl ;
   }
}

//----------------------------------------------------------------------
size_t
PDE_SetOfDiscreteFields:: find_index( std::string const& field_name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDiscreteFields:: find_index" ) ;
   
   size_t result = PEL::bad_index() ;
   for( size_t i=0 ; result==PEL::bad_index() && i<FIELDS->count() ; ++i )
   {
      PDE_DiscreteField* f =
                static_cast<PDE_DiscreteField*>( FIELDS->at(i) ) ;
      if( f->name()==field_name )
      {
         result = i ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_SetOfDiscreteFields:: save_state( PEL_ObjectWriter* writer ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDiscreteFields:: save_state" ) ;
   PEL_CHECK_PRE( save_state_PRE( writer ) ) ;

   writer->start_new_object( "PDE_SetOfDiscreteFields" ) ;
   for( start() ; is_valid() ; go_next() )
   {
      item()->save_state( writer ) ;
   }
   writer->finalize_object() ;

   PEL_CHECK_POST( save_state_POST( writer ) ) ;
}

//----------------------------------------------------------------------
void
PDE_SetOfDiscreteFields:: restore_state( PEL_ObjectReader* reader )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDiscreteFields:: restore_state" ) ;
   PEL_CHECK_PRE( restore_state_PRE( reader ) ) ;

   reader->start_object_retrieval( "PDE_SetOfDiscreteFields" ) ;
   for( start() ; is_valid() ; go_next() )
   {
      item()->restore_state( reader ) ;
   }
   reader->end_object_retrieval() ;

   PEL_CHECK_POST( restore_state_POST( reader ) ) ;
}
