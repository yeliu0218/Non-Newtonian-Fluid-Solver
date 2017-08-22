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

#include <PDE_BasisFunction.hh>

#include <PEL_Vector.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>

#include <PDE_BoundFE.hh>
#include <PDE_CellFE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_FaceFE.hh>
#include <PDE_ReferenceElement.hh>

#include <algorithm>
#include <iostream>
#include <list>
#include <string>

using std::cout ; using std::endl ;
using std::string ;

//-------------------------------------------------------------------------
PDE_BasisFunction:: PDE_BasisFunction( PEL_Object* a_owner )
//-------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , ACTIVE( false )
   , REFINED( false )
   , I_FIELD( PEL::bad_index() )
   , ID_NUMBER(  PEL::bad_index() )
{
}

//-------------------------------------------------------------------------
PDE_BasisFunction:: ~PDE_BasisFunction( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunction:: id_number( void ) const
//-------------------------------------------------------------------------
{
   return( ID_NUMBER ) ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunction:: set_id_number( size_t id )
//-------------------------------------------------------------------------
{
  ID_NUMBER = id ;
}

//-------------------------------------------------------------------------
bool
PDE_BasisFunction:: is_active( void ) const
//-------------------------------------------------------------------------
{
   return( ACTIVE ) ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunction:: set_active( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunction:: set_active" ) ;

   ACTIVE = true ;

   PEL_CHECK_POST( is_active() ) ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunction:: set_inactive( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunction:: set_inactive" ) ;

   ACTIVE = false ;

   PEL_CHECK_POST( !is_active() ) ;
}

//-------------------------------------------------------------------------
bool
PDE_BasisFunction:: is_refined( void ) const
//-------------------------------------------------------------------------
{
   return( REFINED ) ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunction:: set_refined( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunction:: set_refined" ) ;

   REFINED = true ;

   PEL_CHECK_POST( is_refined() ) ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunction:: set_unrefined( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunction:: set_unrefined" ) ;

   REFINED = false ;

   PEL_CHECK_POST( !is_refined() ) ;
}

//-------------------------------------------------------------------------
bool
PDE_BasisFunction:: is_attached_to_valid_DOFs( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunction:: is_attached_to_valid_DOFs" ) ;

   size_t ifield = 0 ;
   bool result = false ;
   for( ; !result && ifield<FIELDS.size() ; ++ifield )
   {
      if( FIELDS[ifield]->node_is_active( FIELD_NODES[ifield] ) )
      {
         result=true ;
      }
   }

   PEL_CHECK_POST( IMPLIES( is_active(), result ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunction:: append_one_field( PDE_DiscreteField* ff, size_t n )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunction:: append_one_field" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;
   PEL_CHECK_PRE( !is_in_basis_function_set_of( ff ) ) ;
   PEL_CHECK_PRE( n != PEL::bad_index() ) ;
   
   size_t const id = ff->id_number() ;
   if( FIELD_IDXS.size()<= id )
   {
      FIELD_IDXS.resize( id+1, PEL::bad_index() ) ;
   }
   FIELD_IDXS[id] = FIELDS.size() ;
   FIELDS.push_back( ff ) ;
   FIELD_NODES.push_back( n ) ;

   PEL_CHECK_POST( is_in_basis_function_set_of( ff ) ) ;
   PEL_CHECK_POST( node_of_DOF( ff ) == n ) ;
}

//-------------------------------------------------------------------------
bool
PDE_BasisFunction:: is_in_basis_function_set_of( 
                                         PDE_DiscreteField const* ff ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunction:: is_in_basis_function_set_of" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;

   size_t const id = ff->id_number() ;
   return( id<FIELD_IDXS.size() && FIELD_IDXS[id] != PEL::bad_index() ) ;
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunction:: node_of_DOF( PDE_DiscreteField const* ff ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunction:: node_of_DOF" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;
   PEL_CHECK_PRE( is_in_basis_function_set_of( ff ) ) ;
   return( FIELD_NODES[ FIELD_IDXS[ff->id_number()] ] ) ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunction:: start_field_iterator( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunction:: start_field_iterator" ) ;

   I_FIELD = 0 ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunction:: go_next_field( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunction:: go_next_field" ) ;
   PEL_CHECK_PRE( valid_field() ) ;

   ++I_FIELD ;
}

//-------------------------------------------------------------------------
bool
PDE_BasisFunction:: valid_field( void ) const
//-------------------------------------------------------------------------
{
   return( I_FIELD<FIELDS.size() ) ;
}

//-------------------------------------------------------------------------
PDE_DiscreteField*
PDE_BasisFunction:: field( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunction:: field" ) ;
   PEL_CHECK_PRE( valid_field() ) ;

   PDE_DiscreteField* result = FIELDS[I_FIELD] ;

   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
size_t
PDE_BasisFunction:: node_of_field( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunction:: node_of_field" ) ;
   PEL_CHECK_PRE( valid_field() ) ;

   size_t result = FIELD_NODES[I_FIELD] ;

   PEL_CHECK_POST( result != PEL::bad_index() ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
PDE_BasisFunction:: print( std::ostream& os, size_t indent_width ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BasisFunction:: print" ) ;

   string bl( indent_width, ' ' ) ;
   for( size_t i=0 ; i<FIELD_NODES.size() ; ++i )
   {
      os << bl << "\"" << FIELDS[i]->name()
         << "\" : node " << FIELD_NODES[ i ] << endl ;
   }
   os << bl << "id_number : " << ID_NUMBER << endl ;
}

//--------------------------------------------------------------------------
bool
PDE_BasisFunction:: is_located_in_cell_PRE( PDE_CellFE const* a_cell ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( a_cell != 0 ) ;
   PEL_ASSERT( refinement_level() <= a_cell->refinement_level() ) ;
   return( true ) ;
}

//--------------------------------------------------------------------------
bool
PDE_BasisFunction:: is_located_on_bound_PRE( PDE_BoundFE const* a_bound ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( a_bound != 0 ) ;
   PEL_ASSERT( refinement_level() <= a_bound->refinement_level() ) ;
   return( true ) ;
}
