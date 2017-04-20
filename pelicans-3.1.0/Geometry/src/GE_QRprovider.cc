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

#include <GE_QRprovider.hh>

#include <PEL_Error.hh>
#include <PEL_Iterator.hh>
#include <PEL_ObjectRegister.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL_Vector.hh>

#include <GE_Customized_QR.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>

#include <iostream>

using std::string ;

//-------------------------------------------------------------------------
GE_QRprovider const*
GE_QRprovider:: object( std::string a_name )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "GE_QRprovider:: object" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;

   GE_QRprovider* result =
      static_cast<GE_QRprovider*>( plugins_map()->item( a_name ) ) ;

   // Build the inner quadrature rules at the first time :
   if( result->QRS->count() == 0 )
   {
      result->build() ;
   }
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( PEL_Root::object() ) ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_QRprovider:: GE_QRprovider( std::string const& a_name )
//-----------------------------------------------------------------------------
   : PEL_Object( plugins_map() )
   , NAME( a_name )
   , QRS( PEL_Vector::create( this, 0 ) )
   , STAT( 0 )
{
   PEL_LABEL( "GE_QRprovider:: GE_QRprovider" ) ;

   plugins_map()->register_item( a_name, this ) ;
   
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_QRprovider:: GE_QRprovider( PEL_Object* a_owner, std::string const& a_name )
//-----------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , NAME( a_name )
   , QRS( PEL_Vector::create( this, 0 ) )
   , STAT( 0 )
{
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_QRprovider:: ~GE_QRprovider( void )
//-----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
size_t
GE_QRprovider:: status( void ) const
//-----------------------------------------------------------------------------
{
   return( STAT ) ;
}

//-----------------------------------------------------------------------------
GE_QuadratureRule const*
GE_QRprovider:: quadrature_rule( GE_ReferencePolyhedron const* poly ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_QRprovider:: quadrature_rule" ) ;
   PEL_CHECK_PRE( poly != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   GE_QuadratureRule const* result = 0 ;
   size_t idx = poly->id_number() ;
   if( QRS->index_limit()>idx )
   {
      result = static_cast<GE_QuadratureRule const*>( QRS->at(idx) ) ;
   }

   if( result == 0 )
   {
      std::string msg = "*** GE_QRprovider : unknown quadrature rule defined\n" ;
      msg += "    provider name : \""+NAME+"\"\n" ;
      msg += "    reference polyhedron name : \""+poly->name()+"\"" ;
      PEL_Error::object()->raise_plain( msg ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( PEL_Root::object() ) ||
                   result->is_under_ownership_of( this ) ) ;
   PEL_CHECK_POST( result->reference_polyhedron() == poly ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Customized_QR const*
GE_QRprovider:: create_subset( PEL_Object* a_owner,
                               GE_Mpolyhedron const* pp,
                               PEL_Sequence const* poly_list ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_QRprovider:: create_subset" ) ;
   PEL_CHECK_PRE( create_subset_PRE( a_owner, pp, poly_list ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Iterator* it_poly = poly_list->create_iterator( 0 ) ;

   // Order :
   size_t order = 0 ;
   {
      it_poly->start() ;
      GE_Mpolyhedron const* m =
                     static_cast<GE_Mpolyhedron const*>( it_poly->item() ) ;
      GE_ReferencePolyhedron const* m_ref = m->reference_polyhedron() ;
      GE_QuadratureRule const* m_qr = quadrature_rule( m_ref ) ;
      order = m_qr->order() ;
   }

   // Building the integration points :
   GE_Customized_QR* result = 
                      GE_Customized_QR::create( a_owner,
                                                pp->reference_polyhedron(),
                                                order ) ;
   {
      GE_Point* pt_real = GE_Point::create( 0, pp->nb_space_dimensions() ) ;
      GE_Point* pp_pt_ref  = GE_Point::create( 0, pp->dimension() ) ;
      GE_ReferencePolyhedron const* pp_ref = pp->reference_polyhedron() ;
      for( it_poly->start() ; it_poly->is_valid() ; it_poly->go_next() )
      {
         PEL_CHECK( dynamic_cast<GE_Mpolyhedron const*>(it_poly->item())!=0 ) ;
         GE_Mpolyhedron const* m =
                     static_cast<GE_Mpolyhedron const*>( it_poly->item() ) ;
         GE_ReferencePolyhedron const* m_ref = m->reference_polyhedron() ;
         double const jac =
         ( m->measure()/m_ref->measure() )*( pp_ref->measure()/pp->measure() ) ;
         GE_QuadratureRule const* m_qr = quadrature_rule( m_ref ) ;
         for( size_t i=0 ; i<m_qr->nb_points() ; ++i )
         {
            GE_Point const* m_pt_ref = m_qr->point( i ) ;
            m->apply_mapping( m_pt_ref, pt_real ) ;
            pp->apply_inverse_mapping( pt_real, pp_pt_ref ) ;
            result->insert_point( pp_pt_ref, jac*m_qr->weight(i) ) ;
         }
      }
      result->finalize() ;
      pp_pt_ref->destroy() ; pp_pt_ref = 0 ;
      pt_real->destroy() ; pt_real = 0 ;
   }
   it_poly->destroy() ; it_poly = 0 ;
 
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( create_subset_POST( result, a_owner, pp, poly_list ) ) ;
   return( result ) ;

}

//-----------------------------------------------------------------------------
void
GE_QRprovider:: add_QR( GE_ReferencePolyhedron const* poly,
                        GE_QuadratureRule const* qr )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_QRprovider:: add_QR" ) ;
   PEL_CHECK_PRE( poly!=0 ) ;
   PEL_CHECK_PRE( qr!=0 ) ;
   PEL_CHECK_PRE( qr->reference_polyhedron() == poly ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_CHECK( QRS!=0 && !QRS->has(qr) ) ;

   size_t idx = poly->id_number() ;
   if( QRS->index_limit()<=idx )
   {
      QRS->resize( idx+1 ) ;
   }
   QRS->set_at( idx, const_cast<GE_QuadratureRule*>(qr) ) ;

   PEL_CHECK_POST( quadrature_rule( poly ) == qr ) ;
}

//-----------------------------------------------------------------------------
void
GE_QRprovider:: notify_status_change( void ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_QRprovider:: notify_status_change" ) ;
   PEL_SAVEOLD( size_t, status, status() ) ;

   ++STAT ;

   PEL_CHECK_POST( status() != OLD( status ) ) ;
}

//----------------------------------------------------------------------------
bool
GE_QRprovider:: create_subset_PRE( PEL_Object const* a_owner,
                                   GE_Mpolyhedron const* pp,
                                   PEL_Sequence const* poly_list ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( pp!=0 ) ;
   PEL_ASSERT( pp->dimension()==pp->nb_space_dimensions() ||
               pp->dimension()==pp->nb_space_dimensions()-1 ) ;
   PEL_ASSERT( poly_list!=0 ) ;
   PEL_ASSERT( poly_list->index_limit()==poly_list->count() ) ;
   PEL_ASSERT(
      FORALL(
         ( size_t i=0 ; i<poly_list->count() ; ++i ),
         dynamic_cast<GE_Mpolyhedron const*>( poly_list->at(i) )!=0 &&
         static_cast<GE_Mpolyhedron const*>( poly_list->at(i) )->nb_space_dimensions()==pp->nb_space_dimensions() &&
         static_cast<GE_Mpolyhedron const*>( poly_list->at(i) )->dimension()==pp->dimension() ) ) ;
   PEL_ASSERT(
      FORALL(
         ( size_t i=1 ; i<poly_list->count() ; ++i ),
         quadrature_rule( static_cast<GE_Mpolyhedron const*>( poly_list->at(i) )->reference_polyhedron() )->order()
         ==quadrature_rule( static_cast<GE_Mpolyhedron const*>( poly_list->at(0) )->reference_polyhedron() )->order() ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------------
bool
GE_QRprovider:: create_subset_POST( GE_Customized_QR const* result,
                                    PEL_Object const* a_owner,
                                    GE_Mpolyhedron const* pp,
                                    PEL_Sequence const* poly_list ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( result->has_been_finalized() ) ;
   PEL_ASSERT( result->reference_polyhedron()==pp->reference_polyhedron() ) ;
   PEL_ASSERT( result->order()==quadrature_rule( static_cast<GE_Mpolyhedron const*>( poly_list->at(0) )->reference_polyhedron() )->order() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------------
bool
GE_QRprovider:: invariant( void ) const
//----------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;

   return true ;
}

//----------------------------------------------------------------------
PEL_ObjectRegister*
GE_QRprovider:: plugins_map( void )
//----------------------------------------------------------------------
{
   static PEL_ObjectRegister* result =
              PEL_ObjectRegister::create( PEL_Root::object(),
                                          "GE_QRprovider descendant" ) ;
   return( result ) ;
}
