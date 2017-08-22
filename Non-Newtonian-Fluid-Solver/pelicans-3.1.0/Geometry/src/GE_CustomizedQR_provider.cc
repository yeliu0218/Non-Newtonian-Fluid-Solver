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

#include <GE_CustomizedQR_provider.hh>

#include <PEL_assertions.hh>
#include <PEL_Vector.hh>

#include <GE_Customized_QR.hh>
#include <GE_ReferenceCube.hh>
#include <GE_ReferencePoint.hh>
#include <GE_ReferenceSegment.hh>
#include <GE_ReferenceSquare.hh>
#include <GE_ReferenceTetrahedron.hh>
#include <GE_ReferenceTriangle.hh>

//-----------------------------------------------------------------------------
GE_CustomizedQR_provider const*
GE_CustomizedQR_provider:: create( PEL_Object* a_owner )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CustomizedQR_provider:: create" ) ;

   GE_CustomizedQR_provider* result = new GE_CustomizedQR_provider( a_owner ) ;
   result->build() ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_CustomizedQR_provider:: GE_CustomizedQR_provider( PEL_Object* a_owner )
//-----------------------------------------------------------------------------
   : GE_QRprovider( a_owner, "GE_CustomizedQR_provider" )
   , CQR( PEL_Vector::create( this, 0 ) )
{
}

//-----------------------------------------------------------------------------
GE_CustomizedQR_provider:: ~GE_CustomizedQR_provider( void )
//-----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//-----------------------------------------------------------------------------
GE_Customized_QR*
GE_CustomizedQR_provider:: quadrature_rule_to_be_customized(
                                     GE_ReferencePolyhedron const* poly ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CustomizedQR_provider:: quadrature_rule_to_be_customized" ) ;
   PEL_CHECK_PRE( poly!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Customized_QR* result = qr( poly ) ;
   if( result->nb_points() != 0 )
   {
      notify_status_change() ;
      result->reset_points() ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->reference_polyhedron() == poly ) ;
   PEL_CHECK_POST( result->order() == 0 ) ;
   PEL_CHECK_POST( result->nb_points() == 0 ) ;
   PEL_CHECK_POST( !result->has_been_finalized() ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
GE_Customized_QR*
GE_CustomizedQR_provider:: qr( GE_ReferencePolyhedron const* poly ) const
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_CustomizedQR_provider:: qr" ) ;
   PEL_CHECK( poly!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_Customized_QR* result = 0 ;

   size_t idx = poly->id_number() ;
   if( CQR->index_limit()>idx )
   {
      result = static_cast<GE_Customized_QR*>( CQR->at(idx) ) ;
   }
   else
   {
      CQR->resize( idx+1 ) ;
   }
   if( result==0 )
   {
      result = GE_Customized_QR::create(
                    const_cast<GE_CustomizedQR_provider*>( this ), poly, 0 ) ;
      CQR->set_at( idx, result ) ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->reference_polyhedron() == poly ) ;
   PEL_CHECK_POST( result->order() == 0 ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------------
void
GE_CustomizedQR_provider:: build( void )
//-----------------------------------------------------------------------------
{
   add_QR( GE_ReferenceSquare::object(),
           qr( GE_ReferenceSquare::object() ) ) ;
   add_QR( GE_ReferenceTriangle::object(),
           qr( GE_ReferenceTriangle::object() ) ) ;
   add_QR( GE_ReferenceSegment::object(),
           qr( GE_ReferenceSegment::object() ) ) ;
   add_QR( GE_ReferenceCube::object(),
           qr( GE_ReferenceCube::object() ) ) ;
   add_QR( GE_ReferenceTetrahedron::object(),
           qr( GE_ReferenceTetrahedron::object() ) ) ;
   add_QR( GE_ReferencePoint::object(),
           qr( GE_ReferencePoint::object() ) ) ;
   PEL_CHECK_INV( invariant() ) ;
}
