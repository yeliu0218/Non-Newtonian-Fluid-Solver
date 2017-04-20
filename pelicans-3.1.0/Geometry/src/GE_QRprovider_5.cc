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

#include <GE_QRprovider_5.hh>

#include <PEL_assertions.hh>

#include <GE_QuadratureRule.hh>
#include <GE_ReferenceCube.hh>
#include <GE_ReferencePoint.hh>
#include <GE_ReferenceSegment.hh>
#include <GE_ReferenceSquare.hh>
#include <GE_ReferenceTetrahedron.hh>
#include <GE_ReferenceTriangle.hh>

GE_QRprovider_5 const* 
GE_QRprovider_5:: REGISTRATOR = new GE_QRprovider_5() ;

//-----------------------------------------------------------------------------
GE_QRprovider_5:: GE_QRprovider_5( void )
//-----------------------------------------------------------------------------
   : GE_QRprovider( "GE_QRprovider_5" )
{
}

//-----------------------------------------------------------------------------
GE_QRprovider_5:: ~GE_QRprovider_5( void )
//-----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
   REGISTRATOR = 0 ;
}

//-----------------------------------------------------------------------------
void
GE_QRprovider_5:: build( void )
//-----------------------------------------------------------------------------
{
   // Point
   // -----
   add_QR( GE_ReferencePoint::object(),
           GE_QuadratureRule::object( "GE_Mpoint_QR1" ) ) ;
   
   // Segment
   // -------
   add_QR( GE_ReferenceSegment::object(),
           GE_QuadratureRule::object( "GE_Segment_QR5" ) ) ;
   
   // Square
   // ------
   add_QR( GE_ReferenceSquare::object(),
           GE_QuadratureRule::object( "GE_Square_QR5" ) ) ;
   
   // Triangle
   // --------
   add_QR( GE_ReferenceTriangle::object(),
           GE_QuadratureRule::object( "GE_Triangle_QR5" ) ) ;
   
   // Cube
   // ----
   add_QR( GE_ReferenceCube::object(),
           GE_QuadratureRule::object( "GE_Cube_QR5" ) ) ;
   
   // Tetrahedron
   // -----------
   add_QR( GE_ReferenceTetrahedron::object(),
           GE_QuadratureRule::object( "GE_Tetrahedron_QR5" ) ) ;

   PEL_CHECK_INV( invariant() ) ;
}



