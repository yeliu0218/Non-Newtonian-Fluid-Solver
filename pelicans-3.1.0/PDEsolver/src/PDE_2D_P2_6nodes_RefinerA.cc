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

#include <GE_ReferenceTriangleWithTriangles.hh>

#include <PDE_2D_P2_6nodes_RefinerA.hh>

//---------------------------------------------------------------------------
PDE_2D_P2_6nodes_RefinerA const*
PDE_2D_P2_6nodes_RefinerA:: object( void )
//---------------------------------------------------------------------------
{
   static PDE_2D_P2_6nodes_RefinerA* result = 0 ;
   if( result == 0 ) 
   {
      result = new PDE_2D_P2_6nodes_RefinerA() ;
   }
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_2D_P2_6nodes_RefinerA:: PDE_2D_P2_6nodes_RefinerA( void )
//---------------------------------------------------------------------------
   : PDE_ReferenceElementRefiner( 
                         "PDE_2D_P2_6nodes",
                         GE_ReferenceTriangleWithTriangles::object_2() )
{
   // *** coarse node 0 ***

   size_t a_parent_node = 0 ; size_t a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 0, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 1, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 5, 3./8., false ) ;

   a_parent_node = 0 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 1, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 5, -1./8., false ) ;

   a_parent_node = 0 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 1, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 5, -1./8., false ) ;

   a_parent_node = 0 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 1, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 5, -1./8., false ) ;

   // *** coarse node 1 ***

   a_parent_node = 1 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 1, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 2, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 3, 1./2., false ) ;

   a_parent_node = 1 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 0, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 1, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 5, 1./2., false ) ;

   a_parent_node = 1 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 1, 1./4., false ) ;
   append_child( a_parent_node, a_subcell, 3, 1./2., false ) ;
   append_child( a_parent_node, a_subcell, 4, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 5, 1./2., false ) ;

   a_parent_node = 1 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 1, 1./4., false ) ;

   // *** coarse node 2 ***

   a_parent_node = 2 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 1, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 3, -1./8., false ) ;

   a_parent_node = 2 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 1, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 2, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 3, 3./8., false ) ;

   a_parent_node = 2 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 1, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 3, -1./8., false ) ;

   a_parent_node = 2 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 1, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 3, -1./8., false ) ;

   // *** coarse node 3 ***

   a_parent_node = 3 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 3, 1./4., false ) ;

   a_parent_node = 3 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 3, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 4, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 5, 1./2., false ) ;

   a_parent_node = 3 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 0, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 1, 1./2., false ) ;
   append_child( a_parent_node, a_subcell, 3, 1./4., false ) ;
   append_child( a_parent_node, a_subcell, 5, 1./2., false ) ;

   a_parent_node = 3 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 1, 1./2., false ) ;
   append_child( a_parent_node, a_subcell, 2, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 3, 3./4., false ) ;

   // *** coarse node 4 ***

   a_parent_node = 4 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 3, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 5, -1./8., false ) ;

   a_parent_node = 4 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 3, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 5, -1./8., false ) ;

   a_parent_node = 4 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 3, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 5, -1./8., false ) ;

   a_parent_node = 4 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 3, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 5, 3./8., false ) ;

   // *** coarse node 5 ***

   a_parent_node = 5 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 3, 1./2., false ) ;
   append_child( a_parent_node, a_subcell, 4, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 5, 3./4., false ) ;

   a_parent_node = 5 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 5, 1./4., false ) ;

   a_parent_node = 5 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 1, 1./2., false ) ;
   append_child( a_parent_node, a_subcell, 2, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 3, 1./2., false ) ;
   append_child( a_parent_node, a_subcell, 5, 1./4., false ) ;

   a_parent_node = 5 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 0, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 1, 1./2., false ) ;
   append_child( a_parent_node, a_subcell, 5, 3./4., false ) ;
}

//---------------------------------------------------------------------------
PDE_2D_P2_6nodes_RefinerA:: ~PDE_2D_P2_6nodes_RefinerA( void )
//---------------------------------------------------------------------------
{
}
