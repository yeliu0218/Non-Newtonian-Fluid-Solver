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

#include <PDE_2D_P0_1node_RefinerA.hh>

//---------------------------------------------------------------------------
PDE_2D_P0_1node_RefinerA const*
PDE_2D_P0_1node_RefinerA:: object( void )
//---------------------------------------------------------------------------
{
   static PDE_2D_P0_1node_RefinerA* result = 0 ;
   if( result == 0 ) 
   {
      result = new PDE_2D_P0_1node_RefinerA() ;
   }
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_2D_P0_1node_RefinerA:: PDE_2D_P0_1node_RefinerA( void )
//---------------------------------------------------------------------------
   : PDE_ReferenceElementRefiner( "PDE_2D_P0_1node",
                                  GE_ReferenceTriangleWithTriangles::object_2() )
{
   size_t a_parent_node = 0 ; size_t a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 0, 1.0, false ) ;

   a_parent_node = 0 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 0, 1.0, false ) ;

   a_parent_node = 0 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 0, 1.0, false ) ;

   a_parent_node = 0 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 0, 1.0, false ) ;
}

//---------------------------------------------------------------------------
PDE_2D_P0_1node_RefinerA:: ~PDE_2D_P0_1node_RefinerA( void )
//---------------------------------------------------------------------------
{
}
