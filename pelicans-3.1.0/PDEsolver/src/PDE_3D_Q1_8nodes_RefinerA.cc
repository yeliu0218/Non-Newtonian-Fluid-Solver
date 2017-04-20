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

#include <GE_ReferenceCubeWithCubes.hh>

#include <PDE_3D_Q1_8nodes_RefinerA.hh>

//---------------------------------------------------------------------------
PDE_3D_Q1_8nodes_RefinerA const*
PDE_3D_Q1_8nodes_RefinerA:: object( void )
//---------------------------------------------------------------------------
{
   static PDE_3D_Q1_8nodes_RefinerA* result = 0 ;
   if( result == 0 ) 
   {
      result = new PDE_3D_Q1_8nodes_RefinerA() ;
   }
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_3D_Q1_8nodes_RefinerA:: PDE_3D_Q1_8nodes_RefinerA( void )
//---------------------------------------------------------------------------
   : PDE_ReferenceElementRefiner( "PDE_3D_Q1_8nodes",
                                  GE_ReferenceCubeWithCubes::object_2() )
{
   // *** coarse node 0 ***

   size_t a_parent_node = 0 ; size_t a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 0, 1.000, true ) ;
   append_child( a_parent_node, a_subcell, 1, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.125, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.250, false ) ;

   a_parent_node = 0 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 0, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.125, false ) ;

   a_parent_node = 0 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 0, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.125, false ) ;

   a_parent_node = 0 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 0, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.125, false ) ;

   a_parent_node = 0 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 0, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.125, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.250, false ) ;

   a_parent_node = 0 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 0, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.125, false ) ;

   a_parent_node = 0 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 0, 0.125, false ) ;

   a_parent_node = 0 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 0, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.125, false ) ;

   // *** coarse node 1 ***
   // deduced from coarse node 0 by changing :
   // subcell:  0 <-> 1     child node:  0 <-> 1
   //           3 <-> 2                  3 <-> 2
   //           7 <-> 6                  4 <-> 5
   //           4 <-> 5                  7 <-> 6

   a_parent_node = 1 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 1, 1.000, true ) ;
   append_child( a_parent_node, a_subcell, 0, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.125, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.250, false ) ;

   a_parent_node = 1 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 1, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.125, false ) ;

   a_parent_node = 1 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 1, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.125, false ) ;

   a_parent_node = 1 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 1, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.125, false ) ;

   a_parent_node = 1 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 1, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.125, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.250, false ) ;

   a_parent_node = 1 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 1, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.125, false ) ;

   a_parent_node = 1 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 1, 0.125, false ) ;

   a_parent_node = 1 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 1, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.125, false ) ;

   // *** coarse node 2 ***
   // deduced from coarse node 1 by changing :
   // subcell:  1 <-> 2     child node:  0 <-> 3
   //           0 <-> 3                  1 <-> 2
   //           5 <-> 6                  4 <-> 7
   //           7 <-> 4                  5 <-> 6

   a_parent_node = 2 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 2, 1.000, true ) ;
   append_child( a_parent_node, a_subcell, 3, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.125, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.250, false ) ;

   a_parent_node = 2 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 2, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.125, false ) ;

   a_parent_node = 2 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 2, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.125, false ) ;

   a_parent_node = 2 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 2, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.125, false ) ;

   a_parent_node = 2 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 2, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.125, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.250, false ) ;

   a_parent_node = 2 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 2, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.125, false ) ;

   a_parent_node = 2 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 2, 0.125, false ) ;

   a_parent_node = 2 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 2, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.125, false ) ;

   // *** coarse node 3 ***
   // deduced from coarse node 0 by changing :
   // subcell:  1 <-> 2     child node:  0 <-> 3
   //           0 <-> 3                  1 <-> 2
   //           5 <-> 6                  4 <-> 7
   //           7 <-> 4                  5 <-> 6

   a_parent_node = 3 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 3, 1.000, true ) ;
   append_child( a_parent_node, a_subcell, 2, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.125, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.250, false ) ;

   a_parent_node = 3 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 3, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.125, false ) ;

   a_parent_node = 3 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 3, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.125, false ) ;

   a_parent_node = 3 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 3, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.125, false ) ;

   a_parent_node = 3 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 3, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.125, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.250, false ) ;

   a_parent_node = 3 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 3, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.125, false ) ;

   a_parent_node = 3 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 3, 0.125, false ) ;

   a_parent_node = 3 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 3, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.125, false ) ;

   // *** coarse node 4 ***
   // deduced from coarse node 0 by changing :
   // subcell:  0 <-> 4     child node:  0 <-> 4
   //           1 <-> 5                  1 <-> 5
   //           2 <-> 6                  2 <-> 6
   //           3 <-> 7                  3 <-> 7

   a_parent_node = 4 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 4, 1.000, true ) ;
   append_child( a_parent_node, a_subcell, 5, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.125, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.250, false ) ;

   a_parent_node = 4 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 4, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.125, false ) ;

   a_parent_node = 4 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 4, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.125, false ) ;

   a_parent_node = 4 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 4, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.125, false ) ;

   a_parent_node = 4 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 4, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.125, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.250, false ) ;

   a_parent_node = 4 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 4, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.125, false ) ;

   a_parent_node = 4 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 4, 0.125, false ) ;

   a_parent_node = 4 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 4, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.125, false ) ;

   // *** coarse node 5 ***
   // deduced from coarse node 1 by changing :
   // subcell:  0 <-> 4     child node:  0 <-> 4
   //           1 <-> 5                  1 <-> 5
   //           2 <-> 6                  2 <-> 6
   //           3 <-> 7                  3 <-> 7

   a_parent_node = 5 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 5, 1.000, true ) ;
   append_child( a_parent_node, a_subcell, 4, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.125, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.250, false ) ;

   a_parent_node = 5 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 5, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.125, false ) ;

   a_parent_node = 5 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 5, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.125, false ) ;

   a_parent_node = 5 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 5, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.125, false ) ;

   a_parent_node = 5 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 5, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.125, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.250, false ) ;

   a_parent_node = 5 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 5, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.125, false ) ;

   a_parent_node = 5 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 5, 0.125, false ) ;

   a_parent_node = 5 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 5, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.125, false ) ;

   // *** coarse node 6 ***
   // deduced from coarse node 2 by changing :
   // subcell:  0 <-> 4     child node:  0 <-> 4
   //           1 <-> 5                  1 <-> 5
   //           2 <-> 6                  2 <-> 6
   //           3 <-> 7                  3 <-> 7

   a_parent_node = 6 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 6, 1.000, true ) ;
   append_child( a_parent_node, a_subcell, 7, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.125, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.250, false ) ;

   a_parent_node = 6 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 6, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.125, false ) ;

   a_parent_node = 6 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 6, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.125, false ) ;

   a_parent_node = 6 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 6, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.125, false ) ;

   a_parent_node = 6 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 6, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.125, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.250, false ) ;

   a_parent_node = 6 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 6, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.125, false ) ;

   a_parent_node = 6 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 6, 0.125, false ) ;

   a_parent_node = 6 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 6, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 7, 0.125, false ) ;

   // *** coarse node 7 ***
   // deduced from coarse node 3 by changing :
   // subcell:  0 <-> 4     child node:  0 <-> 4
   //           1 <-> 5                  1 <-> 5
   //           2 <-> 6                  2 <-> 6
   //           3 <-> 7                  3 <-> 7

   a_parent_node = 7 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 7, 1.000, true ) ;
   append_child( a_parent_node, a_subcell, 6, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.500, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 1, 0.125, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.250, false ) ;

   a_parent_node = 7 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 7, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 0, 0.125, false ) ;

   a_parent_node = 7 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 7, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.125, false ) ;

   a_parent_node = 7 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 7, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 3, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 2, 0.125, false ) ;

   a_parent_node = 7 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 7, 0.50, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 5, 0.125, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.250, false ) ;

   a_parent_node = 7 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 7, 0.25, false ) ;
   append_child( a_parent_node, a_subcell, 4, 0.125, false ) ;

   a_parent_node = 7 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 7, 0.125, false ) ;

   a_parent_node = 7 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 7, 0.250, false ) ;
   append_child( a_parent_node, a_subcell, 6, 0.125, false ) ;
}

//---------------------------------------------------------------------------
PDE_3D_Q1_8nodes_RefinerA:: ~PDE_3D_Q1_8nodes_RefinerA( void )
//---------------------------------------------------------------------------
{
}
