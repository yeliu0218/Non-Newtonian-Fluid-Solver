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

#include <PDE_3D_Q2_27nodes_RefinerA.hh>

//---------------------------------------------------------------------------
PDE_3D_Q2_27nodes_RefinerA const*
PDE_3D_Q2_27nodes_RefinerA:: object( void )
//---------------------------------------------------------------------------
{
   static PDE_3D_Q2_27nodes_RefinerA* result = 0 ;
   if( result == 0 ) 
   {
      result = new PDE_3D_Q2_27nodes_RefinerA() ;
   }
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_3D_Q2_27nodes_RefinerA:: PDE_3D_Q2_27nodes_RefinerA( void )
//---------------------------------------------------------------------------
   : PDE_ReferenceElementRefiner( "PDE_3D_Q2_27nodes",
                                  GE_ReferenceCubeWithCubes::object_2() )
{
   // *** coarse node 0 ***

   size_t a_parent_node = 0 ; size_t a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 0, 1.000, true ) ;
   append_child( a_parent_node, a_subcell, 1, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 5, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 9, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./512., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./64., false ) ;

   a_parent_node = 0 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 1, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./64, false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;

   a_parent_node = 0 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 4, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;

   a_parent_node = 0 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 4, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 5, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./64., false ) ;

   a_parent_node = 0 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 9, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./64., false ) ;

   a_parent_node = 0 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 10, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;

   a_parent_node = 0 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 13, -1./512., false ) ;

   a_parent_node = 0 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;
   append_child( a_parent_node, a_subcell, 14, 1./64., false ) ;

   // *** coarse node 1 ***

   a_parent_node = 1 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 1, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 2, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 3, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 11, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;

   a_parent_node = 1 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 0, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 1, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 5, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 9, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./64., false ) ;

   a_parent_node = 1 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 4, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 5, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./64., false ) ;

   a_parent_node = 1 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 3, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;

   a_parent_node = 1 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 10, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 11, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256, false ) ;

   a_parent_node = 1 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 9, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./64., false ) ;

   a_parent_node = 1 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;
   append_child( a_parent_node, a_subcell, 14, 1./64., false ) ;

   a_parent_node = 1 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 12, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;

   // *** coarse node 2 ***
   // deduced from coarse node 0 by changing :
   // subcell:  0 <-> 1     child node:  8 <-> 6      15 <-> 17
   //           3 <-> 2                  3 <-> 5      14 <-> 12
   //           7 <-> 6                  2 <-> 0       9 <-> 11
   //           4 <-> 5                 

   a_parent_node = 2 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 2, 1.000, true ) ;
   append_child( a_parent_node, a_subcell, 1, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 3, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 11, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./512., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./64., false ) ;

   a_parent_node = 2 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 1, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./64, false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;

   a_parent_node = 2 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 4, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;

   a_parent_node = 2 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 4, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 3, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./64., false ) ;

   a_parent_node = 2 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 11, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./64., false ) ;

   a_parent_node = 2 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 10, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;

   a_parent_node = 2 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 13, -1./512., false ) ;

   a_parent_node = 2 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;
   append_child( a_parent_node, a_subcell, 12, 1./64., false ) ;

   // *** coarse node 3 ***
   // deduced from coarse node 5 by changing :
   // subcell:  0 <-> 1     child node:  8 <-> 6      15 <-> 17
   //           3 <-> 2                  3 <-> 5      14 <-> 12
   //           7 <-> 6                  2 <-> 0       9 <-> 11
   //           4 <-> 5                 

   a_parent_node = 3 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 3, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 8, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 7, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 17, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;

   a_parent_node = 3 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 2, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 3, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 1, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 11, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./64., false ) ;

   a_parent_node = 3 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 4, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 1, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./64., false ) ;

   a_parent_node = 3 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 7, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;

   a_parent_node = 3 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 12, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 17, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256, false ) ;

   a_parent_node = 3 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 11, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./64., false ) ;

   a_parent_node = 3 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;
   append_child( a_parent_node, a_subcell, 10, 1./64., false ) ;

   a_parent_node = 3 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 16, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;

   // *** coarse node 4 ***

   a_parent_node = 4 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 3, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 7, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 8, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 12, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 17, 3./8., false ) ;

   a_parent_node = 4 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 4, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 5, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 6, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 7, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 15, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./32., false ) ;

   a_parent_node = 4 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 0, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 1, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 5, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 9, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./32., false ) ;

   a_parent_node = 4 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 1, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 2, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 3, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 11, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;

   a_parent_node = 4 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 12, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128, false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 17, -1./8., false ) ;

   a_parent_node = 4 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 15, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./32., false ) ;

   a_parent_node = 4 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 9, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;

   a_parent_node = 4 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 10, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 11, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;

   // *** coarse node 5 ***
   // deduced from coarse node 1 by changing :
   // subcell:  1 <-> 3     child node:  2 <-> 6    11 <-> 15    20 <-> 24
   //           5 <-> 7                  3 <-> 7    12 <-> 16    21 <-> 25
   //                                    1 <-> 5    10 <-> 14    19 <-> 23

   a_parent_node = 5 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 5, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 6, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 7, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 15, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;

   a_parent_node = 5 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 0, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 5, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 1, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 9, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./64., false ) ;

   a_parent_node = 5 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 4, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 1, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./64., false ) ;

   a_parent_node = 5 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 7, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;

   a_parent_node = 5 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 14, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 15, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256, false ) ;

   a_parent_node = 5 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 9, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./64., false ) ;

   a_parent_node = 5 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;
   append_child( a_parent_node, a_subcell, 10, 1./64., false ) ;

   a_parent_node = 5 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 16, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;

   // *** coarse node 6 ***
   // deduced from coarse node 0 by changing :
   // subcell:  1 <-> 2     child node:  2 <-> 8      11 <-> 17
   //           0 <-> 3                  1 <-> 7      10 <-> 16
   //           5 <-> 6                  0 <-> 6       9 <-> 15
   //           7 <-> 4                 

   a_parent_node = 6 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 6, 1.000, true ) ;
   append_child( a_parent_node, a_subcell, 7, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 5, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 15, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./512., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./64., false ) ;

   a_parent_node = 6 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 7, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./64, false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;

   a_parent_node = 6 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 4, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;

   a_parent_node = 6 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 4, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 5, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./64., false ) ;

   a_parent_node = 6 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 15, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./64., false ) ;

   a_parent_node = 6 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 16, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;

   a_parent_node = 6 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 13, -1./512., false ) ;

   a_parent_node = 6 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;
   append_child( a_parent_node, a_subcell, 14, 1./64., false ) ;

   // *** coarse node 7 ***
   // deduced from coarse node 1 by changing :
   // subcell:  1 <-> 2     child node:  2 <-> 8      11 <-> 17
   //           0 <-> 3                  1 <-> 7      10 <-> 16
   //           5 <-> 6                  0 <-> 6       9 <-> 15
   //           7 <-> 4                 

   a_parent_node = 7 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 7, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 8, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 3, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 17, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;

   a_parent_node = 7 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 6, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 7, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 5, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 15, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./64., false ) ;

   a_parent_node = 7 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 4, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 5, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./64., false ) ;

   a_parent_node = 7 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 3, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;

   a_parent_node = 7 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 16, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 17, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256, false ) ;

   a_parent_node = 7 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 15, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./64., false ) ;

   a_parent_node = 7 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;
   append_child( a_parent_node, a_subcell, 14, 1./64., false ) ;

   a_parent_node = 7 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 12, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;

   // *** coarse node 8 ***
   // deduced from coarse node 2 by changing :
   // subcell:  1 <-> 2     child node:  2 <-> 8      11 <-> 17
   //           0 <-> 3                  1 <-> 7      10 <-> 16
   //           5 <-> 6                  0 <-> 6       9 <-> 15
   //           7 <-> 4                 

   a_parent_node = 8 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 8, 1.000, true ) ;
   append_child( a_parent_node, a_subcell, 7, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 3, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 17, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./512., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./64., false ) ;

   a_parent_node = 8 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 7, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./64, false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;

   a_parent_node = 8 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 4, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;

   a_parent_node = 8 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 4, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 3, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./64., false ) ;

   a_parent_node = 8 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 17, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./64., false ) ;

   a_parent_node = 8 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 16, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;

   a_parent_node = 8 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 13, -1./512., false ) ;

   a_parent_node = 8 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;
   append_child( a_parent_node, a_subcell, 12, 1./64., false ) ;

   // *** coarse node 9 ***
   // deduced from coarse node 1 by changing :
   // subcell:  1 <-> 4     child node:  1 <-> 9    4  <-> 14    7  <-> 15
   //           2 <-> 7                  2 <-> 18   3  <-> 23    8  <-> 24
   //                                    11<-> 19   12 <-> 22    17 <-> 25

   a_parent_node = 9 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 9, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 18, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 23, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 19, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;

   a_parent_node = 9 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 0, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 9, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 5, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 1, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./64., false ) ;

   a_parent_node = 9 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 14, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 5, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./64., false ) ;

   a_parent_node = 9 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 23, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;

   a_parent_node = 9 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 10, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 19, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256, false ) ;

   a_parent_node = 9 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 1, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./64., false ) ;

   a_parent_node = 9 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;
   append_child( a_parent_node, a_subcell, 4, 1./64., false ) ;

   a_parent_node = 9 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 22, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;

   // *** coarse node 10 ***
   // deduced from coarse node 4 by changing :
   // subcell:  3 <-> 4     child node:  8  <-> 20   7  <-> 19    6  <-> 18
   //           2 <-> 5                  3  <-> 11   4  <-> 10    5  <-> 9
   //                                    17 <-> 21   16 <-> 22    15 <-> 23

   a_parent_node = 10 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 11, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 19, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 20, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 12, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 21, 3./8., false ) ;

   a_parent_node = 10 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 10, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 9, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 18, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 19, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 23, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./32., false ) ;

   a_parent_node = 10 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 0, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 1, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 9, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 5, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./32., false ) ;

   a_parent_node = 10 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 1, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 2, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 11, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 3, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;

   a_parent_node = 10 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 12, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128, false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 21, -1./8., false ) ;

   a_parent_node = 10 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 23, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./32., false ) ;

   a_parent_node = 10 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 5, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;

   a_parent_node = 10 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 4, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 3, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;

   // *** coarse node 11 ***
   // deduced from coarse node 9 by changing :
   // subcell:  0 <-> 1     child node:  8 <-> 6      15 <-> 17     18 <-> 20
   //           3 <-> 2                  3 <-> 5      14 <-> 12     23 <-> 21
   //           7 <-> 6                  2 <-> 0       9 <-> 11     24 <-> 26
   //           4 <-> 5                 

   a_parent_node = 11 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 11, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 20, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 21, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 19, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;

   a_parent_node = 11 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 2, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 11, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 3, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 1, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./64., false ) ;

   a_parent_node = 11 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 12, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 3, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./64., false ) ;

   a_parent_node = 11 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 21, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;

   a_parent_node = 11 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 10, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 19, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256, false ) ;

   a_parent_node = 11 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 1, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./64., false ) ;

   a_parent_node = 11 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;
   append_child( a_parent_node, a_subcell, 4, 1./64., false ) ;

   a_parent_node = 11 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 22, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;

   // *** coarse node 12 ***
   // deduced from coarse node 14 by changing :
   // subcell:  0 <-> 1     child node:  8 <-> 6      15 <-> 17     18 <-> 20
   //           3 <-> 2                  3 <-> 5      14 <-> 12     23 <-> 21
   //           7 <-> 6                  2 <-> 0       9 <-> 11     24 <-> 26
   //           4 <-> 5                 

   a_parent_node = 12 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 21, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 17, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 26, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 22, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 25, 3./8., false ) ;

   a_parent_node = 12 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 12, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 3, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 8, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 17, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 7, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./32., false ) ;

   a_parent_node = 12 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 2, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 11, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 3, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 1, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./32., false ) ;

   a_parent_node = 12 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 11, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 20, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 21, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 19, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;

   a_parent_node = 12 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 22, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128, false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 25, -1./8., false ) ;

   a_parent_node = 12 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 7, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./32., false ) ;

   a_parent_node = 12 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 1, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;

   a_parent_node = 12 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 10, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 19, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;

   // *** coarse node 13 ***

   a_parent_node = 13 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 12, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./64., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 17, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 21, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 25, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 26, 1.0, true ) ;

   a_parent_node = 13 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 13, 27./64., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 15, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 23, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 24, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 25, 3./4., false ) ;

   a_parent_node = 13 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 9, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./64., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 18, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 19, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 23, 3./4., false ) ;

   a_parent_node = 13 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 10, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 11, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./64., false ) ;
   append_child( a_parent_node, a_subcell, 19, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 20, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 21, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./16., false ) ;

   a_parent_node = 13 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 3, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 7, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 8, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 12, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./64., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 17, 3./4., false ) ;

   a_parent_node = 13 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 4, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 5, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 6, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 7, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./64., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 15, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./16., false ) ;

   a_parent_node = 13 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 0, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 1, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 5, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 9, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./64., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./16., false ) ;

   a_parent_node = 13 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 1, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 2, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 3, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 11, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./64., false ) ;

   // *** coarse node 14 ***
   // deduced from coarse node 4 by changing :
   // subcell:  1 <-> 4     child node:  1 <-> 9    4  <-> 14    7  <-> 15
   //           2 <-> 7                  2 <-> 18   3  <-> 23    8  <-> 24
   //                                    11<-> 19   12 <-> 22    17 <-> 25

   a_parent_node = 14 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 23, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 15, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 24, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 22, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 25, 3./8., false ) ;

   a_parent_node = 14 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 14, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 5, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 6, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 15, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 7, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./32., false ) ;

   a_parent_node = 14 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 0, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 9, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 5, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 1, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./32., false ) ;

   a_parent_node = 14 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 9, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 18, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 23, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 19, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;

   a_parent_node = 14 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 22, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128, false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 25, -1./8., false ) ;

   a_parent_node = 14 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 7, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./32., false ) ;

   a_parent_node = 14 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 1, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;

   a_parent_node = 14 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 10, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 19, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;

   // *** coarse node 15 ***
   // deduced from coarse node 7 by changing :
   // subcell:  1 <-> 4     child node:  1 <-> 9    4  <-> 14    7  <-> 15
   //           2 <-> 7                  2 <-> 18   3  <-> 23    8  <-> 24
   //                                    11<-> 19   12 <-> 22    17 <-> 25

   a_parent_node = 15 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 15, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 24, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 23, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 25, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;

   a_parent_node = 15 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 6, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 15, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 5, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 7, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./64., false ) ;

   a_parent_node = 15 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 14, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 5, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./64., false ) ;

   a_parent_node = 15 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 23, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;

   a_parent_node = 15 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 16, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 25, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256, false ) ;

   a_parent_node = 15 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 7, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./64., false ) ;

   a_parent_node = 15 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;
   append_child( a_parent_node, a_subcell, 4, 1./64., false ) ;

   a_parent_node = 15 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 22, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;

   // *** coarse node 16 ***
   // deduced from coarse node 10 by changing :
   // subcell:  1 <-> 2     child node:  2 <-> 8      11 <-> 17     18 <-> 24
   //           0 <-> 3                  1 <-> 7      10 <-> 16     19 <-> 25
   //           5 <-> 6                  0 <-> 6       9 <-> 15     20 <-> 26
   //           7 <-> 4                 

   a_parent_node = 16 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 17, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 25, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 26, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 12, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 21, 3./8., false ) ;

   a_parent_node = 16 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 16, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 15, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 24, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 25, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 23, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./32., false ) ;

   a_parent_node = 16 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 6, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 7, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 15, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 5, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./32., false ) ;

   a_parent_node = 16 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 7, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 8, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 17, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 3, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;

   a_parent_node = 16 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 12, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128, false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 21, -1./8., false ) ;

   a_parent_node = 16 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 23, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./32., false ) ;

   a_parent_node = 16 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 5, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;

   a_parent_node = 16 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 4, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 3, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;

   // *** coarse node 17 ***
   // deduced from coarse node 15 by changing :
   // subcell:  0 <-> 1     child node:  8 <-> 6      15 <-> 17     18 <-> 20
   //           3 <-> 2                  3 <-> 5      14 <-> 12     23 <-> 21
   //           7 <-> 6                  2 <-> 0       9 <-> 11     24 <-> 26
   //           4 <-> 5          
       
   a_parent_node = 17 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 17, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 26, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 21, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 25, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;

   a_parent_node = 17 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 8, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 17, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 3, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 7, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;
   append_child( a_parent_node, a_subcell, 4, 9./64., false ) ;

   a_parent_node = 17 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 12, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 3, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./64., false ) ;

   a_parent_node = 17 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 21, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;

   a_parent_node = 17 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 16, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 25, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256, false ) ;

   a_parent_node = 17 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 7, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 4, -3./64., false ) ;

   a_parent_node = 17 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;
   append_child( a_parent_node, a_subcell, 4, 1./64., false ) ;

   a_parent_node = 17 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 22, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;

   // *** coarse node 18 ***
   // deduced from coarse node 0 by changing :
   // subcell:  0 <-> 4     child node:  0 <-> 18       5 <-> 23
   //           1 <-> 5                  1 <-> 19       6 <-> 24
   //           2 <-> 6                  2 <-> 20       7 <-> 25
   //           3 <-> 7                  3 <-> 21       8 <-> 26
   //                                    4 <-> 22

   a_parent_node = 18 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 18, 1.000, true ) ;
   append_child( a_parent_node, a_subcell, 19, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 23, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 9, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./512., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./64., false ) ;

   a_parent_node = 18 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 19, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./64, false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;

   a_parent_node = 18 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 22, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;

   a_parent_node = 18 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 22, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 23, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./64., false ) ;

   a_parent_node = 18 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 9, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./64., false ) ;

   a_parent_node = 18 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 10, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;

   a_parent_node = 18 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 13, -1./512., false ) ;

   a_parent_node = 18 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;
   append_child( a_parent_node, a_subcell, 14, 1./64., false ) ;

   // *** coarse node 19 ***
   // deduced from coarse node 1 by changing :
   // subcell:  0 <-> 4     child node:  0 <-> 18       5 <-> 23
   //           1 <-> 5                  1 <-> 19       6 <-> 24
   //           2 <-> 6                  2 <-> 20       7 <-> 25
   //           3 <-> 7                  3 <-> 21       8 <-> 26
   //                                    4 <-> 22

   a_parent_node = 19 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 19, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 20, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 21, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 11, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;

   a_parent_node = 19 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 18, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 19, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 23, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 9, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./64., false ) ;

   a_parent_node = 19 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 22, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 23, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./64., false ) ;

   a_parent_node = 19 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 21, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;

   a_parent_node = 19 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 10, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 11, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256, false ) ;

   a_parent_node = 19 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 9, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./64., false ) ;

   a_parent_node = 19 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;
   append_child( a_parent_node, a_subcell, 14, 1./64., false ) ;

   a_parent_node = 19 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 12, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;

   // *** coarse node 20 ***
   // deduced from coarse node 2 by changing :
   // subcell:  0 <-> 4     child node:  0 <-> 18       5 <-> 23
   //           1 <-> 5                  1 <-> 19       6 <-> 24
   //           2 <-> 6                  2 <-> 20       7 <-> 25
   //           3 <-> 7                  3 <-> 21       8 <-> 26
   //                                    4 <-> 22

   a_parent_node = 20 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 20, 1.000, true ) ;
   append_child( a_parent_node, a_subcell, 19, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 21, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 11, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./512., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./64., false ) ;

   a_parent_node = 20 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 19, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./64, false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;

   a_parent_node = 20 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 22, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;

   a_parent_node = 20 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 22, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 21, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./64., false ) ;

   a_parent_node = 20 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 11, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./64., false ) ;

   a_parent_node = 20 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 10, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;

   a_parent_node = 20 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 13, -1./512., false ) ;

   a_parent_node = 20 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;
   append_child( a_parent_node, a_subcell, 12, 1./64., false ) ;

   // *** coarse node 21 ***
   // deduced from coarse node 3 by changing :
   // subcell:  0 <-> 4     child node:  0 <-> 18       5 <-> 23
   //           1 <-> 5                  1 <-> 19       6 <-> 24
   //           2 <-> 6                  2 <-> 20       7 <-> 25
   //           3 <-> 7                  3 <-> 21       8 <-> 26
   //                                    4 <-> 22

   a_parent_node = 21 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 21, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 26, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 25, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 17, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;

   a_parent_node = 21 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 20, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 21, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 19, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 11, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./64., false ) ;

   a_parent_node = 21 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 22, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 19, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./64., false ) ;

   a_parent_node = 21 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 25, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;

   a_parent_node = 21 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 12, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 17, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256, false ) ;

   a_parent_node = 21 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 11, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./64., false ) ;

   a_parent_node = 21 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;
   append_child( a_parent_node, a_subcell, 10, 1./64., false ) ;

   a_parent_node = 21 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 16, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;

   // *** coarse node 22 ***
   // deduced from coarse node 4 by changing :
   // subcell:  0 <-> 4     child node:  0 <-> 18       5 <-> 23
   //           1 <-> 5                  1 <-> 19       6 <-> 24
   //           2 <-> 6                  2 <-> 20       7 <-> 25
   //           3 <-> 7                  3 <-> 21       8 <-> 26
   //                                    4 <-> 22

   a_parent_node = 22 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 21, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 25, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 26, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 12, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 17, 3./8., false ) ;

   a_parent_node = 22 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 22, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 23, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 24, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 25, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 15, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./32., false ) ;

   a_parent_node = 22 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 18, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 19, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 23, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 9, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./32., false ) ;

   a_parent_node = 22 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 19, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 20, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 21, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./16., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 11, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./128., false ) ;

   a_parent_node = 22 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 12, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128, false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 17, -1./8., false ) ;

   a_parent_node = 22 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 15, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./32., false ) ;

   a_parent_node = 22 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 9, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;

   a_parent_node = 22 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 10, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 11, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./128., false ) ;

   // *** coarse node 23 ***
   // deduced from coarse node 5 by changing :
   // subcell:  0 <-> 4     child node:  0 <-> 18       5 <-> 23
   //           1 <-> 5                  1 <-> 19       6 <-> 24
   //           2 <-> 6                  2 <-> 20       7 <-> 25
   //           3 <-> 7                  3 <-> 21       8 <-> 26
   //                                    4 <-> 22

   a_parent_node = 23 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 23, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 24, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 25, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 15, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;

   a_parent_node = 23 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 18, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 23, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 19, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 9, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;
   append_child( a_parent_node, a_subcell, 10, 9./64., false ) ;

   a_parent_node = 23 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 22, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 19, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./64., false ) ;

   a_parent_node = 23 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 25, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;

   a_parent_node = 23 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 14, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 15, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256, false ) ;

   a_parent_node = 23 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 9, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 10, -3./64., false ) ;

   a_parent_node = 23 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;
   append_child( a_parent_node, a_subcell, 10, 1./64., false ) ;

   a_parent_node = 23 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 16, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;

   // *** coarse node 24 ***
   // deduced from coarse node 6 by changing :
   // subcell:  0 <-> 4     child node:  0 <-> 18       5 <-> 23
   //           1 <-> 5                  1 <-> 19       6 <-> 24
   //           2 <-> 6                  2 <-> 20       7 <-> 25
   //           3 <-> 7                  3 <-> 21       8 <-> 26
   //                                    4 <-> 22

   a_parent_node = 24 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 24, 1.000, true ) ;
   append_child( a_parent_node, a_subcell, 25, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 23, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 15, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./512., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./64., false ) ;

   a_parent_node = 24 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 25, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./64, false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;

   a_parent_node = 24 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 22, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;

   a_parent_node = 24 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 22, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 23, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./64., false ) ;

   a_parent_node = 24 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 15, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./64., false ) ;

   a_parent_node = 24 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 16, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;

   a_parent_node = 24 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 13, -1./512., false ) ;

   a_parent_node = 24 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;
   append_child( a_parent_node, a_subcell, 14, 1./64., false ) ;

   // *** coarse node 25 ***
   // deduced from coarse node 7 by changing :
   // subcell:  0 <-> 4     child node:  0 <-> 18       5 <-> 23
   //           1 <-> 5                  1 <-> 19       6 <-> 24
   //           2 <-> 6                  2 <-> 20       7 <-> 25
   //           3 <-> 7                  3 <-> 21       8 <-> 26
   //                                    4 <-> 22

   a_parent_node = 25 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 25, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 26, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 21, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 17, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;

   a_parent_node = 25 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 24, 1.0, true ) ;
   append_child( a_parent_node, a_subcell, 25, 3./4., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 23, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 15, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./256., false ) ;
   append_child( a_parent_node, a_subcell, 14, 9./64., false ) ;

   a_parent_node = 25 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 22, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 23, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./64., false ) ;

   a_parent_node = 25 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 21, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;

   a_parent_node = 25 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 16, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 17, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256, false ) ;

   a_parent_node = 25 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 15, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./32., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./256., false ) ;
   append_child( a_parent_node, a_subcell, 14, -3./64., false ) ;

   a_parent_node = 25 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;
   append_child( a_parent_node, a_subcell, 14, 1./64., false ) ;

   a_parent_node = 25 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 12, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./256., false ) ;

   // *** coarse node 26 ***
   // deduced from coarse node 8 by changing :
   // subcell:  0 <-> 4     child node:  0 <-> 18       5 <-> 23
   //           1 <-> 5                  1 <-> 19       6 <-> 24
   //           2 <-> 6                  2 <-> 20       7 <-> 25
   //           3 <-> 7                  3 <-> 21       8 <-> 26
   //                                    4 <-> 22

   a_parent_node = 26 ; a_subcell = 7 ;
   append_child( a_parent_node, a_subcell, 26, 1.000, true ) ;
   append_child( a_parent_node, a_subcell, 25, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 21, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 17, 3./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, 9./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 27./512., false ) ;
   append_child( a_parent_node, a_subcell, 12, 9./64., false ) ;

   a_parent_node = 26 ; a_subcell = 5 ;
   append_child( a_parent_node, a_subcell, 25, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 22, -3./64, false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;

   a_parent_node = 26 ; a_subcell = 4 ;
   append_child( a_parent_node, a_subcell, 22, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;

   a_parent_node = 26 ; a_subcell = 6 ;
   append_child( a_parent_node, a_subcell, 22, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 21, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./64., false ) ;

   a_parent_node = 26 ; a_subcell = 3 ;
   append_child( a_parent_node, a_subcell, 17, -1./8., false ) ;
   append_child( a_parent_node, a_subcell, 16, -3./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, -9./512., false ) ;
   append_child( a_parent_node, a_subcell, 12, -3./64., false ) ;

   a_parent_node = 26 ; a_subcell = 1 ;
   append_child( a_parent_node, a_subcell, 16, 1./64., false ) ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;

   a_parent_node = 26 ; a_subcell = 0 ;
   append_child( a_parent_node, a_subcell, 13, -1./512., false ) ;

   a_parent_node = 26 ; a_subcell = 2 ;
   append_child( a_parent_node, a_subcell, 13, 3./512., false ) ;
   append_child( a_parent_node, a_subcell, 12, 1./64., false ) ;
}                             

//---------------------------------------------------------------------------
PDE_3D_Q2_27nodes_RefinerA:: ~PDE_3D_Q2_27nodes_RefinerA( void )
//---------------------------------------------------------------------------
{
}
