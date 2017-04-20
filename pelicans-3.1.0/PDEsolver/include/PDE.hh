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

#ifndef PDE_HH
#define PDE_HH

#include <PEL_Object.hh>

class LA_Matrix ;
class LA_SeqMatrix ;
class LA_Vector ;

class PDE_LinkDOF2Unknown ;
class PDE_LocalEquation ;
class PDE_SystemNumbering ;

/*
PUBLISHED
*/

class PEL_EXPORT PDE
{
   public: //-----------------------------------------------------------
      
      static void assemble_in_matrix_vector_0( 
                                     LA_Matrix* matrix,
                                     LA_Vector* vector,
                                     PDE_LocalEquation const* leq,
                                     PDE_SystemNumbering const* r_nmb,
                                     size_t i_link,
                                     PDE_SystemNumbering const* c_nmb,
                                     size_t j_link ) ;
      
      static void assemble_in_matrix_vector_0( 
                                     LA_Matrix* matrix,
                                     LA_Vector* vector,
                                     PDE_LocalEquation const* leq,
                                     PDE_SystemNumbering const* r_nmb,
                                     PDE_SystemNumbering const* c_nmb ) ;
      
      static void assemble_in_matrix_vector_0( 
                                     LA_Matrix* matrix,
                                     LA_Vector* vector,
                                     PDE_LocalEquation const* leq,
                                     PDE_SystemNumbering const* nmb,
                                     size_t i_link,
                                     size_t j_link ) ;
      
      static void assemble_in_matrix_vector_0( 
                                     LA_Matrix* matrix,
                                     LA_Vector* vector,
                                     PDE_LocalEquation const* leq,
                                     PDE_SystemNumbering const* nmb ) ;
      
      static void assemble_in_matrix_0( LA_Matrix* matrix,
                                        PDE_LocalEquation const* leq ,
                                        PDE_SystemNumbering const* r_nmb,
                                        size_t i_link ,
                                        PDE_SystemNumbering const* c_nmb,
                                        size_t j_link ) ;

      static void assemble_in_matrix_0( LA_Matrix* matrix,
                                        PDE_LocalEquation const* leq ,
                                        PDE_SystemNumbering const* r_nmb,
                                        PDE_SystemNumbering const* c_nmb ) ;

      static void assemble_in_matrix_0( LA_Matrix* matrix,
                                        PDE_LocalEquation const* leq ,
                                        PDE_SystemNumbering const* nmb,
                                        size_t i_link ,
                                        size_t j_link ) ;

      static void assemble_in_matrix_0( LA_Matrix* matrix,
                                        PDE_LocalEquation const* leq ,
                                        PDE_SystemNumbering const* nmb ) ;

      static void assemble_in_vector_1( LA_Vector* vector,
                                        PDE_LocalEquation const* leq ,
                                        PDE_SystemNumbering const* nmb ) ;
      
      static void assemble_in_vector_1( LA_Vector* vector,
                                        PDE_LocalEquation const* leq ,
                                        PDE_SystemNumbering const* nmb,
                                        size_t i_link ) ;
      
      static LA_SeqMatrix* create_extracted_block( PEL_Object* a_owner,
                                            LA_Matrix const* matrix,
                                            PDE_SystemNumbering const* nmb,
                                            size_t i_link,
                                            size_t j_link ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------
      
      static void add_line( size_t il, size_t i, size_t ig,
                            LA_Matrix* matrix,
                            LA_Vector* vector,
                            PDE_LocalEquation const* leq,
                            double coef,
                            PDE_SystemNumbering const* c_nmb,
                            size_t j_link ) ;

} ; 

#endif
