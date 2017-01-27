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

#ifndef GE_CUSTOMIZED_QR_HH
#define GE_CUSTOMIZED_QR_HH

#include <GE_QuadratureRule.hh>

/*
Quadrature rules 
   - that initially contain no points ;
   - that have to be defined by the client using `::insert_point'.

PUBLISHED
*/

class PEL_EXPORT GE_Customized_QR : public GE_QuadratureRule
{
   public: //------------------------------------------------------------------

   //-- Instance delivery and initialization

      static GE_Customized_QR* create( PEL_Object* a_owner,
                                       GE_ReferencePolyhedron const* a_poly,
                                       size_t a_order ) ;

   //-- Rule building
      
      void reset_points( void ) ;
      
      void insert_point( GE_Point const* pt, double pt_weight ) ;

      void finalize( void ) ;

      bool has_been_finalized( void ) const ;

   protected: //---------------------------------------------------------------

   private: //-----------------------------------------------------------------

      GE_Customized_QR( void ) ;
     ~GE_Customized_QR( void ) ;
      GE_Customized_QR( GE_Customized_QR const& other ) ;
      GE_Customized_QR const& operator=( GE_Customized_QR const& other ) ;

      GE_Customized_QR( PEL_Object* a_owner,
                        GE_ReferencePolyhedron const* a_poly,
                        size_t a_order ) ;

   //-- Attributes

      bool FINALIZED ;
} ;

#endif
