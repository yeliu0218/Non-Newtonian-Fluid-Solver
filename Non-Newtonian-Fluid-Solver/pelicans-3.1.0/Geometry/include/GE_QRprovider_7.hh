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

#ifndef GE_QR_PROVIDER_7_HH
#define GE_QR_PROVIDER_7_HH

#include <GE_QRprovider.hh>

/*
Providers of quadrature rules of 7-th order.

PUBLISHED
*/

class PEL_EXPORT GE_QRprovider_7 : public GE_QRprovider
{
   public: //------------------------------------------------------------------
      
   protected: //---------------------------------------------------------------

   private: //-----------------------------------------------------------------

     ~GE_QRprovider_7( void ) ;
      GE_QRprovider_7( GE_QRprovider_7 const& other ) ;
      GE_QRprovider_7& operator=( GE_QRprovider_7 const& other ) ;

   //-- Plug in

      GE_QRprovider_7( void ) ;

      virtual void build( void ) ;
      
   //-- Class attributes

      static GE_QRprovider_7 const* REGISTRATOR ;
} ;

#endif
