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

#ifndef CH_P_SEMI_IMPLICIT_HH
#define CH_P_SEMI_IMPLICIT_HH

#include <CH_BulkEnergy.hh>

/*
PUBLISHED
*/

class CH_PSemiImplicit : public CH_BulkEnergy
{
   public: //-----------------------------------------------------------

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~CH_PSemiImplicit( void ) ;
      CH_PSemiImplicit( CH_PSemiImplicit const& other ) ;
      CH_PSemiImplicit& operator=( CH_PSemiImplicit const& other ) ;

      CH_PSemiImplicit( PEL_Object* a_owner,
                     double sigma_1,
                     double sigma_2,
                     double sigma_3,
                     PEL_ModuleExplorer const* exp ) ;

   //-- Plug in   
      
      CH_PSemiImplicit( void ) ;
      
      virtual CH_PSemiImplicit* create_replica( 
                                    PEL_Object* a_owner,
                                    double sigma_1,
                                    double sigma_2,
                                    double sigma_3,
                                    PEL_ModuleExplorer const* exp  ) const ;

   //-- Volumic free energy

      virtual double F( double c1, double c2, double c3 ) const ;

      virtual double diF( double c1, double c2, double c3,
                          double c1_exp, double c2_exp, double c3_exp,
                          size_t i ) const ; 

   //-- Bulk Chemical Potential
      
      virtual double dj_ddiF( double c1, double c2, 
                              double c1_exp, double c2_exp, 
                              size_t i, size_t j ) const ;

   //-- Class attributes
         
      static CH_PSemiImplicit const* PROTOTYPE ;
         
   //-- Attributes
      
      double TLAM ;
} ;

#endif
