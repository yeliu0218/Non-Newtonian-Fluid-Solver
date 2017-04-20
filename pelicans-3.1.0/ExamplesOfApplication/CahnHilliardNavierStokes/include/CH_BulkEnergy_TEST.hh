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

#ifndef CH_BULK_ENERGY_TEST_HH
#define CH_BULK_ENERGY_TEST_HH

#include <PEL_ObjectTest.hh>

class CH_BulkEnergy ;

/*
PUBLISHED
*/

class CH_BulkEnergy_TEST : public PEL_ObjectTest
{
   public: //-----------------------------------------------------------------

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

      CH_BulkEnergy_TEST( void ) ;
     ~CH_BulkEnergy_TEST( void ) ;
      CH_BulkEnergy_TEST( CH_BulkEnergy_TEST const& other ) ;
      CH_BulkEnergy_TEST& operator=( CH_BulkEnergy_TEST const& other ) ;

      virtual void process_one_test( PEL_ModuleExplorer const* exp ) ;

   //-- Internals
      
      void check_diF( CH_BulkEnergy const* energy, 
                       double c1, double c2, double c3, bool& ok ) const ;
      
      void check_F_DDiF( CH_BulkEnergy const* energy,
                         double c1, double c2, bool& ok ) const ;
      
      void check_dj_DDiF( CH_BulkEnergy const* energy, 
                          double c1, double c2, bool& ok ) const ;
      
      void check_dj_ddiF( CH_BulkEnergy const* energy, 
                          double c1, double c2, bool& ok ) const ;

      void display_error( std::string const& mesg,
                          double xx_1, double xx_2 ) const ;

   //-- Class attributes

      static CH_BulkEnergy_TEST const* PROTOTYPE ;

   //-- Attributes

      double D_EPS ;
      double D_MIN ;
      double THICKNESS ;
      double S1 ;
      double S2 ;
      double S3 ;
      double ST ;
      double HH ;

} ;

#endif
