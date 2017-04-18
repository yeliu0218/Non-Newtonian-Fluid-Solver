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

#ifndef CH_BULK_CHEMICAL_POTENTIAL_HH
#define CH_BULK_CHEMICAL_POTENTIAL_HH

#include <PEL_Object.hh>

class PEL_ModuleExplorer ;

class CH_BulkEnergy ;

/*
PUBLISHED
*/

class CH_BulkChemicalPotential : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static CH_BulkChemicalPotential* create( 
                                       PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp  ) ;

   //-- Bulk energy

      double Sigma1( void ) const ;

      double Sigma2( void ) const ;

      double Sigma3( void ) const ;
      
      double F( double c1, double c2, double c3 ) const ;
      
   //-- Bulk chemical potential

      double DDiF( double c1, double c2, 
                   double c1_exp, double c2_exp, 
                   size_t i, double eps ) const ;
      
      double dj_DDiF( double c1, double c2, 
                      double c1_exp, double c2_exp, 
                      size_t i, size_t j, double eps ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      CH_BulkChemicalPotential( void ) ;
     ~CH_BulkChemicalPotential( void ) ;
      CH_BulkChemicalPotential( CH_BulkChemicalPotential const& other ) ;
      CH_BulkChemicalPotential& operator=( 
                                CH_BulkChemicalPotential const& other ) ;

      CH_BulkChemicalPotential( PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) ;
      
   //-- Attributes

      double S1 ;
      double S2 ;
      double S3 ;
      double ST ;
      CH_BulkEnergy const* F0 ;
      CH_BulkEnergy const* P ;
} ;

#endif
