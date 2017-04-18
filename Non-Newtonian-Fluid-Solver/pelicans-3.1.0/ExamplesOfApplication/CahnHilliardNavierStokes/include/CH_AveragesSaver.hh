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

#ifndef CH_AVERAGES_SAVER_HH
#define CH_AVERAGES_SAVER_HH

#include <FE_OneStepIteration.hh>

#include <doubleVector.hh>

class doubleArray2D ;
class size_t_vector ;

class GE_QRprovider ;

class PDE_DiscreteField ;
class PDE_LocalFEcell ;

class CH_BulkChemicalPotential ;

class CH_AveragesSaver : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression

      virtual void do_before_time_stepping( FE_TimeIterator const* t_it ) ;
      
      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

      virtual void do_after_time_adaptation( FE_TimeIterator const* t_it ) ;
      
   //-- Savings for post-processing

      virtual void save_other_than_time_and_fields( 
                                            FE_TimeIterator const* t_it,
                                            PDE_ResultSaver* rs ) ;
   //-- Available averages

      // Total free energy
      bool total_free_energy_is_handled( void ) const;
      
      double total_free_energy( void ) const ;
           
      // Kinetic energy
      bool kinetic_energy_is_handled( void ) const ;
      
      double kinetic_energy( FE_TimeIterator const* t_it ) const ;

      // Volume and coordinates of center
      bool center_volume_is_handled( void ) const ;
      
      void center_volume( doubleVector& volume,
                          doubleArray2D& center,
                          doubleArray2D& velocity,
                          size_t_vector& nb_cells ) const ;

      // Perimeter
      bool perimeter_is_handled( void ) const ;
      
      double perimeter( void ) const ;
      
   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~CH_AveragesSaver( void ) ;
      CH_AveragesSaver( CH_AveragesSaver const& other ) ;
      CH_AveragesSaver& operator=( CH_AveragesSaver const& other ) ;

      CH_AveragesSaver( PEL_Object* a_owner,
                        PDE_DomainAndFields const* dom,
                        FE_SetOfParameters const* prms,
                        PEL_ModuleExplorer* exp ) ;

   //-- Plug in

      CH_AveragesSaver( void ) ;

      virtual CH_AveragesSaver* create_replica(  
                                     PEL_Object* a_owner,
                                     PDE_DomainAndFields const* dom,
                                     FE_SetOfParameters const* prms,
                                     PEL_ModuleExplorer* exp ) const ;
      
   //-- Internals
      
      void read_for_total_free_energy( PDE_DomainAndFields const* dom,
                                       FE_SetOfParameters const* prms,
                                       PEL_ModuleExplorer const* exp ) ;
      
      void read_for_kinetic_energy( PDE_DomainAndFields const* dom,
                                    FE_SetOfParameters const* prms,
                                    PEL_ModuleExplorer const* exp ) ;
      
      void read_for_center_volume( PDE_DomainAndFields const* dom,
                                   FE_SetOfParameters const* prms,
                                   PEL_ModuleExplorer const* exp ) ;
      
      void read_for_perimeter( PDE_DomainAndFields const* dom,
                               FE_SetOfParameters const* prms,
                               PEL_ModuleExplorer const* exp ) ;
      
      void read_for_velocity( PDE_DomainAndFields const* dom,
                              FE_SetOfParameters const* prms,
                              PEL_ModuleExplorer const* exp ) ;
      
      void save_all( FE_TimeIterator const* t_it ) ;

   //-- Class attributes

      static CH_AveragesSaver const* PROTOTYPE ;

   //-- Attributes
      
      std::string FNAME ;
      PDE_LocalFEcell* cFE ;
      
      bool SOMETHING_TO_DO ;

      bool TOTAL_FREE_ENERGY ;
      double THICKNESS ;
      CH_BulkChemicalPotential* BULK_MU ;
      PDE_DiscreteField const* C1 ;
      size_t L_C1 ;
      PDE_DiscreteField const* C2 ;
      size_t L_C2 ;
      GE_QRprovider const* QRP_TFE ;
      
      bool KINETIC_ENERGY ;
      PDE_DiscreteField const* UU;
      size_t L_UU ;
      FE_Parameter* DENS ;
      GE_QRprovider const* QRP_KE ;
      
      bool CENTER ;
      PDE_DiscreteField const* CC ;
      size_t L_CC ;
      PDE_DiscreteField const* VV ;
      size_t L_VV ;
      doubleVector CC_MIN ;
      GE_QRprovider const* QRP_CV ;    
      
      bool PERIMETER ;
      PDE_DiscreteField const* CP ;
      size_t L_CP ;
      GE_QRprovider const* QRP_PER ;    

      double NEXT_SAVING_TIME ;
      doubleVector SAVING_TIMES ;
} ;

#endif
