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

#ifndef CH_MEAN_PARAMETER_HH
#define CH_MEAN_PARAMETER_HH

#include <FE_Parameter.hh>
#include <vector>

#include <doubleVector.hh>

class PDE_DiscreteField ;

/*
PUBLISHED
*/

class CH_MeanParameter : public FE_Parameter
{
   public: //---------------------------------------------------------------

   //-- Instance characteristics

      virtual size_t nb_components( void ) const ;

   //-- Values on cells

      virtual double cell_value_at_IP( FE_TimeIterator const* t_it,
                                       PDE_LocalFEcell const* fe,
                                       size_t ic = 0 ) const ;

      virtual double cell_value_at_pt( FE_TimeIterator const* t_it,
                                       PDE_LocalFEcell const* fe,
                                       size_t ic = 0 ) const;

      virtual double cell_gradient_at_IP( FE_TimeIterator const* t_it,
                                          PDE_LocalFEcell const* fe,
                                          size_t a,
                                          size_t ic = 0 ) const ;

   //-- Values on bounds
         
      virtual double bound_value_at_IP( FE_TimeIterator const* t_it,
                                       PDE_LocalFEbound const* fe,
                                       size_t ic = 0 ) const ;

   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

     ~CH_MeanParameter( void ) ;
      CH_MeanParameter( CH_MeanParameter const& other ) ;
      CH_MeanParameter& operator=( CH_MeanParameter const& other ) ;

      CH_MeanParameter( PEL_Object* a_owner,
                        PDE_DomainAndFields const* dom,
                        PEL_ModuleExplorer const* exp ) ;
      
   //-- Plug in

      CH_MeanParameter( void ) ;

      virtual CH_MeanParameter* create_replica( 
                                 PEL_Object* a_owner,
                                 PDE_DomainAndFields const* dom,
                                 PEL_ModuleExplorer const* exp  ) const ;

   //-- Prerequisites to calculations

      virtual void prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const ;
      virtual void prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const ;
      virtual void prepare_for_value_on_sides( PDE_CursorFEside* fe ) const ;
      virtual void prepare_for_gradient_on_cells( PDE_LocalFEcell* fe ) const ;

   //-- Class attributes

      static CH_MeanParameter const* PROTOTYPE ;

   //-- Internals

      double mean_value( void ) const ;
      
      double grad_mean_value( void ) const ;

      static double smoothed_delta( double x, double eps ) ;

      static double smoothed_Heavyside( double x, double eps ) ;

   //-- Attributes

      size_t NB_PHASES ;
      
      std::vector< PDE_DiscreteField const* > CCs ;
      std::vector< size_t > L_CCs ;
      std::vector< double > P_CCs ;
      
      mutable doubleVector PHs ;
      mutable doubleVector dPHs ;

      int type ;
} ;

#endif
