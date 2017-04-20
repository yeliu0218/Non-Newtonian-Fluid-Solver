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

#ifndef FE_STEADY_STATE_ADAPTER_HH
#define FE_STEADY_STATE_ADAPTER_HH

#include <FE_TimeIteratorAdapter.hh>

#include <doubleVector.hh>
#include <stringVector.hh>

class PDE_DiscreteField ;
class PDE_SetOfDiscreteFields ;

/*
Time adapter which stop calculation when steady state is reached.

Example:

   MODULE PEL_Application
      ...
      MODULE FE_TimeIterator
         ...
      END MODULE FE_TimeIterator
      MODULE FE_TimeIteratorAdapter
         concrete_name = "FE_SteadyStateAdapter"
         initial_level = 1
         current_level = 0
         discrete_fields = < "uu" "pp" >
         minimal_error = < 1.E-2 5.E-2 >
         MODULE post_processing
            banner = true
            file_name = "convergence.txt"
         END MODULE post_processing
      END MODULE FE_TimeIteratorAdapter
      ...
   END MODULE PEL_Application

   The time iterations are stopped when:
      -     || uu(0)-uu(1) ||/ || uu(1), 1. || < 1.E-2
        and || pp(0)-pp(1) ||/ || pp(1), 1. || < 5.E-2
      - or when final time defined in FE_TimeIterator module is reached.

   Rem : ||.|| is the L_infinity norm.

   The errors of uu and pp are stored at each time step in the
   output file "convergence.txt".
   
PUBLISHED
*/

class PEL_EXPORT FE_SteadyStateAdapter : public FE_TimeIteratorAdapter
{

   public: //-----------------------------------------------------------
      
   protected: //--------------------------------------------------------
      
   private: //----------------------------------------------------------

      FE_SteadyStateAdapter( void ) ;
     ~FE_SteadyStateAdapter( void ) ;
      FE_SteadyStateAdapter( FE_SteadyStateAdapter const& other ) ;
      FE_SteadyStateAdapter& operator=(
                             FE_SteadyStateAdapter const& other ) ;

      FE_SteadyStateAdapter( FE_TimeIterator* t_it,
                             PDE_DomainAndFields const* dom,
                             FE_SetOfParameters const* prms,
                             PEL_ModuleExplorer const* exp ) ;
      
   //-- Plug in

      FE_SteadyStateAdapter( std::string const& concrete_name ) ;
      
      virtual FE_TimeIteratorAdapter* create_replica( 
                                     FE_TimeIterator* t_it,
                                     PDE_DomainAndFields const* dom,
                                     FE_SetOfParameters const* prms,
                                     PEL_ModuleExplorer const* exp ) const ;
      
   //-- Adaptation of the associated FE_TimeIterator object

      // IMPLEMENTATION : check steady state 
      virtual void define_parameters_for_next_iteration(
                          bool& finished, bool& restart, double& next_dt ) ;
      
      void initialize_cv_file( bool banner ) const ;
      void save_in_cv_file( doubleVector const& values ) const ;
      
   //-- Static attribute

      static FE_SteadyStateAdapter const* PROTOTYPE ;
      
   //-- Attribute

      PDE_SetOfDiscreteFields const* const FIELDS ;
      size_t const LEVEL1 ;
      size_t const LEVEL0 ;
      stringVector FIELDS_TABLE ;
      doubleVector FIELDS_ERROR ;
      std::string OFILE_NAME ;
} ;

#endif
