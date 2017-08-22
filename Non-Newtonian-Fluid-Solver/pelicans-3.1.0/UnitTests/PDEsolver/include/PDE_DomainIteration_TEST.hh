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

#ifndef PDE_DOMAIN_ITERATION_TEST
#define PDE_DOMAIN_ITERATION_TEST

#include <PEL_ObjectTest.hh>

#include <iosfwd>

class GE_QRprovider ;

class PDE_CursorFEside ;
class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_InterfaceAndFields ;
class PDE_LocalFE ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_SetOfDiscreteFields ;

class PEL_EXPORT PDE_DomainIteration_TEST : public PEL_ObjectTest
{
   public: //---------------------------------------------------------------

   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      PDE_DomainIteration_TEST( void ) ;
     ~PDE_DomainIteration_TEST( void ) ;
      PDE_DomainIteration_TEST( PDE_DomainIteration_TEST const& other ) ;
      PDE_DomainIteration_TEST& operator=( 
                                PDE_DomainIteration_TEST const& other ) ;

      virtual void process_one_test( PEL_ModuleExplorer const* exp ) ;

      void process_domain( PDE_DomainAndFields const* dom,
                           std::ostream& os,
                           bool do_print_IPs,
                           bool do_print_vals_at_IP ) ;

      void process_interface( PDE_InterfaceAndFields const* interf,
                              PDE_DomainAndFields const* dom_0,
                              PDE_DomainAndFields const* dom_1,
                              std::ostream& os,
                              bool do_print_IPs,
                              bool do_print_vals_at_IPs ) ;

      void iterate_over_meshes( PDE_LocalFE* fe_1,
                                PDE_LocalFE* fe_2,
                                std::ostream& os,
                                bool do_print_IPs,
                                bool do_print_vals_at_IPs ) ;

      void iterate_over_cells( PDE_LocalFEcell*  cfe,
                               PDE_CursorFEside* sfe,
                               PDE_LocalFEbound* bfe ) ;
      
      void check_node_location( PDE_LocalFE* fe,
                                bool& ok_node_location ) const ;

      void iterate_over_IPs( PDE_LocalFE* fe, 
                             bool& ok_val, 
                             bool& done_der1, bool& ok_der1, 
                             bool& done_der2, bool& ok_der2,
                             std::ostream& os,
                             bool do_print_IPs,
                             bool do_print_vals_at_IPs ) const ;

      void check_IP_vals( PDE_LocalFE* fe, 
                          PDE_DiscreteField const* ff, 
                          bool& ok ) const ;
 
      void check_IP_der1( PDE_LocalFE* fe, 
                          PDE_DiscreteField const* ff, 
                          bool& ok ) const ;
 
      void check_IP_der2( PDE_LocalFE* fe, 
                          PDE_DiscreteField const* ff, 
                          bool& ok ) const ;

      void iterate_over_sides( PDE_CursorFEside* sfe_1,
                               PDE_CursorFEside* sfe_2,
                               PDE_DomainAndFields const* dom,
                               std::ostream& os ) ;
      
      void check_periodicity( PDE_CursorFEside* sfe, bool& ok_perio ) ;

      void print_adjacent_cell_nodes( size_t i_adj,
                                      PDE_CursorFEside* sfe_1,
                                      PDE_CursorFEside* sfe_2,
                                      PDE_DiscreteField const* ff,
                                      std::ostream& os,
                                      bool& ok ) const ;

      void trace_IPs(  PDE_LocalFE* fe, std::ostream& os ) const ;

      void display_error( std::string const& mesg,
                          std::string const& title_1,
                          std::string const& title_2,
                          double xx_1, double xx_2 ) const ;

   //-- Class attributes

      static PDE_DomainIteration_TEST* registered_test ;

   //-- Attributes

      double D_EPS ;
      double D_MIN ;
      GE_QRprovider const* QRP ;
} ;

#endif
