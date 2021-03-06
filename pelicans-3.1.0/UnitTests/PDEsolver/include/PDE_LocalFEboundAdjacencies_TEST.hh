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

#ifndef PDE_LOCAL_FE_BOUND_ADJACENCIES_TEST
#define PDE_LOCAL_FE_BOUND_ADJACENCIES_TEST

#include <PEL_ObjectTest.hh>

#include <queue>
#include <list>

class PEL_DataWithContext ;
class PEL_ModuleExplorer ;
class PEL_Vector ;
class PEL_ContextSimple ;
class PEL_DoubleVector ;

class GE_Color ;
class GE_Point ;
class GE_QRprovider ;

class PDE_DiscreteField ;
class PDE_LocalFEbound ;
class PDE_SetOfDomains ;

class PEL_EXPORT PDE_LocalFEboundAdjacencies_TEST : public PEL_ObjectTest
{
   public: //---------------------------------------------------------------

   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      PDE_LocalFEboundAdjacencies_TEST( void ) ;
     ~PDE_LocalFEboundAdjacencies_TEST( void ) ;
      PDE_LocalFEboundAdjacencies_TEST( 
                      PDE_LocalFEboundAdjacencies_TEST const& other ) ;
      PDE_LocalFEboundAdjacencies_TEST& operator=( 
                      PDE_LocalFEboundAdjacencies_TEST const& other ) ;

      virtual void process_one_test( PEL_ModuleExplorer const* exp ) ;

      void do_all_checks( PDE_SetOfDomains const* sdoms,
                          PEL_ModuleExplorer const* exp,
                          std::queue< double, std::list<double> >& valeurs,
                          bool learn ) ;

      void iterate_and_check_CP( 
                          PDE_LocalFEbound* fe,
                          GE_Color const* interf_color,
                          PDE_DiscreteField const* ff,
                          PEL_DataWithContext const* val,
                          PEL_DataWithContext const* jac,
                          PEL_DataWithContext const* hes,
                          std::queue< double, std::list<double> >& valeurs,
                          bool learn ) ;

      void iterate_and_check_IP( 
                          PDE_LocalFEbound* fe,
                          GE_Color const* interf_color,
                          PDE_DiscreteField const* ff,
                          PEL_DataWithContext const* val,
                          PEL_DataWithContext const* jac,
                          PEL_DataWithContext const* hes,
                          std::queue< double, std::list<double> >& valeurs,
                          bool learn ) ;

      void fill_context_with_coordinates( GE_Point const* pt ) ;

      void get_theo( double val, double& theo,
                     std::queue< double, std::list<double> >& valeurs,
                     bool learn ) const ;

      void display_error( std::string const& mesg,
                          std::string const& title_1,
                          std::string const& title_2,
                          double xx_1, double xx_2 ) const ;

   //-- Class attributes

      static PDE_LocalFEboundAdjacencies_TEST* registered_test ;

   //-- Attributes

      double D_EPS ;
      double D_MIN ;

      PEL_ContextSimple* CONTEXT ;
      PEL_DoubleVector* COORDS ;

      GE_QRprovider const* QRP ;
} ;

#endif
