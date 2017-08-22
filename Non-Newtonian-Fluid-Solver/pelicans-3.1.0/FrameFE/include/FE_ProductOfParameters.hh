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

#ifndef FE_PRODUCT_OF_PARAMETERS_HH
#define FE_PRODUCT_OF_PARAMETERS_HH

#include <FE_Parameter.hh>

#include <vector>

/*
Parameters the value of which is a product of `FE_Parameter::'.

PUBLISHED
*/

class PEL_EXPORT FE_ProductOfParameters : public FE_Parameter
{
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      virtual void do_the_links( FE_SetOfParameters const* prms ) ;
      
      virtual void reset( PEL_ModuleExplorer const* exp = 0 ) ;

   //-- Instance characteristics

      virtual size_t nb_components( void ) const ;


   //-- Values on cells
      
      virtual double cell_value( FE_TimeIterator const* t_it,
                                 PDE_LocalFEcell const* fe,
                                 size_t ic = 0 ) const ;
      
      virtual double cell_value_at_pt( FE_TimeIterator const* t_it,
                                       PDE_LocalFEcell const* fe,
                                       size_t ic = 0 ) const ;

      virtual double cell_value_at_IP( FE_TimeIterator const* t_it,
                                       PDE_LocalFEcell const* fe,
                                       size_t ic = 0 ) const ;

   //-- Values on bounds
      
     virtual double bound_value( FE_TimeIterator const* t_it,
                                 PDE_LocalFEbound const* fe,
                                 size_t ic = 0 ) const ;
      
      virtual double bound_value_at_pt( FE_TimeIterator const* t_it,
                                        PDE_LocalFEbound const* fe,
                                        size_t ic = 0 ) const ;
      
      virtual double bound_value_at_IP( FE_TimeIterator const* t_it,
                                        PDE_LocalFEbound const* fe,
                                        size_t ic = 0 ) const ;

   //-- Values on sides
      
      virtual double side_value( FE_TimeIterator const* t_it,
                                 PDE_CursorFEside const* fe,
                                 size_t ic = 0 ) const ;
      
      virtual double side_value_at_pt( FE_TimeIterator const* t_it,
                                       PDE_CursorFEside const* fe,
                                       size_t ic = 0 ) const ;
      
      virtual double side_value_at_IP( FE_TimeIterator const* t_it,
                                       PDE_CursorFEside const* fe,
                                       size_t ic = 0 ) const ;
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

     ~FE_ProductOfParameters( void ) ;
      FE_ProductOfParameters( FE_ProductOfParameters const& other ) ;
      FE_ProductOfParameters& operator=( FE_ProductOfParameters const& other ) ;

      FE_ProductOfParameters( PEL_Object* a_owner,
                              PDE_DomainAndFields const* dom,
                              PEL_ModuleExplorer const* exp ) ;
   //-- Plug in

      FE_ProductOfParameters( void ) ;

      virtual FE_ProductOfParameters* create_replica( 
                                      PEL_Object* a_owner,
                                      PDE_DomainAndFields const* dom,
                                      PEL_ModuleExplorer const* exp  ) const ;

   //-- Prerequisites to calculations

      virtual void prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const ;
      virtual void prepare_for_value_on_sides( PDE_CursorFEside* fe ) const ;
      virtual void prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const ;
      
   //-- Internals
      
      void infer_nb_components( FE_Parameter const* prm ) ;

   //-- Class attributes

      static FE_ProductOfParameters const* PROTOTYPE ;

   //-- Attributes

      size_t NBCS ;
      std::vector< FE_Parameter* > PRMS ;
      std::vector< std::string > NAMES ;
} ;

#endif