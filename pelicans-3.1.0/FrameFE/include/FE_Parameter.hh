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

#ifndef FE_PARAMETER_HH
#define FE_PARAMETER_HH

#include <PEL_Object.hh>

#include <map>
#include <string>

class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;

class GE_Point ;

class PDE_CursorFEside ;
class PDE_DomainAndFields ;
class PDE_LocalFE ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_SetOfDomains ;

class FE_TimeIterator ;

class FE_SetOfParameters ;

/*
FRAMEWORK INSTANTIATION
PUBLISHED
*/

class PEL_EXPORT FE_Parameter : public PEL_Object
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      static FE_Parameter* make( PEL_Object* a_owner,
                                 PDE_DomainAndFields const* dom,
                                 PEL_ModuleExplorer const* exp ) ;

      static FE_Parameter* make( PEL_Object* a_owner,
                                 PDE_SetOfDomains const* sdoms,
                                 PEL_ModuleExplorer const* exp ) ;

      // IMPLEMENTATION : do nothing.
      virtual void do_the_links( FE_SetOfParameters const* prms ) ;

      // Reset internal datas (read the new values in `exp').
      // IMPLEMENTATION : raise an error.
      virtual void reset( PEL_ModuleExplorer const* exp = 0 ) ;

   //-- Instance characteristics

      std::string const& name( void ) const ;

      virtual size_t nb_components( void ) const = 0 ;

   //-- Prerequisites to calculations

      static int const Val = 1 ;

      static int const Grad = 2 ;

      void transfer_cell_calculation_requirements( PDE_LocalFEcell* fe,
                                                   int requisite ) ;

      bool ok_for_cell_calculations( PDE_LocalFEcell const* fe,
                                     int requisite ) const ;

      void transfer_side_calculation_requirements( PDE_CursorFEside* fe,
                                                   int requisite ) ;

      bool ok_for_side_calculations( PDE_CursorFEside const* fe,
                                     int requisite ) const ;

      void transfer_bound_calculation_requirements( PDE_LocalFEbound* fe,
                                                    int requisite ) ;
      
      bool ok_for_bound_calculations( PDE_LocalFEbound const* fe,
                                      int requisite ) const ;

   //-- Values on cells
      
      virtual double cell_value( FE_TimeIterator const* t_it,
                                 PDE_LocalFEcell const* fe,
                                 size_t ic = 0 ) const ;
      
      virtual double cell_gradient( FE_TimeIterator const* t_it,
                                    PDE_LocalFEcell const* fe,
                                    size_t a,
                                    size_t ic = 0 ) const ;
      
      virtual double cell_value_at_pt( FE_TimeIterator const* t_it,
                                       PDE_LocalFEcell const* fe,
                                       size_t ic = 0 ) const ;

      virtual double cell_gradient_at_pt( FE_TimeIterator const* t_it,
                                          PDE_LocalFEcell const* fe,
                                          size_t a,
                                          size_t ic = 0 ) const ;

      virtual double cell_value_at_IP( FE_TimeIterator const* t_it,
                                       PDE_LocalFEcell const* fe,
                                       size_t ic = 0 ) const ;

      virtual double cell_gradient_at_IP( FE_TimeIterator const* t_it,
                                          PDE_LocalFEcell const* fe,
                                          size_t a,
                                          size_t ic = 0 ) const ;
      
   //-- Values on bounds
      
     virtual double bound_value( FE_TimeIterator const* t_it,
                                 PDE_LocalFEbound const* fe,
                                 size_t ic = 0 ) const ;
      
     virtual double bound_gradient( FE_TimeIterator const* t_it,
                                    PDE_LocalFEbound const* fe,
                                    size_t a,
                                    size_t ic = 0 ) const ;

      virtual double bound_value_at_pt( FE_TimeIterator const* t_it,
                                        PDE_LocalFEbound const* fe,
                                        size_t ic = 0 ) const ;
      
      virtual double bound_gradient_at_pt( FE_TimeIterator const* t_it,
                                           PDE_LocalFEbound const* fe,
                                           size_t a,
                                           size_t ic = 0 ) const ;

      virtual double bound_value_at_IP( FE_TimeIterator const* t_it,
                                        PDE_LocalFEbound const* fe,
                                        size_t ic = 0 ) const ;

      virtual double bound_gradient_at_IP( FE_TimeIterator const* t_it,
                                           PDE_LocalFEbound const* fe,
                                           size_t a,
                                           size_t ic = 0 ) const ;

   //-- Values on sides
      
      virtual double side_value( FE_TimeIterator const* t_it,
                                 PDE_CursorFEside const* fe,
                                 size_t ic = 0 ) const ;
      
      virtual double side_gradient( FE_TimeIterator const* t_it,
                                    PDE_CursorFEside const* fe,
                                    size_t a,
                                    size_t ic = 0 ) const ;
      
      virtual double side_value_at_pt( FE_TimeIterator const* t_it,
                                       PDE_CursorFEside const* fe,
                                       size_t ic = 0 ) const ;
      
      virtual double side_gradient_at_pt( FE_TimeIterator const* t_it,
                                          PDE_CursorFEside const* fe,
                                          size_t a,
                                          size_t ic = 0 ) const ;

      virtual double side_value_at_IP( FE_TimeIterator const* t_it,
                                       PDE_CursorFEside const* fe,
                                       size_t ic = 0 ) const ;
      
      virtual double side_gradient_at_IP( FE_TimeIterator const* t_it,
                                          PDE_CursorFEside const* fe,
                                          size_t a,
                                          size_t ic = 0 ) const ;

   protected: //--------------------------------------------------------------

   //-- Plug in

      virtual ~FE_Parameter( void ) ;

      // Registration of the concrete subclass under the name `a_name'.
      FE_Parameter( std::string const& a_name ) ;

      FE_Parameter( PEL_Object* a_owner,
                    std::string const& a_name ) ;

      virtual FE_Parameter* create_replica( 
                                   PEL_Object* a_owner,
                                   PDE_DomainAndFields const* dom,
                                   PEL_ModuleExplorer const* exp ) const = 0 ;

      virtual FE_Parameter* create_replica( 
                                   PEL_Object* a_owner,
                                   PDE_SetOfDomains const* sdoms,
                                   PEL_ModuleExplorer const* exp ) const ;

      bool is_a_prototype( void ) const ;

   //-- Prerequisites to calculations

      virtual void prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const ;

      virtual void prepare_for_gradient_on_cells( PDE_LocalFEcell* fe ) const ;

      virtual void prepare_for_value_on_sides( PDE_CursorFEside* fe ) const ;

      virtual void prepare_for_gradient_on_sides( PDE_CursorFEside* fe ) const ;

      virtual void prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const ;

      virtual void prepare_for_gradient_on_bounds( 
                                                PDE_LocalFEbound* fe ) const ;

   //-- Preconditions, Postconditions, Invariant

      bool do_the_links_PRE( FE_SetOfParameters const* prms ) const ;
      
      virtual bool reset_PRE( PEL_ModuleExplorer const* exp ) const ;

      virtual bool prepare_for_value_on_cells_PRE(
                                  PDE_LocalFEcell const* fe ) const ;

      virtual bool prepare_for_gradient_on_cells_PRE(
                                  PDE_LocalFEcell const* fe ) const ;

      virtual bool prepare_for_value_on_sides_PRE(
                                  PDE_CursorFEside const * fe ) const ;

      virtual bool prepare_for_gradient_on_sides_PRE(
                                  PDE_CursorFEside const* fe ) const ;

      virtual bool prepare_for_value_on_bounds_PRE(
                                  PDE_LocalFEbound const* fe ) const ;

      virtual bool prepare_for_gradient_on_bounds_PRE(
                                  PDE_LocalFEbound const* fe ) const;
      
      virtual bool cell_value_PRE( FE_TimeIterator const* t_it,
                                   PDE_LocalFEcell const* fe,
                                   size_t ic ) const ;
      
      virtual bool cell_gradient_PRE( FE_TimeIterator const* t_it,
                                      PDE_LocalFEcell const* fe,
                                      size_t a,
                                      size_t ic ) const ;
      
      virtual bool cell_value_at_pt_PRE( FE_TimeIterator const* t_it,
                                         PDE_LocalFEcell const* fe,
                                         size_t ic ) const ;
      
      virtual bool cell_gradient_at_pt_PRE( FE_TimeIterator const* t_it,
                                            PDE_LocalFEcell const* fe,
                                            size_t a,
                                            size_t ic ) const ;

      virtual bool cell_value_at_IP_PRE( FE_TimeIterator const* t_it,
                                         PDE_LocalFEcell const* fe,
                                         size_t ic ) const ;

      virtual bool cell_gradient_at_IP_PRE( FE_TimeIterator const* t_it,
                                            PDE_LocalFEcell const* fe,
                                            size_t a,
                                            size_t ic ) const ;

      virtual bool bound_value_PRE( FE_TimeIterator const* t_it,
                                    PDE_LocalFEbound const* fe,
                                    size_t ic ) const ;
      
      virtual bool bound_gradient_PRE( FE_TimeIterator const* t_it,
                                       PDE_LocalFEbound const* fe,
                                       size_t a,
                                       size_t ic ) const ;
      
      virtual bool bound_value_at_pt_PRE( FE_TimeIterator const* t_it,
                                          PDE_LocalFEbound const* fe,
                                          size_t ic ) const ;
      
      virtual bool bound_gradient_at_pt_PRE( FE_TimeIterator const* t_it,
                                             PDE_LocalFEbound const* fe,
                                             size_t a,
                                             size_t ic ) const ;

      virtual bool bound_value_at_IP_PRE( FE_TimeIterator const* t_it,
                                          PDE_LocalFEbound const* fe,
                                          size_t ic ) const ;

      virtual bool bound_gradient_at_IP_PRE( FE_TimeIterator const* t_it,
                                             PDE_LocalFEbound const* fe,
                                             size_t a,
                                             size_t ic ) const ;

      virtual bool side_value_PRE( FE_TimeIterator const* t_it,
                                   PDE_CursorFEside const* fe,
                                   size_t ic ) const ;

      virtual bool side_gradient_PRE( FE_TimeIterator const* t_it,
                                      PDE_CursorFEside const* fe,
                                      size_t a,
                                      size_t ic ) const ;
      
      virtual bool side_value_at_pt_PRE( FE_TimeIterator const* t_it,
                                         PDE_CursorFEside const* fe,
                                         size_t ic ) const ;
      
      virtual bool side_gradient_at_pt_PRE( FE_TimeIterator const* t_it,
                                            PDE_CursorFEside const* fe,
                                            size_t a,
                                            size_t ic ) const ;

      virtual bool side_value_at_IP_PRE( FE_TimeIterator const* t_it,
                                         PDE_CursorFEside const* fe,
                                         size_t ic ) const ;
      
      virtual bool side_gradient_at_IP_PRE( FE_TimeIterator const* t_it,
                                            PDE_CursorFEside const* fe,
                                            size_t a,
                                            size_t ic ) const ;

      virtual bool create_replica_PRE( PEL_Object* a_owner,
                                       PDE_DomainAndFields const* dom,
                                       PEL_ModuleExplorer const* exp ) const ;

      virtual bool create_replica_POST( FE_Parameter const* result,
                                        PEL_Object* a_owner,
                                        PDE_DomainAndFields const* dom,
                                        PEL_ModuleExplorer const* exp ) const ;

      virtual bool create_replica_PRE( PEL_Object* a_owner,
                                       PDE_SetOfDomains const* sdoms,
                                       PEL_ModuleExplorer const* exp ) const ;

      virtual bool create_replica_POST( FE_Parameter const* result,
                                        PEL_Object* a_owner,
                                        PDE_SetOfDomains const* sdoms,
                                        PEL_ModuleExplorer const* exp ) const ;

   private: //----------------------------------------------------------------

      FE_Parameter( void ) ;
      FE_Parameter( FE_Parameter const& other ) ;
      FE_Parameter& operator=( FE_Parameter const& other ) ;

      static PEL_ObjectRegister* plugins_map( void ) ;
      
   //-- Attributes

      bool const IS_PROTO ;
      std::string const NAME ;
      std::map< PDE_LocalFEcell const*, bool > VAL_ON_CELLS ;
      std::map< PDE_LocalFEcell const*, bool > GRAD_ON_CELLS ;
      std::map< PDE_LocalFEbound const*, bool > VAL_ON_BOUNDS ;
      std::map< PDE_LocalFEbound const*, bool > GRAD_ON_BOUNDS ;
} ;

#endif
