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

#ifndef PDE_DOF_S_INITIALIZER_HH
#define PDE_DOF_S_INITIALIZER_HH

#include <PEL_Object.hh>

class PEL_Context ;
class PEL_DoubleVector ;
class PEL_ModuleExplorer ;
class PEL_Vector ;
class boolArray2D ;
class doubleArray2D ;
class doubleVector ;

class GE_Point ;
class GE_SetOfPoints ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_InterfaceAndFields ;
class PDE_LocalFE ;
class PDE_ResultReader ;

/*
Initializers of `::PDE_DiscreteField' objects in case of finite element
discretization.
*/

class PEL_EXPORT PDE_DOFsInitializer : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_DOFsInitializer const* object( void ) ;

   //-- Discrete fields initialization

      // Set the values and the imposed values of the DOFs of the discrete
      // fields of `df' linked to `df_exp'.
      void initialize_discrete_fields(
                              PDE_DomainAndFields const* df,
                              PEL_ModuleExplorer const* df_exp ) const ;

      void initialize_discrete_fields(
                              PDE_InterfaceAndFields const* interf,
                              PEL_ModuleExplorer const* interf_exp ) const ;

      void set_discrete_fields_DOFs_values(
                              PDE_DomainAndFields const* df,
                              PEL_ModuleExplorer const* df_exp,
                              bool verbose ) const ;

      void set_discrete_fields_imposed_DOFs(
                              PDE_DomainAndFields const* df,
                              PEL_ModuleExplorer const* df_exp,
                              bool verbose ) const ;

      void set_discrete_fields_imposed_DOFs(
                              PDE_InterfaceAndFields const* interf,
                              PEL_ModuleExplorer const* interf_exp,
                              bool verbose ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_DOFsInitializer( void ) ;
     ~PDE_DOFsInitializer( void ) ;
      PDE_DOFsInitializer( PDE_DOFsInitializer const& other ) ;
      PDE_DOFsInitializer& operator=( PDE_DOFsInitializer const& other ) ;
      
   //-- Discrete fields initialization

      void apply_DOFs_values( PDE_DiscreteField* f,
                              PDE_DomainAndFields const* df,
                              PEL_ModuleExplorer const* df_exp,
                              PEL_ModuleExplorer const* f_exp,
                              PEL_Vector const* readers,
                              bool is_boundary_field,
                              bool verbose ) const ;

      void apply_DOFs_imposed_value( PDE_DiscreteField* f,
                                     PDE_DomainAndFields const* df,
                                     PEL_ModuleExplorer const* f_exp,
                                     bool is_boundary_field,
                                     bool verbose ) const ;

      void apply_DOFs_imposed_value( PDE_DiscreteField* f,
                                     PDE_InterfaceAndFields const* df,
                                     PEL_ModuleExplorer const* f_exp,
                                     bool verbose ) const ;

      void impose_on_bounds( PDE_DiscreteField* f,
                             PDE_DomainAndFields const* df,
                             PEL_ModuleExplorer const* exp,
                             bool verbose,
                             boolArray2D& already_done ) const ;

      void impose_in_region( PDE_DiscreteField* f,
                             GE_SetOfPoints const* vertices,
                             PDE_LocalFE* fe,
                             PEL_ModuleExplorer const* exp,
                             bool verbose,
                             boolArray2D& already_done ) const ;

      static void read_components( PEL_ModuleExplorer const* exp,
                                   PDE_DiscreteField const* ff,
                                   size_t& nbc, int& comp ) ;

      static void modify_field_DOFs( PDE_DiscreteField* ff,
                                     size_t node,
                                     GE_Point const* pt_node,
                                     int comp,
                                     doubleVector const& val,
                                     boolArray2D& already_done ) ;

      void initialize_by_L2_projection( PDE_DiscreteField* f,
                                        PDE_DomainAndFields const* df,
                                        PEL_ModuleExplorer const* df_exp,
                                        PEL_ModuleExplorer const* exp,
                                        PDE_ResultReader const* reader,
                                        bool is_boundary_field,
                                        bool verbose ) const ;
      
      void initialize_by_value_at_node_location(
                                        PDE_DiscreteField* f,
                                        PDE_DomainAndFields const* df,
                                        PEL_ModuleExplorer const* df_exp,
                                        PEL_ModuleExplorer const* exp,
                                        PDE_ResultReader const* reader,
                                        bool is_boundary_field,
                                        bool verbose ) const ;

   //-- Attributes

      PEL_Context* CTX ;
      PEL_DoubleVector* COORDS ;
} ;


#endif
