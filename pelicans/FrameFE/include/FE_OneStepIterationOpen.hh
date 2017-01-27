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

#ifndef FE_ONE_STEP_ITERATION_OPEN_HH
#define FE_ONE_STEP_ITERATION_OPEN_HH

#include <FE_OneStepIteration.hh>

#include <map>
#include <string>

class LA_Matrix ;
class LA_SeqMatrix ;
class LA_SeqVector ;
class LA_Vector ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LinkDOF2Unknown ;
class PDE_SystemNumbering ;

/*
FRAMEWORK INSTANTIATION
PUBLISHED
*/

class PEL_EXPORT FE_OneStepIterationOpen : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      static FE_OneStepIterationOpen* object( std::string const& a_name ) ;
      
      static size_t nb_objects( void ) ;
      
      static FE_OneStepIterationOpen* object( size_t i ) ;

   //-- Instance characteristics

      std::string const& name( void ) const ;
      
      PDE_DomainAndFields const* domain( void ) const ;

   //-- Unknowns

      virtual size_t nb_unknowns( void ) const = 0 ;

      virtual PDE_DiscreteField* field( size_t i_unk ) const = 0 ;
      
      size_t index_of_field( std::string const& fname ) const ;

      virtual size_t level_of_field( size_t i_unk ) const = 0 ;

      virtual PDE_LinkDOF2Unknown const* link_DOF_2_unknown( 
                                              size_t i_unk ) const = 0 ;
      
   //-- Jacobian

      virtual void build_function_and_jacobian( FE_TimeIterator const* t_it ) ;

      virtual LA_SeqVector const* create_function( PEL_Object* a_owner,
                                                   size_t i_unk ) const ;

      virtual LA_SeqMatrix const* create_jacobian( PEL_Object* a_owner,
                                                   size_t i_eq, 
                                                   size_t j_unk ) const ;

   //-- Discretization
      
      virtual void assemble_contribution( FE_TimeIterator const* t_it,
                                          LA_Matrix* matrix,
                                          LA_Vector* vector,
                                          PDE_SystemNumbering const* nmb,
                                          size_t i_link ) const ;

      virtual void update_DOFs( LA_Vector* vector,
                                PDE_SystemNumbering const* nmb,
                                size_t i_link ) const ;
      
   protected: //--------------------------------------------------------------

      virtual ~FE_OneStepIterationOpen( void ) ;

      FE_OneStepIterationOpen( PEL_Object* a_owner, 
                               PDE_DomainAndFields const* dom,
                               PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      FE_OneStepIterationOpen( std::string const& a_name ) ;
      
   //-- Errors
      
      void check_assembled_system( PDE_SystemNumbering const* nmb,
                                   std::string const& function_name ) const ;

   //-- Preconditions, Postconditions, Invariant

      bool field_PRE( size_t i_unk ) const ;
      bool field_POST( PDE_DiscreteField const* result,
                       size_t i_unk ) const ;

      bool level_of_field_PRE( size_t i_unk ) const;
      bool level_of_field_POST( size_t result, size_t i_unk ) const ;

      bool link_DOF_2_unknown_PRE( size_t i_unk ) const ;
      bool link_DOF_2_unknown_POST( PDE_LinkDOF2Unknown const* result, 
                                    size_t i_unk ) const ;

      bool create_function_PRE( PEL_Object* a_owner, size_t i_unk ) const ;
      bool create_function_POST( LA_SeqVector const* result, 
                                 PEL_Object* a_owner,
                                 size_t i_unk ) const ;

      bool create_jacobian_PRE( PEL_Object* a_owner,
                                size_t i_eq, 
                                size_t j_unk ) const ;
      bool create_jacobian_POST( LA_SeqMatrix const* result,
                                 PEL_Object* a_owner,
                                 size_t i_eq, 
                                 size_t j_unk ) const ;
      
      bool assemble_contribution_PRE( FE_TimeIterator const* t_it,
                                      LA_Matrix* matrix,
                                      LA_Vector* vector,
                                      PDE_SystemNumbering const* nmb,
                                      size_t i_link ) const ;
      
      bool update_DOFs_PRE( LA_Vector* vector,
                            PDE_SystemNumbering const* nmb,
                            size_t i_link ) const ;

   private: //----------------------------------------------------------------

      FE_OneStepIterationOpen( void ) ;
      FE_OneStepIterationOpen( FE_OneStepIterationOpen const& other ) ;
      FE_OneStepIterationOpen& operator=( 
                               FE_OneStepIterationOpen const& other ) ;

   //-- Class attributes

      static std::map< std::string, FE_OneStepIterationOpen* > OBJS ;

   //-- Attributes

      std::string NN ;
      PDE_DomainAndFields const* DOM ;
} ;

#endif
