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

#ifndef FE_GRID_HIERARCHY_BUILDER_HH
#define FE_GRID_HIERARCHY_BUILDER_HH

#include <FE_OneStepIteration.hh>

class PDE_AdapterCHARMS ;

/*
HIGHLY UNSTABLE CLASS

PUBLISHED
*/

class PEL_EXPORT FE_GridHierarchyBuilder : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------

   //-- Substeps of the step by step progression

      virtual void do_before_time_stepping( FE_TimeIterator const* t_it ) ;

      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;
      
   //-- Elapsed times

      virtual void print_additional_times( std::ostream& os,
                                           size_t indent_width ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~FE_GridHierarchyBuilder( void ) ;
      FE_GridHierarchyBuilder( FE_GridHierarchyBuilder const& other ) ;
      FE_GridHierarchyBuilder& operator=( 
                               FE_GridHierarchyBuilder const& other ) ;

      FE_GridHierarchyBuilder( PEL_Object* a_owner, 
                               PDE_DomainAndFields const* dom,
                               FE_SetOfParameters const* prms,
                               PEL_ModuleExplorer const* exp ) ;

    //-- Plug in

      FE_GridHierarchyBuilder( void ) ;

      virtual FE_GridHierarchyBuilder* create_replica( 
	                         PEL_Object* a_owner,
	                         PDE_DomainAndFields const* dom,
	                         FE_SetOfParameters const* prms,
	                         PEL_ModuleExplorer* exp ) const ;

    //-- Internals

      void do_saving( void ) ;

    //-- Class attributes

      static FE_GridHierarchyBuilder const* PROTOTYPE ;
  
    //-- Attributes

      PDE_DomainAndFields const* DOM ;
      PDE_AdapterCHARMS* DA  ;
} ;

#endif
