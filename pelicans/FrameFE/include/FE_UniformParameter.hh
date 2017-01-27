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

#ifndef FE_UNIFORM_PARAMETER_HH
#define FE_UNIFORM_PARAMETER_HH

#include <FE_Parameter.hh>

#include <doubleVector.hh>

/*
PUBLISHED
*/

class PEL_EXPORT FE_UniformParameter : public FE_Parameter
{
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static FE_UniformParameter* create( PEL_Object* a_owner,
		   		          std::string const& a_name,
				          doubleVector const& values ) ;

      virtual void reset( PEL_ModuleExplorer const* exp=0 ) ;

   //-- Instance characteristics

      virtual size_t nb_components( void ) const ;

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

    //-- Persistence

      virtual void save_state( PEL_ObjectWriter* writer ) const ;

      virtual void restore_state( PEL_ObjectReader* reader ) ;
      
   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

     ~FE_UniformParameter( void ) ;
      FE_UniformParameter( FE_UniformParameter const& other ) ;
      FE_UniformParameter& operator=( FE_UniformParameter const& other ) ;

      FE_UniformParameter( PEL_Object* a_owner,
                           PDE_DomainAndFields const* dom,
                           PEL_ModuleExplorer const* exp ) ;

      FE_UniformParameter( PEL_Object* a_owner,
                           std::string const& a_name,
                           doubleVector const& values ) ;

   //-- Plug in

      FE_UniformParameter( void ) ;

      virtual FE_UniformParameter* create_replica( 
                                 PEL_Object* a_owner,
                                 PDE_DomainAndFields const* dom,
                                 PEL_ModuleExplorer const* exp  ) const ;

   //-- Class attributes

      static FE_UniformParameter const* PROTOTYPE ;

   //-- Attributes

      doubleVector VAL ; 
} ;

#endif
