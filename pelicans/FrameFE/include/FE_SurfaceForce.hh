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

#ifndef FE_SURFACE_FORCE_HH
#define FE_SURFACE_FORCE_HH

#include <FE_OneStepIteration.hh>

#include <doubleVector.hh>

#include <map>

class GE_QRprovider ;
class GE_Color ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LocalFEbound ;

class FE_Parameter ;

class FE_SurfaceForce : public FE_OneStepIteration
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      static FE_SurfaceForce* object( std::string const& a_name ) ;

   //-- Instance characteristics

      virtual std::string const& name( void ) const ;

   //-- Substeps of the step by step progression

      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

   //-- Savings for post-processing
 
      void compute_force( FE_TimeIterator const* t_it ) ;

      virtual void save_other_than_time_and_fields( 
                                             FE_TimeIterator const* t_it, 
                                             PDE_ResultSaver* rs ) ;

   //-- Results derived from the computations

      double force( size_t d ) const ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~FE_SurfaceForce( void ) ;
      FE_SurfaceForce( FE_SurfaceForce const& other ) ;
      FE_SurfaceForce& operator=( FE_SurfaceForce const& other ) ;

      FE_SurfaceForce( PEL_Object* a_owner,
                       PDE_DomainAndFields const* dom,
                       FE_SetOfParameters const* prms,
                       PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      FE_SurfaceForce( void ) ;

      virtual FE_SurfaceForce* create_replica( 
                                   PEL_Object* a_owner,
                                   PDE_DomainAndFields const* dom,
                                   FE_SetOfParameters const* prms,
                                   PEL_ModuleExplorer* exp ) const ;

   //-- Class attributes

      static std::map< std::string, FE_SurfaceForce* > OBJS ;

      static FE_SurfaceForce const* PROTOTYPE ;

   //-- Attributes

      std::string NAME ;

      PDE_DiscreteField const* UU ;
      size_t L_UU ;
      PDE_DiscreteField const* PP ;
      size_t L_PP ;

      PDE_LocalFEbound* bFE ;
      GE_QRprovider const* QRP ;

      FE_Parameter* MU ;
      GE_Color const* SURFCOL ;

      doubleVector FORCE ;
      doubleVector FORCE_P ;
      doubleVector FORCE_VISC ;
} ;

#endif
