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

#ifndef CH_RESIDUAL_NS_CELL_INDICATOR_HH
#define CH_RESIDUAL_NS_CELL_INDICATOR_HH

#include <FE_AdaptationIndicator.hh>

#include <doubleArray2D.hh>

class GE_Point ;
class GE_QRprovider ;

class PDE_CursorFEside ;
class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_SetOfBCs ;

class FE_SetOfParameters ;
class FE_Parameter ;
class FE_TimeIterator ;

class CH_ResidualNSCellIndicator : public FE_AdaptationIndicator
{
   public: //-----------------------------------------------------------
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~CH_ResidualNSCellIndicator( void ) ;
      CH_ResidualNSCellIndicator( CH_ResidualNSCellIndicator const& other ) ;
      CH_ResidualNSCellIndicator& operator=( 
                                  CH_ResidualNSCellIndicator const& other ) ;

      CH_ResidualNSCellIndicator( PEL_Object* a_owner,
                                  PDE_DomainAndFields const* dom,
                                  FE_SetOfParameters const* prms,
                                  PEL_ModuleExplorer const* exp,
                                  size_t a_verbose_level ) ;

   //-- Plug in

      CH_ResidualNSCellIndicator( void ) ;

      virtual CH_ResidualNSCellIndicator* create_replica( 
                                            PEL_Object* a_owner,
                                            PDE_DomainAndFields const* dom,
                                            FE_SetOfParameters const* prms,
                                            PEL_ModuleExplorer const* exp,
                                            size_t a_verbose_level ) const ;

   //-- Indicator calculation

      virtual void reset_self( void ) ;

      virtual void build_self( FE_TimeIterator const* t_it ) ;

   //-- Indicator results

      virtual double cell_indicator( size_t cell_id ) const ;

      virtual bool to_be_refined( double bf_indicator,
                                  GE_Mpolyhedron const* poly,
                                  PDE_ReferenceElement const* elm,
                                  size_t local_node ) const ;

      virtual bool to_be_unrefined( double bf_indicator,
                                    GE_Mpolyhedron const* poly,
                                    PDE_ReferenceElement const* elm,
                                    size_t local_node ) const ;

   //-- Class attributes

      static CH_ResidualNSCellIndicator const* PROTOTYPE ;

   //-- Attributes

      size_t NB_REFS ;
      mutable GE_Point* PT ;
      size_t ICALL ;

      PDE_DiscreteField const* UU ;
      PDE_DiscreteField const* PP ; 

      size_t L_UU ; 
      size_t L_PP ;
      
      FE_Parameter* DENS ;
      FE_Parameter* MU ;
      FE_Parameter* RHS ;
      
      PDE_SetOfBCs const* BCs ;
      PDE_CursorFEside* sFE ;
      PDE_LocalFEbound* bFE ;
      PDE_LocalFEcell* cFE ;

      GE_QRprovider const* QRP ;
      size_t NB_DIMS ;
      doubleVector CELL_ERRORS ;
      double CELL_ERRORS_MAX ;
      double MAX_ERR ;
      double MIN_ERR ;
      bool BUILD_OK ;
} ;

#endif
