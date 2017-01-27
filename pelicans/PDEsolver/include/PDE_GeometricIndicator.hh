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

#ifndef PDE_GEOMETRIC_INDICATOR_HH
#define PDE_GEOMETRIC_INDICATOR_HH

#include <PDE_AdaptationIndicator.hh>

#include <doubleVector.hh>

class GE_Point ;

class PEL_Context ;
class PEL_DataWithContext ;
class PEL_DoubleVector ;
class PEL_Int ;
class PEL_ModuleExplorer ;

class PEL_EXPORT PDE_GeometricIndicator : public PDE_AdaptationIndicator
{
   public: //-----------------------------------------------------------

   //-- Indicator calculation

      virtual void reset( void ) ;

      virtual void build( void ) ;

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

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~PDE_GeometricIndicator( void ) ;
      PDE_GeometricIndicator( PDE_GeometricIndicator const& other ) ;
      PDE_GeometricIndicator& operator=( 
                               PDE_GeometricIndicator const& other ) ;

      PDE_GeometricIndicator( PEL_Object* a_owner,
                               PDE_DomainAndFields const* dom,
                               PEL_ModuleExplorer const* exp,
                               size_t a_verbose_level ) ;

   //-- Plug in

      PDE_GeometricIndicator( void ) ;

      virtual PDE_GeometricIndicator* create_replica( 
                                              PEL_Object* a_owner,
                                              PDE_DomainAndFields const* dom,
                                              PEL_ModuleExplorer const* exp,
                                              size_t a_verbose_level ) const ;

   //-- Class attributes

      static PDE_GeometricIndicator const* PROTOTYPE ;

   //-- Attributes

      PEL_Context* CTX ;
      mutable PEL_DoubleVector* COORDS ;
      mutable PEL_Int* ITER ;
      PEL_DataWithContext* R_INDIC ;
      PEL_DataWithContext* U_INDIC ;
      size_t NB_REFS ;
      mutable GE_Point* PT ;
      size_t ICALL ;
      size_t ICALL_MAX ;
} ;

#endif
