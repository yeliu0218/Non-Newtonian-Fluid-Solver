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

#ifndef CH_INTERFACE_INDICATOR_HH
#define CH_INTERFACE_INDICATOR_HH

#include <PDE_AdaptationIndicator.hh>

#include <doubleVector.hh>

class GE_QRprovider ;

class PDE_DiscreteField ;
class PDE_LocalFEcell ;

class CH_InterfaceIndicator : public PDE_AdaptationIndicator
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

     ~CH_InterfaceIndicator( void ) ;
      CH_InterfaceIndicator( CH_InterfaceIndicator const& other ) ;
      CH_InterfaceIndicator& operator=( 
                               CH_InterfaceIndicator const& other ) ;

      CH_InterfaceIndicator( PEL_Object* a_owner,
                               PDE_DomainAndFields const* dom,
                               PEL_ModuleExplorer const* exp,
                               size_t a_verbose_level ) ;

   //-- Plug in

      CH_InterfaceIndicator( void ) ;

      virtual CH_InterfaceIndicator* create_replica( 
                                            PEL_Object* a_owner,
                                            PDE_DomainAndFields const* dom,
                                            PEL_ModuleExplorer const* exp,
                                            size_t a_verbose_level ) const ;

   //-- Class attributes

      static CH_InterfaceIndicator const* PROTOTYPE ;

   //-- Attributes

      PDE_DiscreteField const* C1 ;
      PDE_DiscreteField const* C2 ;
      size_t LL ;
      PDE_LocalFEcell* cFE ;
      GE_QRprovider const* QRP ;
      doubleVector CELL_INTERF ;
      double H_INTERF ;
      double VAL_REFI ;
      double VAL_UNREFI ;
} ;

#endif
