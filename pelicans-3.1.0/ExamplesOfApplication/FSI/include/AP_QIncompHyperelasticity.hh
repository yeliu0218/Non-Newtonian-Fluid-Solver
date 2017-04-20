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

#ifndef AP_Q_INCOMP_HYPERELASTICITY_HH
#define AP_Q_INCOMP_HYPERELASTICITY_HH

#include <FE_OneStepIterationOpen.hh>

#include <doubleVector.hh>

class doubleArray2D ;
class doubleArray4D ;

class LA_Matrix ;
class LA_SeqVector ;
class LA_Solver ;
class LA_Vector ;

class GE_QRprovider ;

class PDE_DiscreteField ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalEquation ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_SetOfBCs ;

class FE_Parameter ;
class FE_SetOfParameters ;

class AP_ConstitutiveLaw ;
class AP_KinematicState ;
class AP_LoadCalculator ;
class AP_StructureSystem ;

/*
PUBLISHED
*/

class AP_QIncompHyperelasticity : public FE_OneStepIterationOpen
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression
      
      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

   //-- Jacobian

      virtual size_t nb_unknowns( void ) const ;

      virtual PDE_DiscreteField* field( size_t i_unk ) const ;

      virtual size_t level_of_field( size_t i_unk ) const ;

      virtual PDE_LinkDOF2Unknown const* link_DOF_2_unknown( size_t i_unk ) const ;

      virtual void build_function_and_jacobian( FE_TimeIterator const* t_it ) ;

      virtual LA_SeqVector const* create_function( PEL_Object* a_owner,
                                                   size_t i_unk ) const ;

      virtual LA_SeqMatrix const* create_jacobian( PEL_Object* a_owner,
                                                   size_t i_eq,
                                                   size_t j_unk ) const ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~AP_QIncompHyperelasticity( void ) ;
      AP_QIncompHyperelasticity( AP_QIncompHyperelasticity const& other ) ;
      AP_QIncompHyperelasticity& operator=( 
                                 AP_QIncompHyperelasticity const& other ) ;

      AP_QIncompHyperelasticity( PEL_Object* a_owner, 
                                 PDE_DomainAndFields const* dom,
                                 FE_SetOfParameters const* prms,
                                 PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      AP_QIncompHyperelasticity( void ) ;

      virtual AP_QIncompHyperelasticity* create_replica( 
                                       PEL_Object* a_owner,
                                       PDE_DomainAndFields const* dom,
                                       FE_SetOfParameters const* prms,
                                       PEL_ModuleExplorer* exp ) const ;

   //-- Discrete system building

      void loop_on_cells( FE_TimeIterator const* t_it ) ;

      void loop_on_bounds( FE_TimeIterator const* t_it ) ;

   //-- Volumetric term tools

      double cell_dilatation( void ) ;

      static void add_volumetric_part( AP_KinematicState const* st,
                                       size_t nb_dims,
                                       double pp,
                                       doubleArray2D& Pio2,
                                       doubleArray4D& dPio2dGL ) ;

      double d_gamma( double x ) const ;

      double d2_gamma( double x ) const ;

      double hh( double x ) const ;

      double d_hh( double x ) const ;

      double d2_hh( double x ) const ;

   //-- Class attributes

      static AP_QIncompHyperelasticity const* PROTOTYPE ;

   //-- Attributes

      PDE_DiscreteField* DISP ;
      size_t L_UPDATE ;
      size_t L_EXPLICIT ;

      AP_KinematicState* ST ;
      AP_ConstitutiveLaw* LAW ;

      PDE_SetOfBCs const* BCs ;

      PDE_LocalEquation* ELEMENT_EQ ;
      PDE_LocalEquation* LEQ_2 ;
      PDE_LocalEquation* LEQ_3 ;
      GE_QRprovider const* QRP ;
      PDE_LocalFEcell* cFE ;
      PDE_LocalFEbound* bFE ;

      AP_LoadCalculator* LC ;
      LA_SeqVector* LOAD ;

      FE_Parameter* EXT_LOAD ;

      PDE_SystemNumbering* NMB ;
      LA_Matrix* A ;
      LA_Vector* F ;
      LA_Vector* X ;
      LA_SeqVector* X_LOC ;

      LA_Solver* SOLVER ;
      
      double PENALTY ;
      size_t ITER ;
      size_t MAX_ITER ; //??? NB_ITER_MAX
      double EPS_DISP ; //??? TOL cf CahnHilliard
      double TOL_AL ;

      doubleVector LAMBDA ;
} ;

#endif
