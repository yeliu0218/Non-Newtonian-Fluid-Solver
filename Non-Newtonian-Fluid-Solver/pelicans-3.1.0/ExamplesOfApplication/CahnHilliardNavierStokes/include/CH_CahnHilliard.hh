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

#ifndef CH_CAHN_HILLIARD_HH
#define CH_CAHN_HILLIARD_HH

#include <FE_OneStepIterationOpen.hh>

#include <doubleVector.hh>

class LA_Matrix ;
class LA_SeqVector ;
class LA_Solver ;
class LA_Vector ;

class GE_QRprovider ;

class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LinkDOF2Unknown ;
class PDE_LocalFEcell ;
class PDE_LocalEquation ;
class PDE_ResultSaver ;
class PDE_SystemNumbering ;

class FE_Parameter ;

class CH_BulkChemicalPotential ;

/*
PUBLISHED
*/

class CH_CahnHilliard : public FE_OneStepIterationOpen
{
   public: //-----------------------------------------------------------

   //-- Jacobian

      virtual size_t nb_unknowns( void ) const  ;

      virtual PDE_DiscreteField* field( size_t i_unk ) const ;

      virtual size_t level_of_field( size_t i_unk ) const ;

      virtual PDE_LinkDOF2Unknown const* link_DOF_2_unknown( 
                                              size_t i_unk ) const ;

      virtual void build_function_and_jacobian( FE_TimeIterator const* t_it ) ;

      virtual LA_SeqVector const* create_function( PEL_Object* a_owner,
                                                   size_t i_unk ) const ;

      virtual LA_SeqMatrix const* create_jacobian( PEL_Object* a_owner,
                                                   size_t i_eq,
                                                   size_t j_unk ) const ;

   //-- Substeps of the step by step progression

      virtual void do_one_inner_iteration( FE_TimeIterator const* t_it ) ;

   //-- Savings for post-processing

      virtual void save_other_than_time_and_fields(
                                               FE_TimeIterator const* t_it,
                                               PDE_ResultSaver* rs ) ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
                
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~CH_CahnHilliard( void ) ;
      CH_CahnHilliard( CH_CahnHilliard const& other ) ;
      CH_CahnHilliard& operator=( CH_CahnHilliard const& other ) ;

      CH_CahnHilliard( PEL_Object* a_owner,
                       PDE_DomainAndFields const* dom,
                       FE_SetOfParameters const* prms,
                       PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      CH_CahnHilliard( void ) ;

      virtual CH_CahnHilliard* create_replica(
                                      PEL_Object* a_owner,
                                      PDE_DomainAndFields const* dom,
                                      FE_SetOfParameters const* prms,
			                             PEL_ModuleExplorer* exp ) const ;

   //-- Discrete system building

      bool convergence_achieved( void ) const ;

      void loop_on_cells( FE_TimeIterator const* t_it ) ;

      void loop_Mi_Ci( PDE_DiscreteField const* ci,
                       PDE_DiscreteField const* ci_exp,
                       PDE_DiscreteField const* mi,
                       double& vol,
                       double dt ) ;

      void loop_Mi_Mi( PDE_DiscreteField const* mi,
                       double sigma ) ;
      
      void loop_Mi_Mi_MobMod( PDE_DiscreteField const* mi,
                              PDE_DiscreteField const* mi_exp,
                              double sigma ) ;

      void loop_Ci_Mi( PDE_DiscreteField const* mi ) ;

      void loop_Ci_Ci( size_t i, double cap, double theta ) ;

      void loop_Ci_Cj( size_t i, size_t j ) ;

    //-- Class attributes

      static CH_CahnHilliard const* PROTOTYPE ;

   //-- Attributes

      PDE_DiscreteField* C1 ;
      PDE_DiscreteField* C1_EXP ;
      PDE_DiscreteField* M1 ;
      PDE_DiscreteField* M1_EXP ;
      PDE_DiscreteField* C2 ;
      PDE_DiscreteField* C2_EXP ;
      PDE_DiscreteField* M2 ;
      PDE_DiscreteField* M2_EXP ;

      size_t L_UPDATE ;
      size_t L_EXPLICIT ;

      PDE_DiscreteField* AA ;
      size_t L_AA ;

      double EPS ;
      double CAP1 ;
      double CAP2 ;
      double THETA ;
      size_t NB_PRE_EULER_STEP ;
      
      double MOB_CST ;
      double MOB_DEG ;
      bool MOB_EXP ;
      bool MOB_MOD ;
      double MOB_MAX ;
      
      CH_BulkChemicalPotential* BULK_MU ;

      PDE_LocalEquation* ELEMENT_EQ ;
      GE_QRprovider const* QRP ;
      PDE_LocalFEcell* cFE ;

      double VOL1 ;
      double VOL2 ;

      size_t ITER ;
      size_t NB_ITER_MAX ;
      double TOL ;

      size_t idx_C1 ;
      size_t idx_M1 ;
      size_t idx_C2 ;
      size_t idx_M2 ;

      PDE_SystemNumbering* NMB ;

      LA_Matrix* LHS_cst ;
      LA_Matrix* LHS ;
      LA_Vector* RHS ;
      LA_Vector* UNK ;
      LA_SeqVector* UNK_LOC ;

      LA_Solver* SOLVER ;

      doubleVector PRINTED_RES ;
} ;

#endif

