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

#ifndef PDE_ALGEBRAIC_COARSENER_HH
#define PDE_ALGEBRAIC_COARSENER_HH

#include <PEL_Object.hh>

#include <boolVector.hh>

#include <PDE_AdapterCHARMS.hh>

#include <map>
#include <vector>

class LA_Matrix ;
class LA_Vector ;

class PDE_BasisFunctionCell ;
class PDE_DiscreteField ;
class PDE_LinkDOF2Unknown ;
class PDE_GridFE ;
class PDE_SystemNumbering ;

/*
HIGHLY UNSTABLE CLASS
*/

class PEL_EXPORT PDE_AlgebraicCoarsener : public PEL_Object
{
   public: //-----------------------------------------------------------

      void reset( void ) ;

      void prepare_for_coarsening( PDE_SystemNumbering const* nmb ) ;

      size_t nb_levels( void ) const ;

      bool coarsening_is_possible( void ) const ;

      void do_one_coarsening( void ) ;

      void build_current_prolongation_matrix( LA_Matrix* mat ) const ;

      void build_smoothing_lines( LA_Vector* smoothing_lines ) const ;

      void build_current_level_unknowns(
                                 LA_Vector* current_level_unknowns ) const ;

      size_t current_fine_level( void ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_AlgebraicCoarsener( void ) ;
     ~PDE_AlgebraicCoarsener( void ) ;
      PDE_AlgebraicCoarsener( PDE_AlgebraicCoarsener const& other ) ;
      PDE_AlgebraicCoarsener& operator=(
                              PDE_AlgebraicCoarsener const& other ) ;

      friend PDE_AlgebraicCoarsener*
             PDE_AdapterCHARMS:: algebraic_coarsener( void ) const ;

      PDE_AlgebraicCoarsener( PEL_Object* a_owner, PDE_GridFE* a_grid ) ;

   //-- Attributes

      PDE_GridFE* GRID ;

      size_t NB_BLOCKS ;
      std::vector< PDE_DiscreteField const* > FF ;
      std::vector< std::vector< PDE_BasisFunctionCell*> > BFS ;
      std::vector< boolVector > NN_FINE ;
      std::vector< boolVector > NN_COAR ;
      size_t FINE_LEVEL ;
      size_t LEVEL_MAX ;

      PDE_SystemNumbering* ROW_NMB ;
      PDE_SystemNumbering* COL_NMB ;

} ;

#endif
