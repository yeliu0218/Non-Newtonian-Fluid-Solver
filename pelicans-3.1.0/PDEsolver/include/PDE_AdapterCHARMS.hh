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

#ifndef PDE_ADAPTER_CHARMS_HH
#define PDE_ADAPTER_CHARMS_HH

#include <PEL_Object.hh>

#include <PDE_DomainAndFields.hh>
#include <boolVector.hh>
#include <size_t_array2D.hh>
#include <size_t_vector.hh>

#include <list>
#include <map>
#include <set>
#include <vector>

class PEL_List ;
class PEL_ListIterator ;
class PEL_ModuleExplorer ;
class size_t_vector ;

class GE_Point ;
class GE_ReferencePolyhedronRefiner ;
class GE_SetOfPoints ;

class PDE_Activator ;
class PDE_AdaptationIndicator ;
class PDE_AdaptationRequest ;
class PDE_AlgebraicCoarsener ;
class PDE_BasisFunctionCell ;
class PDE_CellFE ;
class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_GridFE ;
class PDE_MeshingCoarsener ;
class PDE_RefinementPatternProvider ;
class PDE_SetOfBasisFunctions ;
class PDE_SetOfDiscreteFields ;
class PDE_FaceFE ;

/*
HIGHLY UNSTABLE CLASS
*/

class PEL_EXPORT PDE_AdapterCHARMS : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Configuration

      void append_indicator( PDE_AdaptationIndicator* indic ) ;

      void add_excluded_field( PDE_DiscreteField const* ff ) ;

      bool is_excluded( PDE_DiscreteField const* ff ) const ;

   //-- Adaptation

      void reset( void ) ;

      void adapt( void ) ;

      void unsplit_meshes( void ) ;

      bool something_changed( void ) const ;

   //-- Coarsening

      PDE_AlgebraicCoarsener* algebraic_coarsener( void ) const ;

      PDE_MeshingCoarsener* meshing_coarsener( void ) const ;

   //-- Input - Output

      virtual void print_statistics( std::ostream& os,
                                     size_t indent_width ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_AdapterCHARMS( void ) ;
     ~PDE_AdapterCHARMS( void ) ;
      PDE_AdapterCHARMS( PDE_AdapterCHARMS const& other ) ;
      PDE_AdapterCHARMS& operator=( PDE_AdapterCHARMS const& other ) ;

      friend PDE_AdapterCHARMS*
             PDE_DomainAndFields:: adapter_CHARMS( void ) const ;

      PDE_AdapterCHARMS( PEL_Object* a_owner,
                         PDE_DomainAndFields const* dom,
                         PDE_GridFE* a_grid,
                         PDE_SetOfBasisFunctions* bf_set,
                         PEL_ModuleExplorer const* exp ) ;

      void possibly_unsplit_cell( PDE_CellFE* ccell,
                                  std::list< PDE_CellFE* >& cccs ) ;

      void split_meshes( PDE_AdaptationRequest* adap ) ;

      void build_new_nodes_and_bfs( PDE_AdaptationRequest* adap ) ;

      void reset_counters( void ) ;

    //-- Attributes

      bool HB ;
      size_t MAX_REF_LEVEL ;
      PDE_Activator* ACTIVATOR ;
      PDE_DomainAndFields const* DOM ;
      PDE_GridFE* GRID ;
      PDE_SetOfDiscreteFields const* DFs ;
      GE_Point* PT ;
      PDE_RefinementPatternProvider const* RPP ;
      PEL_List* INDICS ; // list of PDE_AdaptationIndicator*
      PEL_ListIterator* INDICS_IT ;
      PEL_ModuleExplorer const* EE ;
      bool SOMETHING_CHANGED ;
      size_t VERB_LEVEL ;

      size_t NB_ACTIVE_CELLS ;
      size_t NB_ACT_CELLS ;
      size_t NB_DEACT_CELLS ;
      
      mutable PDE_AlgebraicCoarsener* A_COARSENER ;
      mutable PDE_MeshingCoarsener*   M_COARSENER ;

      PDE_CrossProcessBFNumbering* BF_NUMBERING ;
      PDE_SetOfBasisFunctions* BF_SET ;
} ;

#endif
