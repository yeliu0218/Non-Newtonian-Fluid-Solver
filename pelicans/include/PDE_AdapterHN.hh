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

#ifndef PDE_ADAPTER_HN_HH
#define PDE_ADAPTER_HN_HH

#include <PEL_Object.hh>

#include <PDE_DomainAndFields.hh>

class PEL_Context ;
class PEL_DataWithContext ;
class PEL_DoubleVector ;
class PEL_ModuleExplorer ;

class GE_Point ;
class GE_ReferencePolyhedronRefiner ;
class GE_SetOfPoints ;

class PDE_CellFE ;
class PDE_RefinementPatternProvider ;
class PDE_SetOfDiscreteFields ;
class PDE_SetOfBasisFunctions ;
class PDE_DomainAndFields ;
class PDE_GridFE ;

/*
HIGHLY UNSTABLE CLASS
*/

class PEL_EXPORT PDE_AdapterHN : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Adaptation

      void reset( void ) ;

      void adapt( void ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_AdapterHN( void ) ;
     ~PDE_AdapterHN( void ) ;
      PDE_AdapterHN( PDE_AdapterHN const& other ) ;
      PDE_AdapterHN& operator=( PDE_AdapterHN const& other ) ;

      friend PDE_AdapterHN* PDE_DomainAndFields:: adapter_HN( void ) const ;

      PDE_AdapterHN( PEL_Object* a_owner,
                     PDE_DomainAndFields const* dom,
                     PDE_GridFE* a_grid,
                     PDE_SetOfBasisFunctions* bf_set,
                     PEL_ModuleExplorer const* exp ) ;
      
    //-- Internals
         
      void refine_cell( PDE_CellFE* ccell ) ;
      
      void build_new_nodes_and_bfs( PDE_CellFE const* ccell, size_t ee ) ;

      void find_cells_to_refine( PEL_Vector* cells ) ;

   //-- Attributes

      PDE_DomainAndFields const* DOM ;
      PDE_GridFE* GRID ;
      PDE_SetOfDiscreteFields const* DFs ;
      PDE_RefinementPatternProvider const* RPP ;
      size_t VERB_LEVEL ;
      
      PDE_CrossProcessBFNumbering* BF_NUMBERING ;
      PDE_SetOfBasisFunctions* BF_SET ;
      
      PEL_Context* CTX ;
      mutable PEL_DoubleVector* COORDS ;
      PEL_DataWithContext* R_INDIC ;
} ;

#endif
