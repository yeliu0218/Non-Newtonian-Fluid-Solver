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

#ifndef PDE_ADAPTER_CHARMS_TEST
#define PDE_ADAPTER_CHARMS_TEST

#include <PEL_ObjectTest.hh>

#include <iosfwd>
#include <stringVector.hh>

class PEL_ContextSimple ;
class PEL_DataWithContext ;
class PEL_DoubleVector ;

class GE_Color ;
class GE_Point ;
class GE_QRprovider ;
class GE_SetOfPoints ;

class PDE_CursorFEside ;
class PDE_DiscreteField ;
class PDE_DomainAndFields ;
class PDE_LocalFE ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;
class PDE_MeshingCoarsener ;
class PDE_SetOfDiscreteFields ;

class PEL_EXPORT PDE_AdapterCHARMS_TEST : public PEL_ObjectTest
{
   public: //---------------------------------------------------------------

   protected: //------------------------------------------------------------

   private: //--------------------------------------------------------------

      PDE_AdapterCHARMS_TEST( void ) ;
     ~PDE_AdapterCHARMS_TEST( void ) ;
      PDE_AdapterCHARMS_TEST( PDE_AdapterCHARMS_TEST const& other ) ;
      PDE_AdapterCHARMS_TEST& operator=( 
                              PDE_AdapterCHARMS_TEST const& other ) ;

      virtual void process_one_test( PEL_ModuleExplorer const* exp ) ;

      void check_meshing_coarsening( PDE_DomainAndFields const* dom,
                                     PDE_MeshingCoarsener* coar ) ;
      
      void print_discretization( PDE_DomainAndFields const* dom,
                                 std::string const& fname,
                                 size_t level ) ;
      
      void comp_discretizations( std::string const& fname_1,
                                    std::string const& fname_2,
                                    size_t level,
                                    bool& ok ) ;
      
      void iterate_1( size_t ref_step,
                      PDE_SetOfDiscreteFields* sdf,
                      PDE_LocalFE* fe,
                      std::ostream& os ) ;

      void process_domain( size_t ref_step,
                           PDE_DomainAndFields const* dom,
                           std::ostream& os ) ;

      void iterate( size_t ref_step,
                    PDE_SetOfDiscreteFields* sdf,
                    PDE_LocalFE* fe,
                    std::ostream& os ) ;

      void check_node_location( size_t ref_step,
                                PDE_SetOfDiscreteFields* sdf,
                                PDE_LocalFE* fe,
                                bool& ok_node_loc ) const ;

      void check_interpolation( size_t ref_step, PDE_LocalFE* fe ) ;
      
      void check_imposed_values( size_t ref_step, PDE_LocalFEbound* bfe ) ;

      void check_values_at_vertices( size_t ref_step,
                                     PDE_DiscreteField* ff,
                                     GE_SetOfPoints const* vertices,
                                     PDE_LocalFE* fe ) ;
      
      void check_cell_and_faces( size_t ref_step, 
                                 PDE_DomainAndFields const* dom ) ;

      void check_covering_of_cells_boundary( 
                                 size_t ref_step, 
                                 PDE_DomainAndFields const* dom  ) ;
      
      void fill_context_with_coordinates( GE_Point const* pt ) ;

      void display_error( PDE_LocalFE const* fe,
                          double theo, double xx ) const ;

   //-- Class attributes

      static PDE_AdapterCHARMS_TEST const* REGISTRATOR ;

   //-- Attributes

      std::string TEST_NAME ;

      GE_QRprovider const* QRP ;
      PEL_DataWithContext* VAL_CHECK ;
      PDE_DiscreteField const* UU_CHECK ;
      double EPS_INTERP ;
      double MIN_INTERP ;
      stringVector FFS_VERTS ;
      double EPS_VERTS ;
      double MIN_VERTS ;
      PEL_ContextSimple* CONTEXT ;
      PEL_DoubleVector* COORDS ;
      bool CHECK_BD_CELLS ;
      PEL_DataWithContext* VAL_BDIMP ;
      PDE_DiscreteField const* UU_BDIMP ;
      GE_Color const* COLOR_BDIMP ;
      double EPS_BDIMP ;
      double MIN_BDIMP ;
} ;

#endif
