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

#ifndef PDE_GEOMETRIC_MULTILEVEL_PC_HH
#define PDE_GEOMETRIC_MULTILEVEL_PC_HH

#include <LA_Preconditioner.hh>

#include <map>

class PEL_Vector ;

class PEL_DistributedPartition ;

class PDE_AlgebraicCoarsener ;
class PDE_DomainAndFields ;
class PDE_LinkDOF2Unknown ;
class PDE_SystemNumbering ;
class size_t_vector ;
class intVector ;

/*
HIGHLY UNSTABLE CLASS

FRAMEWORK INSTANTIATION
PUBLISHED
*/

class PEL_EXPORT PDE_GeometricMultilevel_PC : public LA_Preconditioner
{
   public: //----------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_GeometricMultilevel_PC* object( std::string const& a_name ) ;

   //-- Instance characteristics

      std::string const& name( void ) const ;

   //-- Status

      virtual bool is_valid( void ) const ;

      virtual size_t dimension( void ) const ;

   //-- Building

      void set_discretization_scene( PDE_DomainAndFields const* dom,
                                     PDE_SystemNumbering const* nmb ) ;

      virtual void build( LA_Matrix const* mat ) ;

      virtual void unbuild( void ) ;

   //-- Linear system solution

      virtual size_t nb_cycles_performed( void ) const = 0 ;

   protected: //-------------------------------------------------------

   //-- Derivation

      virtual ~PDE_GeometricMultilevel_PC( void ) ;

      PDE_GeometricMultilevel_PC( std::string const& class_name ) ;

      PDE_GeometricMultilevel_PC( PEL_Object* a_owner,
                                  PEL_ModuleExplorer const* exp ) ;

      PDE_GeometricMultilevel_PC( PEL_Object* a_owner,
                                  PDE_GeometricMultilevel_PC const* other ) ;

   //-- Discretization scene

      PDE_SystemNumbering const* system_numbering( void ) const ;

      size_t nb_fine_unknowns( void ) const ;

      size_t nb_local_fine_unknowns( void ) const ;

      PDE_AlgebraicCoarsener* algebraic_coarsener( void ) const ;

   //-- Internal data

      size_t nb_levels( void ) const ;

      LA_Matrix* mat_of_level( size_t level ) const ;

      LA_Matrix const* finest_mat( void ) const ;

      LA_Vector* rhs_of_level( size_t level ) const ;

      LA_Vector* res_of_level( size_t level ) const ;

      LA_Vector* sol_of_level( size_t level ) const ;

      LA_Matrix const* coarse_to_fine( size_t coarse_level ) const ;

      LA_Vector const* smoothing_lines_of_level( size_t level ) const ;

      LA_Vector const* unknowns_of_level( size_t level ) const ;

   //-- Smoothing

      static void smooth_GaussSeidel( size_t nb_steps,
                                      LA_Matrix const* mat,
                                      LA_Vector const* rhs,
                                      LA_Vector* sol ) ;

      static void smooth_GaussSeidel( size_t nb_steps,
                                      LA_Vector const* to_be_smoothed,
                                      LA_Matrix const* mat,
                                      LA_Vector const* rhs,
                                      LA_Vector* sol ) ;

      static void  priority_rows( PEL_DistributedPartition const* row_dist,
                                  intVector const& rows_to_send,
                                  size_t_vector& priority_rows ) ;

      static void  extra_columns( LA_Matrix const* mat,
                                  intVector& ext_cols ) ;


      static void compute_residual( LA_Matrix const* mat,
                                    LA_Vector const* rhs,
                                    LA_Vector const* sol,
                                    LA_Vector* residual ) ;

      static void print_residuals( std::string const& indent, size_t n,
                            double norm_res, double norm_res_0 ) ;

   private: //---------------------------------------------------------

      PDE_GeometricMultilevel_PC( void ) ;
      PDE_GeometricMultilevel_PC( PDE_GeometricMultilevel_PC const& other ) ;
      PDE_GeometricMultilevel_PC& operator=(
                                  PDE_GeometricMultilevel_PC const& other ) ;

   //-- Internals

      void get_prolongation_matrices( PDE_AlgebraicCoarsener* coar ) ;

   //-- Class attributes

      static std::map< std::string, PDE_GeometricMultilevel_PC* > OBJS ;

   //-- Attributes

      std::string NAME ;
      PDE_SystemNumbering const* NMB ;
      PDE_AlgebraicCoarsener* COAR ;

      LA_Matrix const* MAT_PROTO ;
      PEL_Vector* AA ;               // vector of LA_Matrix*
      LA_Matrix const* FINEST_A ;
      PEL_Vector* B ;                // vector of LA_Vector*
      PEL_Vector* X ;                // vector of LA_Vector*
      PEL_Vector* RES ;              // vector of LA_Vector*

      PEL_Vector* COARSE_TO_FINE ;   // vector of LA_Matrix*
      PEL_Vector* SMOO_LINE ;        // vector of LA_Vector*
      PEL_Vector* UNKNOWNS_LEVEL ;
      size_t NB_LEVELS ;

      bool BUILD_OK ;
} ;

#endif


