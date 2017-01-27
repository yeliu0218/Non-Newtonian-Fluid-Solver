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

#ifndef PDE_SYSTEM_NUMBERING_HH
#define PDE_SYSTEM_NUMBERING_HH

#include <PEL_Object.hh>
#include <size_t_vector.hh>
#include <size_t_array2D.hh>

#include <vector>

class PEL_Vector ;
class size_t_vector ;

class size_t_array2D ;
class intVector ;
class boolVector ;

class LA_Scatter ;
class LA_Vector ;

class PDE_LinkDOF2Unknown ;
class PDE_LocalEquation ;

class PEL_EXPORT PDE_SystemNumbering : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_SystemNumbering* create( PEL_Object* a_owner,
                                          PEL_Vector const* dof2unks,
                                          std::string const& ordering,
                                          size_t a_verbose_level = 0 ) ;

      static PDE_SystemNumbering* create( PEL_Object* a_owner,
                                          PDE_LinkDOF2Unknown* dof2unk,
                                          size_t a_verbose_level = 0 ) ;

      virtual PDE_SystemNumbering* create_clone( PEL_Object* a_owner ) const ;
      
      void reset( void ) ;

      void reset( std::vector<boolVector> const& observed_nodes ) ;
      
      void reset( boolVector const& observed_nodes ) ;
      
   //-- Characteristics

      bool is_distributed( void ) const ;

   //-- Blocks of the global algebraic system

      // number of attached `PDE_LinkDOF2Unknown::' objects
      size_t nb_links( void ) const ;

      // `i_link'-th attached `PDE_LinkDOF2Unknown::' object
      PDE_LinkDOF2Unknown const* link( size_t i_link ) const ;

      // only attached `PDE_LinkDOF2Unknown::' object
      PDE_LinkDOF2Unknown const* link( void ) const ;

   //-- Items of the global algebraic system

      size_t nb_global_unknowns( void ) const ;

      size_t nb_unknowns_on_current_process( void ) const ;//handled by

      size_t nb_unknowns_on_process( size_t rank ) const ;

      size_t global_unknown_for_DOF( size_t n,
                                     size_t ic,
                                     size_t i_link ) const ;

      size_t global_unknown_for_DOF( size_t n,
                                     size_t ic ) const ;

   //-- Scatters

      void define_scatters( LA_Vector const* vec ) ;

      bool scatters_are_defined( void ) const ;

      LA_Scatter const* scatter( size_t i_link ) const ;

      LA_Scatter const* scatter( void ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_SystemNumbering( void ) ;
     ~PDE_SystemNumbering( void ) ;
      PDE_SystemNumbering( PDE_SystemNumbering const& other ) ;
      PDE_SystemNumbering& operator=( PDE_SystemNumbering const& other ) ;

      PDE_SystemNumbering( PEL_Object* a_owner,
                           PEL_Vector const* dof2unks,
                           std::string const& ordering,
                           size_t a_verbose_level ) ;

      PDE_SystemNumbering( PEL_Object* a_owner,
                           PDE_SystemNumbering const* other ) ;

      PDE_SystemNumbering( PEL_Object* a_owner,
                           PDE_LinkDOF2Unknown* dof2unk,
                           size_t a_verbose_level ) ;

   //-- Internals

      size_t_vector& idx_local( size_t i ) const ;

      size_t_vector& idx_global( size_t i ) const ;

      void reset_sizes( size_t verbose_level ) ;
      void reset_numbering( void ) ;

      size_t global_unknown_index( size_t rank,
                                   size_t i_link,
                                   size_t idx_unk ) const ;

      void set_ordering_option( std::string const& option ) ;

      enum BlockOrderingType{ SequenceOfDiscreteFields, SequenceOfUnknowns  } ;

   //-- Attributes

      size_t const NB_LINKS ;
      PEL_Vector* const LINKS ; // PDE_LinkDOF2Unknown*

      BlockOrderingType ORDERING ;

      size_t_array2D SHIFT ;

      size_t_vector** IDX_LOCS ;
      size_t_vector** IDX_GLOBS ;

      bool OK_SCATTERS ;
      PEL_Vector* SCATTERS ; // LA_Scatter*

      bool DISTRIBUTED ;
      size_t_vector NB_HANDLED_UNK ;
      size_t MAT_SIZE ;

      size_t VERB ;
} ;

#endif
