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

#ifndef PDE_CROSS_PROCESS_UNKNOWN_NUMBERING_HH
#define PDE_CROSS_PROCESS_UNKNOWN_NUMBERING_HH

#include <PEL_Object.hh>

#include <size_t_array2D.hh>
#include <size_t_vector.hh>
#include <intVector.hh>

class PEL_Communicator ;
class PDE_CrossProcessNodeNumbering ;
class PDE_LinkDOF2Unknown ;

/*
PUBLISHED
*/

class PEL_EXPORT PDE_CrossProcessUnknownNumbering : public PEL_Object
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_CrossProcessUnknownNumbering* create(
                                         PEL_Object* a_owner,
                                         PDE_LinkDOF2Unknown const* a_link ) ;

      void reset( void ) ;

   //-- Characteristics

      PDE_LinkDOF2Unknown const* link( void ) const ;

      PEL_Communicator const* communicator( void ) const ;

   //-- On-process to cross-process mapping

      // number of unknowns handled by the current process
      size_t nb_unknowns_of_current_process( void ) const ;

      // number of unknowns handled by the process `rank'
      size_t nb_unknowns_on_process( size_t rank ) const ;

      // number of global (cross-process) unknowns
      size_t nb_global_unknowns( void ) const ;

      // global (cross-process) unknown index of the unknown
      // whose local (on-process) index is `i'
      size_t global_unknown_index( size_t i ) const ;

      // rank of process owning associated discrete node
      size_t rank_of_process_handling( size_t i ) const ;

      // global (cross-process) index of the unknown associated
      // to the DOF defined by `n' and component `ic' (or `PEL::bad_index()'
      // if there is no unknown for that DOF)
      size_t global_unknown_linked_to_DOF( size_t n,
                                           size_t ic ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PDE_CrossProcessUnknownNumbering( void ) ;
     ~PDE_CrossProcessUnknownNumbering( void ) ;
      PDE_CrossProcessUnknownNumbering(
                            PDE_CrossProcessUnknownNumbering const& other ) ;
      PDE_CrossProcessUnknownNumbering& operator=(
                            PDE_CrossProcessUnknownNumbering const& other ) ;

      PDE_CrossProcessUnknownNumbering( PEL_Object* a_owner,
                                        PDE_LinkDOF2Unknown const* a_link ) ;

      void globalize( size_t_array2D& unknown ) ;
      void raise_globalize_error( std::string const& mes,
                                  size_t_array2D const& unknown ) const ;

  //-- Internals

      enum OrderingType{ sequenceOfTheComponents, sequenceOfTheNodes } ;

      void set_ordering_option( std::string const& option ) ;

  //-- Attributes

      PDE_LinkDOF2Unknown const* LINK ;
      OrderingType ORDERING ;
      PDE_CrossProcessNodeNumbering const* GLOB_NODE ;

      PEL_Communicator const* COMM ;

      intVector NB_UNKNOWNS_ON_PROC ;
      size_t NB_UNKNOWNS ;
      size_t NB_HANDLED_UNK ;
      size_t_vector LOC_2_GLOB ;
      size_t_vector OWNER ;
} ;

#endif
