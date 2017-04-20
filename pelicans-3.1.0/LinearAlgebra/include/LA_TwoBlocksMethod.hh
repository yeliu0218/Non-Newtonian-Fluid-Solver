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

#ifndef LA_TWO_BLOCKS_METHOD_HH
#define LA_TWO_BLOCKS_METHOD_HH

#include <PEL_Object.hh>

#include <LA.hh>

#include <string>

class PEL_ModuleExplorer ;
class PEL_ObjectRegister ;

class LA_Implementation ;
class LA_Matrix ;
class LA_Vector ;

/*
Estimators for the solution two block linear systems of the form

       |  A   tB  |  | U |   | F |
       |          |  |   | = |   |
       |  B    C  |  | P |   | G |

PUBLISHED
*/

class PEL_EXPORT LA_TwoBlocksMethod : public PEL_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization
      
      static LA_TwoBlocksMethod* make( PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp ) ;
      
   //-- System profile
      
      void set_matrix_prototype( LA_Matrix const* mat ) ;
      
      LA::DistributionStrategy distribution_strategy( void ) const ;
      
      LA_Implementation const* implementation( void ) const ;
      
      void re_initialize_internals( size_t nv_glob, 
                                    size_t np_glob,
                                    size_t nv_loc, 
                                    size_t np_loc ) ;
      
      size_t nb_local_rows_U( void ) const ;
      
      size_t nb_local_rows_P( void ) const ;
      
   //-- Auxiliary items
      
      virtual bool dtinv_is_required( void ) const ;
      
      void set_dtinv( double value ) ;
         
      virtual bool S_is_required( void ) const ;
      
      virtual void set_S( LA_Vector* a_S ) ;
      
      virtual bool L_is_required( void ) const ;
      
      virtual void set_L( LA_Matrix* a_L ) ;
      
      virtual bool K_is_required( void ) const ;
      
      virtual void set_K( LA_Vector* a_K ) ;
      
      virtual bool MV_is_required( void ) const ;
      
      virtual void set_MV( LA_Matrix* a_MV ) ;
      
   //-- System setting
      
      void set_system( LA_Matrix* a_A, LA_Matrix* a_B,
                       LA_Vector* a_F, LA_Vector* a_G,
                       LA_Matrix* a_C = 0 ) ;
      
      bool system_is_set( void ) const ;
      
      virtual void unset_system( void ) ;
      
   //-- Estimation
      
      void estimate_unknowns( bool has_init_U, LA_Vector* U, 
                              bool has_init_P, LA_Vector* P ) ;
      
      bool successful_estimation( void ) const ;

   //-- Input - Output

      void set_indent( std::string const& an_indent ) ;

      virtual void print_times( size_t indent_width ) const ;

   protected: //--------------------------------------------------------

   //-- Plug in
      
      virtual ~LA_TwoBlocksMethod( void ) ;

      LA_TwoBlocksMethod( std::string const& a_name ) ;

      LA_TwoBlocksMethod( PEL_Object* a_owner,
                          PEL_ModuleExplorer const* exp ) ;
                  
      virtual LA_TwoBlocksMethod* create_replica(
                                PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) const = 0 ;

      bool is_a_prototype( void ) const ;
      
   //-- System profile
      
      virtual void set_matrix_prototype_sub( LA_Matrix const* mat ) = 0 ;

      virtual void re_initialize_internals_sub( size_t nv_glob, 
                                                size_t np_glob,
                                                size_t nv_loc, 
                                                size_t np_loc,
                                                size_t& nv_loc_final, 
                                                size_t& np_loc_final ) = 0 ;

   //-- Auxiliary items
      
      double dtinv( void ) const ;
      
   //-- System setting
      
      virtual void set_system_sub( LA_Matrix* a_A, LA_Matrix* a_B,
                                   LA_Vector* a_F, LA_Vector* a_G,
                                   LA_Matrix* a_C ) = 0 ;
      
      virtual void unset_system_sub( void ) = 0 ;
      
      void check_zero_C( LA_Matrix const* a_C ) const ;
      
      void raise_invalid_usage( void ) ;   
               
   //-- Estimation
         
      virtual void estimate_unknowns_sub( bool has_init_U, LA_Vector* U, 
                                          bool has_init_P, LA_Vector* P ) = 0 ;
      
      void notify_success( bool is_successful ) ;
               
   //-- Input - Output
      
      size_t verbose_level( void ) const ;
      
      std::string const& indent( void ) const ;

   //-- Preconditions, Postconditions, Invariant
      
      bool set_S_PRE( LA_Vector* a_S ) const ;
      
      bool set_L_PRE( LA_Matrix* a_L ) const ;
      
      bool set_K_PRE( LA_Vector* a_K ) const ;
      
      bool set_MV_PRE( LA_Matrix* a_MV ) const ;
      
      bool create_replica_PRE( PEL_Object* a_owner, 
                               PEL_ModuleExplorer const* exp ) const ;
      
      bool create_replica_POST( LA_TwoBlocksMethod const* result,
                                PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) const ;

   private: //----------------------------------------------------------

      LA_TwoBlocksMethod( void ) ;
      LA_TwoBlocksMethod( LA_TwoBlocksMethod const& other ) ;
      LA_TwoBlocksMethod& operator=( LA_TwoBlocksMethod const& other ) ;

      static PEL_ObjectRegister* plugins_map( void ) ;
      
   //-- Attributes
      
      bool const IS_PROTO ;
      LA::DistributionStrategy DIST_STRAT ;
      LA_Implementation const* IMPL ;
      
      size_t NB_LOC_U ;
      size_t NB_LOC_P ;
      
      bool IS_SET ;
      bool SUCCESS ;
      
      size_t VERBOSE ;
      std::string INDENT ;
      
      double BoverDT ;
} ; 

#endif
