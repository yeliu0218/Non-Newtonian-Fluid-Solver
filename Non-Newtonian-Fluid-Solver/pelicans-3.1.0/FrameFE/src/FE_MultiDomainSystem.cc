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

#include <FE_MultiDomainSystem.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>
#include <PEL_assertions.hh>

#include <LA_Matrix.hh>
#include <LA_Solver.hh>
#include <LA_Vector.hh>

#include <PDE_InterfaceAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_SetOfDomains.hh>
#include <PDE_SystemNumbering.hh>

#include <FE_OneStepIterationOpen.hh>
#include <FE_TimeIterator.hh>

#include <FE_MortarInterfaceDiscretizer.hh>

#include <iostream>
#include <sstream>

using std::cout ;
using std::endl ;
using std::ostringstream ;
using std::string ;

FE_MultiDomainSystem const* 
FE_MultiDomainSystem:: PROTOTYPE = new FE_MultiDomainSystem() ;

struct FE_MultiDomainSystem_ERROR {
   static void n0( void ) ;
   static void n1( std::string const& domain_name ) ;
   static void n2( void ) ;
   static void n3( std::string const& interf_name ) ;
   static void n4( void ) ;
} ;

//-------------------------------------------------------------------------
FE_MultiDomainSystem:: FE_MultiDomainSystem( void )
//-------------------------------------------------------------------------
   : FE_OneStepIteration( "FE_MultiDomainSystem" )
{
}

//-------------------------------------------------------------------------
FE_MultiDomainSystem*
FE_MultiDomainSystem:: create_replica( PEL_Object* a_owner,
                                       PDE_DomainAndFields const* dom,
                                       FE_SetOfParameters const* prms,
                                       PEL_ModuleExplorer* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_MultiDomainSystem:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_MultiDomainSystem_ERROR::n0() ; 
   
   FE_MultiDomainSystem* result = 0 ;
   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_MultiDomainSystem*
FE_MultiDomainSystem:: create_replica( PEL_Object* a_owner,
                                       PDE_SetOfDomains const* sdoms,
                                       FE_SetOfParameters const* prms,
                                       PEL_ModuleExplorer* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_MultiDomainSystem:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, sdoms, prms, exp ) ) ;

   FE_MultiDomainSystem* result = 
      new FE_MultiDomainSystem( a_owner, sdoms, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, sdoms, prms, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_MultiDomainSystem:: FE_MultiDomainSystem( PEL_Object* a_owner,
                                             PDE_SetOfDomains const* sdoms,
                                             FE_SetOfParameters const* prms,
                                             PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, sdoms, exp )
   , SDOMS( sdoms )
   , D_PBS( PEL_Vector::create( this, sdoms->nb_domains() ) )
   , I_PBS( PEL_Vector::create( this, sdoms->nb_interfaces() ) )
   , NMB( 0 )
   , LHS( 0 )
   , RHS( 0 )
   , UNK( 0 )
{
   PEL_LABEL( "FE_MultiDomainSystem:: FE_MultiDomainSystem" ) ;
      
   PEL_ModuleExplorer* e = exp->create_subexplorer( 0, "domain_discretizers" ) ;
   e->start_module_iterator() ;
   for( ; e->is_valid_module() ; e->go_next_module() )
   {
      PEL_ModuleExplorer* ee = e->create_subexplorer( 0 ) ; 
      FE_OneStepIteration* one_it = 
                           FE_OneStepIteration::make( D_PBS, sdoms, prms, ee ) ;
      ee->destroy() ; ee = 0 ;
   }
   e->destroy() ; e = 0 ;
   
   PEL_Vector* dof2Unknowns = PEL_Vector::create( 0, 0 ) ;
   size_t nb_sol = FE_OneStepIterationOpen::nb_objects() ;
   for( size_t i=0 ; i<nb_sol ; ++i )
   {
      FE_OneStepIterationOpen* one_it = FE_OneStepIterationOpen::object( i ) ;
      PDE_DomainAndFields const* dom = one_it->domain() ;
      for( size_t j=0 ; j<SDOMS->nb_domains() ; ++j )
      {
         if( dom == SDOMS->domain( j ) )
         {
            if( D_PBS->at( j ) != 0 ) 
               FE_MultiDomainSystem_ERROR::n1( dom->name() ) ;
            D_PBS->set_at( j, one_it ) ;            
            for( size_t i_unk=0 ; i_unk<one_it->nb_unknowns() ; ++i_unk )
            {
               PDE_LinkDOF2Unknown* lnk = 
                       one_it->link_DOF_2_unknown( i_unk )->create_clone( 0 ) ;
               dof2Unknowns->append( lnk ) ;
            }
         }
      }
   }
   if( D_PBS->count() != D_PBS->index_limit() )
      FE_MultiDomainSystem_ERROR::n2() ;
   
   e = exp->create_subexplorer( 0, "interface_discretizers" ) ;
   e->start_module_iterator() ;
   for( ; e->is_valid_module() ; e->go_next_module() )
   {
      PEL_ModuleExplorer* ee = e->create_subexplorer( 0 ) ; 
      PDE_InterfaceAndFields const* interf = 
                      SDOMS->interface( ee->string_data( "interface" ) ) ;
      FE_MortarInterfaceDiscretizer* ipb = 
         FE_MortarInterfaceDiscretizer::create( I_PBS, interf, D_PBS, ee ) ;
      for( size_t j=0 ; j<SDOMS->nb_interfaces() ; ++j )
      {
         if( interf == SDOMS->interface( j ) )
         {
            if( I_PBS->at( j ) != 0 ) 
               FE_MultiDomainSystem_ERROR::n3( interf->name() ) ;
            I_PBS->set_at( j, ipb ) ;
            for( size_t i_unk=0 ; i_unk<ipb->nb_unknowns() ; ++i_unk )
            {
               PDE_LinkDOF2Unknown* lnk = 
                       ipb->create_link_DOF_2_unknown( 0, i_unk ) ;
               dof2Unknowns->append( lnk ) ;
            }
         }
      }
      ee->destroy() ; ee = 0 ;
   }
   if( I_PBS->count() != I_PBS->index_limit() ) 
      FE_MultiDomainSystem_ERROR::n4() ;
   e->destroy() ; e = 0 ;
   
   NMB = PDE_SystemNumbering::create( this, dof2Unknowns, 
                                      "sequence_of_the_discrete_fields" ) ;
   
   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "LA_Matrix" ) ;
   LHS = LA_Matrix::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
   
   RHS = LHS->create_vector( this ) ;
   UNK = LHS->create_vector( this ) ;
   
   ee = exp->create_subexplorer( 0, "LA_Solver" ) ;
   SOLVER = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   dof2Unknowns->destroy() ; //???????? changer ce nom
}

//-------------------------------------------------------------------------
FE_MultiDomainSystem:: ~FE_MultiDomainSystem( void ) 
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
FE_MultiDomainSystem:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_MultiDomainSystem:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;
   
   start_total_timer( "FE_MultiDomainSystem:: do_one_inner_iteration" ) ;
   
   NMB->reset() ;
   
   size_t dim = NMB->nb_global_unknowns() ;
   LHS->re_initialize( dim, dim ) ;
   RHS->re_initialize( dim ) ;
   UNK->re_initialize( dim ) ;
   
   NMB->define_scatters( RHS ) ;
   
   start_assembling_timer() ;
   for( size_t i=0 ; i<SDOMS->nb_domains() ; ++i )
   {
      FE_OneStepIterationOpen* sd = subdomain_discretizer( i ) ;
      size_t domain_shift = i ;
      sd->assemble_contribution( t_it, LHS, RHS, NMB, domain_shift ) ;
   }

   for( size_t i=0 ; i<SDOMS->nb_interfaces() ; ++i )
   {
      FE_MortarInterfaceDiscretizer* sd = interface_discretizer( i ) ;

      size_t interf_shift = SDOMS->nb_domains() + i ;
      size_t dom_0_shift = SDOMS->index_of_interface_adjacent_domain( i, 0 ) ;
      size_t dom_1_shift = SDOMS->index_of_interface_adjacent_domain( i, 1 ) ;
      
      sd->assemble_contribution( LHS, RHS, NMB,
                                 interf_shift, dom_0_shift, dom_1_shift ) ;
   }
   stop_assembling_timer() ;
   
   start_solving_timer() ;
   
   SOLVER->set_matrix( LHS ) ;
   SOLVER->solve( RHS, UNK ) ;
   SOLVER->unset_matrix() ;
   
   stop_solving_timer() ;
   
   for( size_t i=0 ; i<SDOMS->nb_domains() ; ++i )
   {
      FE_OneStepIterationOpen* sd = subdomain_discretizer( i ) ;
      size_t domain_shift = i ;
      sd->update_DOFs( UNK, NMB, domain_shift ) ;
   }
   
   stop_total_timer() ;
}

//--------------------------------------------------------------------------
FE_OneStepIterationOpen*
FE_MultiDomainSystem:: subdomain_discretizer( size_t i ) const
//--------------------------------------------------------------------------
{
   PEL_CHECK( i < SDOMS->nb_domains() ) ;

   FE_OneStepIterationOpen* result = 
                     static_cast<FE_OneStepIterationOpen*>( D_PBS->at(i) ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
FE_MortarInterfaceDiscretizer*
FE_MultiDomainSystem:: interface_discretizer( size_t i ) const
//--------------------------------------------------------------------------
{
   PEL_CHECK( i < SDOMS->nb_interfaces() ) ;

   FE_MortarInterfaceDiscretizer* result = 
                     static_cast<FE_MortarInterfaceDiscretizer*>( I_PBS->at(i) ) ;
   return( result ) ;
}

//internal-------------------------------------------------------------------
void
FE_MultiDomainSystem_ERROR:: n0( void )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "FE_MultiDomainSystem:" << endl ;
   mesg << "   can only be used if there is more than one domain" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
FE_MultiDomainSystem_ERROR:: n1( std::string const& domain_name )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "FE_MultiDomainSystem:" << endl ;
   mesg << "   there is more than one submodule of" << endl ; 
   mesg << "      MODULE domain_discretizers" << endl ;
   mesg << "   associated with the domain of name" << endl ;
   mesg << "      \"" << domain_name << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
FE_MultiDomainSystem_ERROR:: n2( void )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "FE_MultiDomainSystem:" << endl ;
   mesg << "   it exists at least one domain for which" << endl ;
   mesg << "   there is no associated submodule of" << endl ;
   mesg << "      MODULE domain_discretizers" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
FE_MultiDomainSystem_ERROR:: n3( std::string const& interf_name )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "FE_MultiDomainSystem:" << endl ;
   mesg << "   there is more than one submodule of" << endl ;
   mesg << "      MODULE interface_discretizers" << endl ;
   mesg << "   associated with the interface of name" << endl ;  
   mesg << "      \"" << interf_name << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
FE_MultiDomainSystem_ERROR:: n4( void )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "FE_MultiDomainSystem:" << endl ;
   mesg << "   it exists at least one interface for which" << endl ;
   mesg << "   there is no associated submodule of" << endl ;
   mesg << "      MODULE interface_discretizers" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

