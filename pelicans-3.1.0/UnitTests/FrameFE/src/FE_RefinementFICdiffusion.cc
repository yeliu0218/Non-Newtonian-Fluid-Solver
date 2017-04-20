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

#include <FE_RefinementFICdiffusion.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <boolVector.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>
#include <LA_Vector.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QRprovider.hh>
#include <GE_SetOfPoints.hh>
#include <GE_Vector.hh>

#include <PDE_AdapterCHARMS.hh>
#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_HelperFIC.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_MeshingCoarsener.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfBCs.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <FE.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using std::endl ;
using std::cout ;
using std::string ;
using std::vector ;

FE_RefinementFICdiffusion const* 
FE_RefinementFICdiffusion::PROTOTYPE = new FE_RefinementFICdiffusion() ;

struct FE_RefinementFICdiffusion_ERROR
{
   static void n1( PEL_ModuleExplorer const* exp,
                   std::string const& keyword, 
                   size_t size) ;
   static void n2( void ) ;
} ;

//----------------------------------------------------------------------
FE_RefinementFICdiffusion:: FE_RefinementFICdiffusion( void )
//----------------------------------------------------------------------
   : FE_OneStepIteration( "FE_RefinementFICdiffusion" )
   , HAS_R( 0 )
{
}

//----------------------------------------------------------------------
FE_RefinementFICdiffusion*
FE_RefinementFICdiffusion:: create_replica( PEL_Object* a_owner,
                                            PDE_DomainAndFields const* dom,
                                            FE_SetOfParameters const* prms,
                                            PEL_ModuleExplorer* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementFICdiffusion:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_RefinementFICdiffusion* result = 
                new FE_RefinementFICdiffusion( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_RefinementFICdiffusion:: FE_RefinementFICdiffusion( 
                                         PEL_Object* a_owner,
                                         PDE_DomainAndFields const* dom,
                                         FE_SetOfParameters const* prms,
                                         PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , UU( dom->set_of_discrete_fields()->item( 
                                    exp->string_data( "unknown_field" ) ) )
   , L_UPDATE( exp->int_data( "level_to_update" ) )
   , ALPHA( exp->double_data( "coef_unsteady" ) )
   , KAPPA( exp->double_data( "coef_diffusion" ) )
   , PI( prms->item( exp->string_data( "param_source" ) ) )
   , BCs( dom->set_of_boundary_conditions() )
   , DA( dom->adapter_CHARMS() )
   , RESID( LA_SeqVector::create( this, 0 ) )
   , HAS_R( 0 )
   , sFE( dom->create_CursorFEside( this ) )
   , bFE( dom->create_LocalFEbound( this ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , VERTS( dom->set_of_vertices() )
   , PT( GE_Point::create( this, dom->nb_space_dimensions() ) )
   , QRP_PI( GE_QRprovider::object( 
             exp->string_data( "quadrature_rule_provider_for_source" ) ) )
   , X_LOC( LA_SeqVector::create( this, 0 ) )
   , d_U0( LA_SeqVector::create( this, 0 ) )
   , TOL( exp->double_data( "tolerance" ) )
{
   PEL_LABEL( "FE_RefinementFICdiffusion:: FE_RefinementFICdiffusion" ) ;

   sFE->require_field_calculation( UU, PDE_LocalFE::node ) ;
   bFE->require_field_calculation( UU, PDE_LocalFE::node ) ;
   cFE->require_field_calculation( UU, PDE_LocalFE::node ) ;
   
   PI->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;
   
   PDE_LinkDOF2Unknown* uu_link = PDE_LinkDOF2Unknown::create( 0, UU, 
                                          "sequence_of_the_components", 
                                          true ) ;
   NMB = PDE_SystemNumbering::create( this, uu_link ) ;  
   
   PEL_ModuleExplorer const* ee = exp->create_subexplorer( 0, "LA_Matrix" ) ;
   A = LA_Matrix::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;

   F = A->create_vector( this ) ;
   X = A->create_vector( this ) ;

   ee = exp->create_subexplorer( 0, "LA_Solver" ) ;
   SOLVER = LA_Solver::make( this, ee ) ;
   ee->destroy() ; ee = 0 ;
}

//----------------------------------------------------------------------
FE_RefinementFICdiffusion:: ~FE_RefinementFICdiffusion( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
FE_RefinementFICdiffusion:: do_one_inner_iteration( 
                                         FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementFICdiffusion:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;
   
   start_total_timer( "FE_RefinementFICdiffusion:: do_one_inner_iteration" ) ;
   // --------------
   
   PDE_MeshingCoarsener* coar = DA->meshing_coarsener() ;
   coar->prepare_for_coarsening() ;
   size_t nb_levels = coar->nb_levels() ;
   size_t level_max = nb_levels-1 ;
   OBSERVED_NODES.resize( nb_levels, boolVector( 0 ) ) ;
   RESID->re_initialize( cFE->nb_meshes() ) ;
   HAS_R.re_initialize( cFE->nb_meshes() ) ;
   HAS_R.set( true ) ;
   
   increase_indent() ;

   size_t level = level_max ;
   if( verbose_level() >= 2 )
      PEL::out() << indent() << "level: " << level << endl ;
   prepare_for_multigrid( level, level_max ) ;
   for( ; level != 0 ; --level )
   {
      coar->do_one_coarsening() ;
      size_t current_level = level-1 ;
      if( verbose_level() >= 2 )
         cout << indent() << "level: " << current_level << endl ;
      prepare_for_multigrid( current_level, level_max ) ;
   }
   
   if( verbose_level() >= 2 ) 
      PEL::out() << indent() << "initialization" << endl ;
   increase_indent() ;
   if( verbose_level() >= 2 )
      PEL::out() << indent() << "level: " << level << endl ;
   build_and_solve( t_it, level ) ;
   d_U0->re_initialize( NMB->link()->unknown_vector_size() ) ;
   LA_Scatter const* sca = NMB->scatter() ;
   sca->get( X, d_U0 ) ;
   decrease_indent() ;   
   
   size_t i_cycle = 0 ;
   do
   {
      ++i_cycle ;
      if( verbose_level() >= 2 ) 
         PEL::out() << indent() << "cycle: " << i_cycle << endl ;
      increase_indent() ;
      PEL_ASSERT( level == 0 ) ;
      if( verbose_level() >= 2 )
      {
         PEL::out() << indent() << "level: " << level << endl ;
         PEL::out() << indent() 
                    << "   coarse contrib to the flux residual" << endl ;
      }
      loop_on_sides( t_it, level, 1 ) ;
      for( ; level != level_max ; ++level )
      {
         coar->do_one_uncoarsening() ;
         size_t current_level = level+1 ;
         if( verbose_level() >= 2 )
            PEL::out() << indent() << "level: " << current_level << endl ;
         build_and_solve( t_it, current_level ) ;
         if( current_level != level_max )
         {
            if( verbose_level() >= 2 )
               PEL::out() << indent() 
                          << "   coarse contrib to the flux residual" 
                          << endl ;
            loop_on_sides( t_it, current_level, 1 ) ;
            //??? contributions des bounds
         }
      }
      for( ; level != 0 ; --level )
      {
         if( verbose_level() >= 2 )
            PEL::out() << indent() 
                       << "   fine contrib to the flux residual of level " 
                       << level-1 << endl ;
         loop_on_sides( t_it, level, 2 ) ;
         //??? contributions des bounds
         coar->do_one_coarsening() ;
         size_t current_level = level-1 ;
         if( verbose_level() >= 2 )
            PEL::out() << indent() << "level: " << current_level << endl ;
         build_and_solve( t_it, current_level ) ;
      }
      decrease_indent() ;
   } while( !converged() ) ;

   // return to the composite discretization
   for( ; level != level_max ; ++level )
   {
      coar->do_one_uncoarsening() ;
      size_t current_level = level+1 ;
      if( verbose_level() >= 2 )
         PEL::out() << indent() << "level: " << current_level << endl ;
   }
   
   decrease_indent() ;
   
   stop_total_timer() ;
   // -------------
}

//-----------------------------------------------------------------------
void
FE_RefinementFICdiffusion:: prepare_for_multigrid( size_t level, 
                                                   size_t level_max )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementFICdiffusion:: prepare_for_multigrid" ) ;
   
   boolVector& observed_nodes = OBSERVED_NODES[level] ;
   observed_nodes.re_initialize( UU->nb_nodes() ) ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      if( level==level_max || 
          ( (level != 0) && (cFE->refinement_level() == level-1) ) )
      {
         RESID->set_item( cFE->mesh_id(), PEL::bad_double() ) ;
         HAS_R( cFE->mesh_id() ) = false ;
      }
      if( cFE->refinement_level() != level ) continue ;
      //                                     --------
      
      if( cFE->nb_local_nodes( UU ) != 1 ) 
                                    FE_RefinementFICdiffusion_ERROR::n2() ;

      observed_nodes( cFE->global_node( UU, 0 ) ) = true ;
   } 
}

//-----------------------------------------------------------------------
void
FE_RefinementFICdiffusion:: build_and_solve( FE_TimeIterator const* t_it ,
                                             size_t level )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementFICdiffusion:: build_and_solve" ) ;
   
   NMB->reset( OBSERVED_NODES[ level ] ) ;
   
   size_t n_glob = NMB->nb_global_unknowns() ;
   
   A->re_initialize( n_glob, n_glob ) ;
   F->re_initialize( n_glob ) ;
   X->re_initialize( n_glob ) ;
   
   X_LOC->re_initialize( NMB->link()->unknown_vector_size() ) ;
   
   NMB->define_scatters( X ) ;
   
   start_assembling_timer() ;
   // -------------------

   loop_on_cells( t_it, level ) ;
   loop_on_sides( t_it, level, 0 ) ;
   loop_on_bounds( t_it, level ) ;
   
   stop_assembling_timer() ;

   start_solving_timer() ;
   // ----------------

   A->synchronize() ;
   F->synchronize() ;
   SOLVER->set_matrix( A ) ;
   SOLVER->solve( F, X ) ;
   SOLVER->unset_matrix() ;
   
   stop_solving_timer() ;
   // ---------------
   
   if( ! SOLVER->solution_is_achieved() )
   {
      PEL_Error::object()->display_info(
         "*** FE_RefinementFICdiffusion error\n"
         "    convergence failure when solving the discrete problem" ) ;
      notify_inner_iterations_stage_failure() ;
   }
   else
   {
      if( verbose_level() >= 2 )
      {
         PEL::out() << indent() << "   update of " << UU->name()
                    << "(" << L_UPDATE << ")" << endl ;
      }
      LA_Scatter const* sca = NMB->scatter() ;
      sca->get( X, X_LOC ) ;
      UU->update_free_DOFs_value( L_UPDATE, X_LOC, NMB->link() ) ;
   }  
}

//------------------------------------------------------------------------
void
FE_RefinementFICdiffusion:: loop_on_sides( FE_TimeIterator const* t_it,
                                           size_t level,
                                           int mode )
//------------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementFICdiffusion:: loop_on_sides" ) ;
   
   PDE_LocalFEcell const* fe_K = sFE->adjacent_localFEcell( 0 ) ;
   PDE_LocalFEcell const* fe_L = sFE->adjacent_localFEcell( 1 ) ;
   
   doubleVector values( 3 ) ;
   for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
   {
      if( sFE->refinement_level() != level ) continue ;
      //                                     --------
      
      size_t n_K = fe_K->global_node( UU, 0 ) ;
      size_t n_L = fe_L->global_node( UU, 0 ) ;
      
      double area = FE::side_measure( sFE ) ;

      size_t u_K = PEL::bad_index() ;
      if( NMB->link()->DOF_is_unknown( n_K, 0 ) )
      {
         PEL_ASSERT( fe_K->refinement_level() == level ) ;
         u_K = NMB->global_unknown_for_DOF( n_K, 0 ) ;
      }
      size_t u_L = PEL::bad_index() ;
      if( NMB->link()->DOF_is_unknown( n_L, 0 ) )
      {
         PEL_ASSERT( fe_L->refinement_level() == level ) ;
         u_L = NMB->global_unknown_for_DOF( n_L, 0 ) ;
      }

      if( u_K != PEL::bad_index() && u_L != PEL::bad_index() )
      {
         double d_KL = sFE->distance_to_adjacent_finite_volume_center( 0 ) +
                       sFE->distance_to_adjacent_finite_volume_center( 1 ) ;
         double xx_0 =   KAPPA * area / d_KL ;
         double xx_1 = - KAPPA * area / d_KL ;
         
         if( mode==1 )
         {
            // xx = - (flux grossier)
            double xx = - xx_0 * ( UU->DOF_value( L_UPDATE, n_K, 0 ) -
                                   UU->DOF_value( L_UPDATE, n_L, 0 ) ) ;
            size_t id_K = fe_K->mesh_id() ;
            if( HAS_R( id_K ) ) RESID->add_to_item( id_K,  xx ) ;
            size_t id_L = fe_L->mesh_id() ;
            if( HAS_R( id_L ) ) RESID->add_to_item( id_L, -xx ) ;
         }
         else if( mode==2 )
         {            
            PDE_HelperFIC* hf = sFE->helper_FIC() ;

            // xx = +( flux fin )
            double xx = xx_0 * ( UU->DOF_value( L_UPDATE, n_K, 0 ) -
                                 UU->DOF_value( L_UPDATE, n_L, 0 ) ) ;
            size_t p_id_K = hf->parent_cell_id( 0 ) ;
            size_t p_id_L = hf->parent_cell_id( 1 ) ;
            if( p_id_K != p_id_L )
            {
               RESID->add_to_item( p_id_K,  xx ) ;
               RESID->add_to_item( p_id_L, -xx ) ;
            }
         }
         else
         {
            PEL_ASSERT( mode == 0 ) ;
            
            A->add_to_item( u_K, u_K,  xx_0 ) ;
            A->add_to_item( u_L, u_K, -xx_0 ) ;

            A->add_to_item( u_K, u_L,  xx_1 ) ;
            A->add_to_item( u_L, u_L, -xx_1 ) ;
         }
      }
      else if( u_K == PEL::bad_index() )
      {
         //??? il semble que l'on ne passe jamais ici
         PEL_ASSERT( false ) ;
      }
      else
      {
         PEL_ASSERT( u_L == PEL::bad_index() ) ;
         PEL_ASSERT( u_K != PEL::bad_index() ) ;
         
         PDE_HelperFIC* hf = sFE->helper_FIC() ;
         
         GE_Point const* cell_pt = fe_K->polyhedron()->finite_volume_center() ;
         PEL_ASSERT( cell_pt != 0 ) ;
         GE_Point const* side_pt = sFE->polyhedron()->center() ;
         PT->set_coordinate( 0, 
               2.0*side_pt->coordinate( 0 ) - cell_pt->coordinate( 0 ) ) ;
         PT->set_coordinate( 1, 
               2.0*side_pt->coordinate( 1 ) - cell_pt->coordinate( 1 ) ) ;
         double dd = PT->distance( cell_pt ) ;
         double d_sigma = sFE->distance_to_adjacent_finite_volume_center( 0 ) ;
         PEL_ASSERT( PEL::abs( dd - 2.0*d_sigma ) < 1.e-8 ) ;
         double xx_0 =   KAPPA * area / dd ;
         double xx_1 = - KAPPA * area / dd ;
         hf->prepare_for_interpolation() ;
         for( size_t i=0 ; i!=3 ; ++i )
         {
            values( i ) = UU->DOF_value( L_UPDATE, hf->node( i, UU ), 0 ) ;
         }
         hf->set_calculation_point( PT ) ;
         double value_u_L = hf->interpolated_value( values ) ;
         if( mode == 1 )
         {
            size_t id_K = fe_K->mesh_id() ;
            if( HAS_R( id_K ) )
            {
               // xx = -( flux grossier )
               double xx = - xx_0 * ( UU->DOF_value( L_UPDATE, n_K, 0 ) -
                                      value_u_L ) ;
               RESID->add_to_item( id_K,  xx ) ;
            }
         }
         else if( mode == 2 )
         {
            // xx = +( flux fin )
            double xx = xx_0 * ( UU->DOF_value( L_UPDATE, n_K, 0 ) -
                                 value_u_L ) ;
            size_t p_id_K = hf->parent_cell_id( 0 ) ;
            RESID->add_to_item( p_id_K,  xx ) ;
         }
         else
         {
            PEL_ASSERT( mode == 0 ) ;
            F->add_to_item( u_K, -xx_1*value_u_L ) ; 
            A->add_to_item( u_K, u_K, xx_0 ) ;
         }
      }
   }
}

//------------------------------------------------------------------------
void
FE_RefinementFICdiffusion:: loop_on_bounds( FE_TimeIterator const* t_it, 
                                            size_t level )
//------------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementFICdiffusion:: loop_on_bounds" ) ;

   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      if( bFE->refinement_level() != level ) continue ;
      //                                     --------
      
      GE_Color const* color = bFE->color() ;

      PEL_ModuleExplorer const* ee = BCs->BC_explorer( color, UU ) ;

      size_t n_K = bFE->global_node( UU, 0 ) ;
      size_t u_K = NMB->global_unknown_for_DOF( n_K, 0 ) ;

      double area = FE::bound_measure( bFE ) ;

      string const& type = ee->string_data( "type" ) ;
      if( type == "Dirichlet" )
      {
         doubleVector const& uimp = ee->doubleVector_data( "imposed_value" ) ;
         if( uimp.size() != 1 )
            FE_RefinementFICdiffusion_ERROR::n1( ee, "imposed_value", 1 ) ;

         double h = bFE->distance_to_adjacent_finite_volume_center() ;

         double xx_l = KAPPA * area / h ;
         double xx_r = xx_l * uimp( 0 ) ;

         A->add_to_item( u_K, u_K, xx_l ) ;
         F->add_to_item( u_K, xx_r ) ; 
      }
      else if( type == "Neumann" )
      {
         doubleVector const& fimp = ee->doubleVector_data( "flux_value" ) ;
         if( fimp.size() != 1 )
            FE_RefinementFICdiffusion_ERROR::n1( ee, "flux_value", 1 ) ;
         F->add_to_item( u_K, area*fimp( 0 ) ) ;
      }
      else
         raise_bad_BC_type( type, "\"Dirichlet\"\n\"Neumann\"", UU ) ;
   }
}

//------------------------------------------------------------------------
void
FE_RefinementFICdiffusion:: loop_on_cells( FE_TimeIterator const* t_it,
                                           size_t level )
//------------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementFICdiffusion:: loop_on_cells" ) ;
   
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      if( cFE->refinement_level() != level ) continue ;
      //                                     --------
      
      size_t n_K = cFE->global_node( UU, 0 ) ;
      size_t u_K = NMB->global_unknown_for_DOF( n_K, 0 ) ;

      double int_source = 0.0 ;
      cFE->start_IP_iterator( QRP_PI ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         int_source += cFE->weight_of_IP() * 
                       PI->cell_value_at_IP( t_it, cFE ) ;
      }
      double xx_r = int_source ;
      
      size_t id_K = cFE->mesh_id() ;
      if( HAS_R( id_K ) ) xx_r += RESID->item( id_K ) ;

      double vol = FE::cell_measure( cFE ) ;
      double xx_l = vol * ALPHA ;

      A->add_to_item( u_K, u_K,  xx_l ) ;
      F->add_to_item( u_K, xx_r ) ;
   }
}

//------------------------------------------------------------------------
bool
FE_RefinementFICdiffusion:: converged( void ) const
//------------------------------------------------------------------------
{
   PEL_LABEL( "FE_RefinementFICdiffusion:: converged" ) ;
   
   LA_Scatter const* sca = NMB->scatter() ;
   sca->get( X, X_LOC ) ;
   d_U0->sum( X_LOC, -1.0 ) ;
   double xx = d_U0->max_norm() ;
   
   bool result = ( xx < TOL ) ;
   if( verbose_level() >= 2 )
      PEL::out() << indent() << "   max|dU| = " << xx << endl ;
   
   d_U0->set( X_LOC ) ;

   return( result ) ;
}

//internal--------------------------------------------------------------
void
FE_RefinementFICdiffusion_ERROR:: n1( PEL_ModuleExplorer const* exp,
                                      std::string const& keyword,
                                      size_t size )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << std::endl << "*** FE_RefinementFICdiffusion:" 
       << std::endl << std::endl ;
   msg << "In module: " << exp->absolute_path_name() << endl ;
   msg << "the data of keyword: " << keyword << endl ;
   msg << "should be a DoubleVector of size: " << size ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void
FE_RefinementFICdiffusion_ERROR:: n2( void )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << std::endl << "*** FE_RefinementFICdiffusion:" << std::endl ;
   msg << "      there should be only one active node on each cell" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

