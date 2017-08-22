#include <FE_InitAfterAdaptation.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <GE_QRprovider.hh>

#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>

#include <PDE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <FE.hh>
#include <FE_TimeIterator.hh>

#include <iomanip>
#include <iostream>
#include <sstream>

using std::string ;

FE_InitAfterAdaptation const* 
FE_InitAfterAdaptation::PROTOTYPE = new FE_InitAfterAdaptation() ;

//---------------------------------------------------------------------------
FE_InitAfterAdaptation:: FE_InitAfterAdaptation( void )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( "FE_InitAfterAdaptation" )
{
}

//---------------------------------------------------------------------------
FE_InitAfterAdaptation*
FE_InitAfterAdaptation:: create_replica( PEL_Object* a_owner,
                                         PDE_DomainAndFields const* dom,
                                         FE_SetOfParameters const* prms,
                                         PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_InitAfterAdaptation:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_InitAfterAdaptation* result = 
                      new FE_InitAfterAdaptation( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
FE_InitAfterAdaptation:: FE_InitAfterAdaptation( 
                                     PEL_Object* a_owner,
                                     PDE_DomainAndFields const* dom,
                                     FE_SetOfParameters const* prms,
                                     PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , QRP( 0 )
   , cFE( dom->create_LocalFEcell( this ) )
   , A( 0 )
   , U_LOC( LA_SeqVector::create( this, 0 ) )
{
   PEL_LABEL( "FE_InitAfterAdaptation:: FE_InitAfterAdaptation" ) ;

   bool has_L2_proj = false ;

   PDE_SetOfDiscreteFields const* dfs = dom->set_of_discrete_fields() ;
   PEL_ModuleExplorer* se = 
                  exp->create_subexplorer( 0, "list_of_PDE_DiscreteField" ) ;
   se->start_module_iterator() ;
   for( ; se->is_valid_module() ; se->go_next_module() )
   {
      PEL_ModuleExplorer const* sse = se->create_subexplorer( 0 ) ;
      PDE_DiscreteField* ff = dfs->item( sse->string_data( "current" ) ) ;
      FIELDs.push_back( ff ) ;
      L_FIELDs.push_back( sse->int_data( "level_of_current" ) ) ;
      std::string f_type =  sse->string_data( "type" )  ;
      if( f_type == "L2_projection_of_explicit" )
      {
         TYPEs.push_back( f_type ) ;
         PDE_DiscreteField const* ff_exp = 
                                  dfs->item( sse->string_data( "explicit" ) ) ;
         FIELD_EXPs.push_back( ff_exp ) ;
         L_FIELD_EXPs.push_back( sse->int_data( "level_of_explicit" ) ) ;
         cFE->require_field_calculation( ff_exp, PDE_LocalFE::N ) ;
         has_L2_proj = true ;
      }
      else if( f_type == "nullify_new_DOFs" )
      {  
         TYPEs.push_back( f_type ) ;
         FIELD_EXPs.push_back( 0 ) ;
         L_FIELD_EXPs.push_back( 0 ) ;
      }
      else
      {
         PEL_Error::object()->raise_bad_data_value( exp, "type", 
				    "    \"L2_projection_of_explicit\"\n"
				    "    \"nullify_new_DOFs\"\n" ) ;
      }
      sse->destroy() ;
      cFE->require_field_calculation( ff, PDE_LocalFE::N ) ;
   }
   se->destroy() ;

   if( exp->has_module( "L2_projection" ) )
   {
      se = exp->create_subexplorer( 0, "L2_projection" ) ;
      QRP = GE_QRprovider::object( se->string_data( "QRprovider_name" ) ) ;
      
      for( size_t i=0 ; i<FIELDs.size() ; ++i )
      {
         PDE_LinkDOF2Unknown* ll = PDE_LinkDOF2Unknown::create( 0, FIELDs[i],
                                       "sequence_of_the_components",
                                       true ) ;
         NMBs.push_back( PDE_SystemNumbering::create( this, ll ) ) ;
      }
      
      PEL_ModuleExplorer const* ee = se->create_subexplorer( 0, "LA_Matrix" ) ;
      A = LA_Matrix::make( this, ee ) ;
      ee->destroy() ; ee = 0 ;

      F = A->create_vector( this ) ;
      U = A->create_vector( this ) ;

      ee = se->create_subexplorer( 0, "LA_Solver" ) ;
      SOLVER = LA_Solver::make( this, ee ) ;
      ee->destroy() ; ee = 0 ;
      
      se->destroy() ; se = 0 ;
   }

   if( has_L2_proj && A==0 )
   {
      PEL_Error::object()->raise_plain( "missing MODULE L2_projection" ) ;
   }
}

//---------------------------------------------------------------------------
FE_InitAfterAdaptation:: ~FE_InitAfterAdaptation( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void 
FE_InitAfterAdaptation:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_InitAfterAdaptation:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "FE_InitAfterAdaptation:: do_one_inner_iteration" ) ;
   
   for( size_t i=0 ; i<FIELDs.size() ; ++i )
   {
      if( TYPEs[i] == "L2_projection_of_explicit" )
      {
         NMBs[i]->reset() ;
         
         size_t n_glob = NMBs[i]->nb_global_unknowns() ;
         size_t n_loc  = NMBs[i]->nb_unknowns_on_current_process() ;

         A->re_initialize( n_glob, n_glob, n_loc, n_loc ) ;
         F->re_initialize( n_glob, n_loc ) ;
         U->re_initialize( n_glob, n_loc ) ;
         
         U_LOC->re_initialize( NMBs[i]->link()->unknown_vector_size() ) ;
         
         NMBs[i]->define_scatters( U ) ;
         
         start_assembling_timer() ;
         // ---------------------
         
         loop_on_cells( t_it, i ) ;
         
         stop_assembling_timer() ;
         start_solving_timer() ;
         // ------------------
         
         A->synchronize() ;
         F->synchronize() ;
         SOLVER->set_matrix( A ) ;
         SOLVER->solve( F, U ) ; 
         SOLVER->unset_matrix() ;
         
         stop_solving_timer() ;
   
         if( ! SOLVER->solution_is_achieved() )
         {
            PEL_Error::object()->display_info(
               "*** FE_InitAfterAdaptation error\n"
               "    No convergence of the linear solver" ) ;
            notify_inner_iterations_stage_failure() ;
         }
         else
         {
            LA_Scatter const* sca = NMBs[i]->scatter() ;
            sca->get( U, U_LOC ) ;
            FIELDs[i]->update_free_DOFs_value( L_FIELDs[i], 
                                               U_LOC, 
                                               NMBs[i]->link() ) ;
         }
      }
      else if( TYPEs[i] == "nullify_new_DOFs" )
      {
         for( size_t n=0 ; n<FIELDs[i]->nb_nodes() ; ++n )
         {
            if( FIELDs[i]->DOF_value( L_FIELDs[i], n ) == PEL::bad_double() )
            {
               FIELDs[i]->set_DOF_value( 0, n, 0.0 ) ;
            }
         }
      }
   }
   stop_total_timer() ;
}

//------------------------------------------------------------------------
void
FE_InitAfterAdaptation:: loop_on_cells( FE_TimeIterator const* t_it,
                                        size_t i )
//------------------------------------------------------------------------
{
   PEL_LABEL( "FE_InitAfterAdaptation:: loop_on_cells" ) ;

   size_t nbc = FIELDs[i]->nb_components() ;
   
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->set_row_and_col_fields( FIELDs[i], FIELDs[i] ) ;
      ELEMENT_EQ->initialize( cFE->row_field_node_connectivity(), nbc,
                              cFE->col_field_node_connectivity(), nbc ) ;

      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         doubleVector ff_exp( nbc ) ;
         for( size_t d=0 ; d<nbc ; ++d )
         {
            ff_exp( d )  = cFE->value_at_IP( FIELD_EXPs[i], L_FIELD_EXPs[i], d )  ;
         }
         FE::add_row_col_S( ELEMENT_EQ, cFE, 1.0 ) ;

         FE::add_row( ELEMENT_EQ, cFE, ff_exp ) ;
      }

      PDE::assemble_in_matrix_vector_0( A, F, ELEMENT_EQ, NMBs[i] ) ;
   }
}

