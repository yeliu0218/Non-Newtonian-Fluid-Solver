#include <CH_AveragesSaver.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <doubleArray2D.hh>
#include <size_t_vector.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QRprovider.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_ResultSaver.hh>

#include <FE.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>

#include <CH_BulkChemicalPotential.hh>

#include <fstream>
#include <ios>
#include <iostream>
#include <iomanip>
#include <sstream>

using std::endl ;
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::ostringstream ;
using std::string ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

struct CH_AveragesSaver_ERROR
{
   static void n0( std::string const& fname ) ;
   static void n1( std::string const& module_name, 
                   size_t invalid_dim,
                   size_t_vector valid_dims ) ;
} ;

CH_AveragesSaver const* 
CH_AveragesSaver::PROTOTYPE = new CH_AveragesSaver() ;

//---------------------------------------------------------------------------
CH_AveragesSaver:: CH_AveragesSaver( void )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( "CH_AveragesSaver" )
   , CC_MIN( 0 )
   , SAVING_TIMES( doubleVector( 0 ) )
{
}

//---------------------------------------------------------------------------
CH_AveragesSaver*
CH_AveragesSaver:: create_replica( PEL_Object* a_owner,
                                   PDE_DomainAndFields const* dom,
                                   FE_SetOfParameters const* prms,
                                   PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_AveragesSaver:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   CH_AveragesSaver* result = new CH_AveragesSaver( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
CH_AveragesSaver:: CH_AveragesSaver( PEL_Object* a_owner,
                                     PDE_DomainAndFields const* dom,
                                     FE_SetOfParameters const* prms,
                                     PEL_ModuleExplorer* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , FNAME( exp->string_data( "output_file" ) )
   , cFE( dom->create_LocalFEcell( this ) ) 
   , SOMETHING_TO_DO( false )
   , TOTAL_FREE_ENERGY( false )
   , THICKNESS( PEL::bad_double() )
   , BULK_MU( 0 )
   , C1( 0 )
   , L_C1( PEL::bad_index() )
   , C2( 0 )
   , L_C2( PEL::bad_index() )
   , QRP_TFE( 0 )
   , KINETIC_ENERGY( false )
   , UU( 0 )
   , L_UU( PEL::bad_index() )
   , DENS( 0 )
   , QRP_KE( 0 )
   , CENTER( 0 )
   , CC( 0 )
   , L_CC( PEL::bad_index() )
   , VV( 0 )
   , L_VV( PEL::bad_index() )
   , CC_MIN( 0 )
   , QRP_CV( 0 )
   , PERIMETER( false )
   , CP( 0 )
   , L_CP( PEL::bad_index() )
   , QRP_PER( 0 )
   , NEXT_SAVING_TIME( PEL::bad_double() )
   , SAVING_TIMES( doubleVector( 0 ) )
{
   PEL_LABEL( "CH_AveragesSaver:: CH_AveragesSaver" ) ;
   
   if( exp->has_entry( "saving_times" ) )
   {
      SAVING_TIMES = exp->doubleVector_data( "saving_times" ) ;  
   }
   
   PEL_ModuleExplorer* se = 0 ;
   if( exp->has_module( "total_free_energy" ) )
   {
      se = exp->create_subexplorer( 0, "total_free_energy" ) ;
      TOTAL_FREE_ENERGY = true ;
      read_for_total_free_energy( dom, prms, se ) ;  
      se->destroy() ; se = 0 ;
      SOMETHING_TO_DO = true ;
   }
   if( exp->has_module( "kinetic_energy" ) ) 
   {
      se = exp->create_subexplorer( 0, "kinetic_energy" ) ;
      KINETIC_ENERGY = true ;
      read_for_kinetic_energy( dom, prms, se ) ;
      se->destroy() ; se = 0 ;
      SOMETHING_TO_DO = true ;
   }
   if( exp->has_module( "center_volume" ) )
   {
      se = exp->create_subexplorer( 0, "center_volume" ) ;
      CENTER = true ;
      read_for_center_volume( dom, prms, se ) ;
      se->destroy() ; se = 0 ;
      SOMETHING_TO_DO = true ;
   } 
   if( exp->has_module( "perimeter" ) )
   {
      se = exp->create_subexplorer( 0, "perimeter" ) ;
      if( dom->nb_space_dimensions() != 2 )
      {
         size_t_vector valid_dims( 1 ) ;
         valid_dims( 0 ) = 2 ;
         CH_AveragesSaver_ERROR:: n1( se->name(),
                                      dom->nb_space_dimensions(),
                                      valid_dims ) ;
      }
      PERIMETER = true ;
      read_for_perimeter( dom, prms, se ) ;
      se->destroy() ; se = 0 ;
      SOMETHING_TO_DO = true ;
   }
}

//---------------------------------------------------------------------------
CH_AveragesSaver:: ~CH_AveragesSaver( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
CH_AveragesSaver:: do_before_time_stepping( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_AveragesSaver:: do_before_time_stepping" ) ;
   PEL_CHECK_PRE( do_before_time_stepping_PRE( t_it ) ) ;
   
   if( SOMETHING_TO_DO && communicator()->rank() == 0 )
   {
      std::ofstream ofs( FNAME.c_str(), std::ios::out | std::ios::trunc ) ;
      if( !ofs ) CH_AveragesSaver_ERROR::n0( FNAME ) ;
      ofs << "#" << endl ;
      ofs << "# CH_AveragesSaver generated file" << endl ;
      ofs << "# meaning of the columns:" << endl ;
      ofs << "#    TIME : current time" << endl ;
      if( TOTAL_FREE_ENERGY )
         ofs << "#    TFE  : total free energy" << endl ;
      if( KINETIC_ENERGY )
         ofs << "#    KE   : kinetic energy" << endl ;
      if( CENTER )
      {
         for( size_t i=0 ; i<CC_MIN.size() ; ++i )
         {
            ofs << "#    VOL(" << i << ") : volume " 
                << "computed from \"" << CC->name() << "\" " 
                << "and threshold "<< CC_MIN( i ) << endl;
         }
         if ( FE::geometry() == FE::axisymmetrical )
         {
            for( size_t i=0 ; i<CC_MIN.size() ; ++i )
            {
               ofs << "#    X1(" << i << ") : 1-coordinate of center "
               << "computed from \"" << CC->name() << "\" " 
               << "and threshold "<< CC_MIN( i ) << endl;
            }         
            for( size_t i=0 ; i<CC_MIN.size() ; ++i )
            {
               ofs << "#    V1(" << i << ") : 1-coordinate of center velocity "
               << "computed from \"" << VV->name() << "\" " 
               << "and threshold "<< CC_MIN( i ) << endl;
            }         
         }
         else
         {
            for( size_t d=0 ; d<cFE->nb_space_dimensions() ; ++d )
            {
               for( size_t i=0 ; i<CC_MIN.size() ; ++i )
               {
                  ofs << "#    X" << d << "(" << i << ") : " << d 
                      << "-coordinate of center "
                      << "computed from \"" << CC->name() << "\" " 
                      << "and threshold "<< CC_MIN( i ) << endl;
               }
            }
            for( size_t d=0 ; d<cFE->nb_space_dimensions() ; ++d )
            {
               for( size_t i=0 ; i<CC_MIN.size() ; ++i )
               {
                  ofs << "#    V" << d << "(" << i << ") : " << d 
                      << "-coordinate of center velocity"
                      << "computed from \"" << VV->name() << "\" " 
                      << "and threshold "<< CC_MIN( i ) << endl;
               }
            }
         }
      }
      if( PERIMETER )
      {
         ofs << "#    PER  : perimeter" << endl  ;
      }
      ofs << "#" << endl ;
      ofs << "#       TIME" ;
      if( TOTAL_FREE_ENERGY )
         ofs << "          TFE" ;
      if( KINETIC_ENERGY )
         ofs << "           KE" ;
      if( CENTER )
      {
         for( size_t i=0 ; i<CC_MIN.size() ; ++i )
         {
            ofs << "       VOL(" << i << ")" ;
         }
         if ( FE::geometry() == FE::axisymmetrical )
         {
            for( size_t i=0 ; i<CC_MIN.size() ; ++i )
            {
               ofs << "        X1(" << i << ")" ;
            }         
            for( size_t i=0 ; i<CC_MIN.size() ; ++i )
            {
               ofs << "        V1(" << i << ")" ;
            }         
         }
         else
         {
            for( size_t d=0 ; d<cFE->nb_space_dimensions() ; ++d )
            {
               for( size_t i=0 ; i<CC_MIN.size() ; ++i )
               {
                  ofs << "        X" << d << "(" << i << ")" ;
               }
            }
            for( size_t d=0 ; d < cFE->nb_space_dimensions() ; ++d )
            {
               for( size_t i=0 ; i < CC_MIN.size() ; ++i )
               {
                  ofs << "        V" << d << "(" << i << ")" ;
               }
            }         
         }         
      }
      if( PERIMETER )
      {
         ofs << "          PER" ;
      }
      ofs << endl ;
      ofs.close() ;
   }
   
   if( !t_it->table_of_times_is_valid( SAVING_TIMES ) )
   {
      t_it->raise_invalid_table_of_times( "CH_AveragesSaver", "saving_times", 
                                          SAVING_TIMES ) ;
   }
   if( SAVING_TIMES.has( t_it->time() ) )
   {
      save_all( t_it ) ;
   }
   NEXT_SAVING_TIME = t_it->next_time_in_table( SAVING_TIMES ) ;
}

//---------------------------------------------------------------------------
void
CH_AveragesSaver:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_AveragesSaver:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;
}

//---------------------------------------------------------------------------
void
CH_AveragesSaver:: do_after_time_adaptation( FE_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_AveragesSaver:: do_after_time_adaptation" ) ;
   PEL_CHECK_PRE( do_after_time_adaptation_PRE( t_it ) ) ;

   if( SOMETHING_TO_DO && ( SAVING_TIMES.size() != 0 ) )
   {
      if( FE_TimeIterator::greater_or_equal( t_it->time(),
                                             NEXT_SAVING_TIME ) )
      {
         start_total_timer( "CH_AveragesSaver:: do_one_inner_iteration" ) ;
         
         save_all( t_it ) ;
         
         NEXT_SAVING_TIME = t_it->next_time_in_table( SAVING_TIMES ) ;
         stop_total_timer() ;      
      }
   }
}

//---------------------------------------------------------------------------
void
CH_AveragesSaver:: save_other_than_time_and_fields( 
                                                 FE_TimeIterator const* t_it, 
                                                 PDE_ResultSaver* rs )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_AveragesSaver:: save_other_than_time_and_fields" ) ;
   PEL_CHECK_PRE( save_other_than_time_and_fields_PRE( t_it, rs ) ) ;

   if( SOMETHING_TO_DO && ( SAVING_TIMES.size() == 0 ) )
   {
      start_total_timer( "CH_AveragesSaver:: save_other_than_time_and_fields" ) ;
      
      save_all( t_it ) ;
      
      stop_total_timer() ;
   }      
}

//----------------------------------------------------------------------------
void
CH_AveragesSaver:: read_for_total_free_energy( PDE_DomainAndFields const* dom,
                                               FE_SetOfParameters const* prms,
                                               PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "CH_AveragesSaver:: read_for_total_free_energy" ) ;
   
   PDE_SetOfDiscreteFields const* dfs = dom->set_of_discrete_fields() ;

   THICKNESS = exp->double_data( "thickness" ) ;
   QRP_TFE = GE_QRprovider::object( 
                            exp->string_data( "quadrature_rule_provider" ) ) ;
   
   C1 = dfs->item( exp->string_data( "phase_field_1" ) ) ;
   L_C1 = exp->int_data( "level_1" ) ;
   
   C2 = dfs->item( exp->string_data( "phase_field_2" ) ) ;
   L_C2 = exp->int_data( "level_2" ) ;
   
   PEL_ModuleExplorer const* se = 
                      exp->create_subexplorer( 0, "CH_BulkChemicalPotential" ) ;
   BULK_MU = CH_BulkChemicalPotential::create( this, se ) ;
   se->destroy() ; se = 0 ;

   cFE->require_field_calculation( C1, PDE_LocalFE::N )  ;   
   cFE->require_field_calculation( C1, PDE_LocalFE::dN ) ; 
   cFE->require_field_calculation( C2, PDE_LocalFE::N )  ;   
   cFE->require_field_calculation( C2, PDE_LocalFE::dN ) ; 
}

//----------------------------------------------------------------------------
void
CH_AveragesSaver:: read_for_kinetic_energy( PDE_DomainAndFields const* dom,
                                            FE_SetOfParameters const* prms,
                                            PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "CH_AveragesSaver:: read_for_kinetic_energy" ) ;
   
   UU = dom->set_of_discrete_fields()->item( exp->string_data( "velocity" ) ) ;
   L_UU = exp->int_data( "level_of_velocity" ) ;
   DENS = prms->item( exp->string_data( "param_density" ) ) ;
   
   QRP_KE = GE_QRprovider::object( 
                            exp->string_data( "quadrature_rule_provider" ) ) ;
   
   cFE->require_field_calculation( UU, PDE_LocalFE::N  ) ;
   DENS->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;

   //?????????? le faire pour les autres ?????
   check_field_nb_components( UU, cFE->nb_space_dimensions() ) ;
}

//----------------------------------------------------------------------------
void
CH_AveragesSaver:: read_for_center_volume( PDE_DomainAndFields const* dom,
                                           FE_SetOfParameters const* prms,
                                           PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "CH_AveragesSaver:: read_for_center_volume" ) ;
   
   CC = dom->set_of_discrete_fields()->item( 
                                       exp->string_data( "phase_field" ) ) ;
   L_CC = exp->int_data( "level_of_phase_field" ) ;
   VV = dom->set_of_discrete_fields()->item(
                                       exp->string_data( "velocity" ) ) ;
   L_VV = exp->int_data( "level_of_velocity" ) ;
   CC_MIN = exp->doubleVector_data( "thresholds" ) ;
   QRP_CV = GE_QRprovider::object( 
                            exp->string_data( "quadrature_rule_provider" ) ) ;
   
   cFE->require_field_calculation( CC, PDE_LocalFE::N  ) ;
   cFE->require_field_calculation( VV, PDE_LocalFE::N  ) ;

   //?????????? le faire pour les autres ?????
   check_field_nb_components( CC, 1 ) ;
   check_field_nb_components( VV, dom->nb_space_dimensions() ) ;
}

//----------------------------------------------------------------------------
void
CH_AveragesSaver:: read_for_perimeter( PDE_DomainAndFields const* dom,
                                       FE_SetOfParameters const* prms,
                                       PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "CH_AveragesSaver:: read_for_perimeter" ) ;
   
   CP = dom->set_of_discrete_fields()->item( 
                                       exp->string_data( "phase_field" ) ) ;
   L_CP = exp->int_data( "level_of_phase_field" ) ;
   
   QRP_PER = GE_QRprovider::object( 
                            exp->string_data( "quadrature_rule_provider" ) ) ;
   
   cFE->require_field_calculation( CP, PDE_LocalFE::dN  ) ;

   //?????????? le faire pour les autres ?????
   check_field_nb_components( CP, 1 ) ;
}

//----------------------------------------------------------------------------
void
CH_AveragesSaver:: save_all( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "CH_AveragesSaver:: save_all" ) ;
   
   size_t rank = communicator()->rank() ;
   std::ofstream ofs ;
   
   if( rank == 0 )
   {      
      ofs.open( FNAME.c_str(), std::ios::out | std::ios::app ) ;
      ofs.setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
      ofs << std::setprecision( 5 ) ;
      ofs << setw( 12 ) << t_it->time() << " " ;
   }
   if( TOTAL_FREE_ENERGY )
   {
      double tfe = total_free_energy() ;
      if( rank == 0 )
      {
         ofs << setw( 12 ) << tfe << " " ;
      }
   }
   if( KINETIC_ENERGY )
   {
      double ke =  kinetic_energy( t_it ) ;
      if( rank == 0 )
      {
         ofs << setw( 12 ) << ke << " " ;
      }
   }
   if( CENTER )
   {
      doubleVector volume( 0 ) ;
      doubleArray2D center( 0, 0 ) ;
      doubleArray2D velocity( 0, 0 ) ;
      size_t_vector nb_cells( 0 ) ;
      center_volume( volume, center, velocity, nb_cells ) ;
      if( rank == 0 ) 
      {
         for( size_t i=0 ; i<CC_MIN.size() ; ++i )
         {
            ofs << setw( 12 ) << volume( i ) << " " ;
         }
         if ( FE::geometry() == FE::axisymmetrical )
         {
            for( size_t i=0 ; i<CC_MIN.size() ; ++i )
            {
               ofs << setw( 12 ) << center( i, 1 ) << " " ;
            }         
            for( size_t i=0 ; i<CC_MIN.size() ; ++i )
            {
               ofs << setw( 12 ) << velocity( i, 1 ) << " " ;
            }         
         }
         else
         {
            for( size_t d=0 ; d<cFE->nb_space_dimensions() ; ++d )
            {
               for( size_t i=0 ; i<CC_MIN.size() ; ++i )
               {
                  ofs << setw( 12 ) << center( i, d ) << " " ;
               }
            }
            for( size_t d=0 ; d<cFE->nb_space_dimensions() ; ++d )
            {
               for( size_t i=0 ; i<CC_MIN.size() ; ++i )
               {
                  ofs << setw( 12 ) << velocity( i, d ) << " " ;
               }
            }
         }
      }
   }
   if( PERIMETER )
   {
      double per = perimeter() ;
      if( rank == 0 )
      {
         ofs << setw( 12 ) << per << " " ;
      }
   }
   if( rank == 0 )
   {
      ofs << endl ;
      ofs.close() ;
   }
}

//----------------------------------------------------------------------------
bool 
CH_AveragesSaver:: total_free_energy_is_handled( void ) const
//----------------------------------------------------------------------------
{
  return( TOTAL_FREE_ENERGY ) ;
}

//----------------------------------------------------------------------------
double
CH_AveragesSaver:: total_free_energy( void ) const
//----------------------------------------------------------------------------
{
   PEL_CHECK_PRE( total_free_energy_is_handled() ) ;   
   
   size_t nb_dims = cFE->nb_space_dimensions() ; 
   
   double s1 = BULK_MU->Sigma1() ;
   double s2 = BULK_MU->Sigma2() ;
   double s3 = BULK_MU->Sigma3() ;

   double result = 0.0 ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->start_IP_iterator( QRP_TFE ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         double weight = cFE->weight_of_IP() ;
         if( FE::geometry() == FE::axisymmetrical )
         {
            weight *= 2.0*PEL::pi()*cFE->coordinates_of_IP()->coordinate( 0 ) ;
         }

         double cc1 = 0.0 ;
         double cc2 = 0.0 ;
         double cc3 = 0.0 ;
         for( size_t d=0; d<nb_dims ; ++d )
         {
            double grad_c = cFE->gradient_at_IP( C1, L_C1, d ) ;
            double grad_d = cFE->gradient_at_IP( C2, L_C2, d ) ;

            cc1 += grad_c * grad_c ;
            cc2 += grad_d * grad_d ;
            cc3 += ( grad_c+grad_d ) * ( grad_c+grad_d ) ;
         }

         double c = cFE->value_at_IP( C1, L_C2 )  ;
         double d = cFE->value_at_IP( C2, L_C2 ) ;

         double free_e = 12.*BULK_MU->F(c,d,1.-c-d)/THICKNESS;

         result += weight * ( 3.*THICKNESS*( s1*cc1 + s2*cc2 + s3*cc3 )/8. 
                              + free_e ) ; 
      }
   } 

   result = communicator()->sum( result ) ;
   
   if( verbose_level() >= 2 )
   {
      PEL::out() << indent() << "   CH_energy = " << result << endl ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------------
bool 
CH_AveragesSaver:: kinetic_energy_is_handled( void ) const
//----------------------------------------------------------------------------
{
  return( KINETIC_ENERGY ) ;
}

//----------------------------------------------------------------------------
double
CH_AveragesSaver:: kinetic_energy( FE_TimeIterator const* t_it ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "CH_AveragesSaver:: kinetic_energy" ) ;
   PEL_CHECK_PRE( kinetic_energy_is_handled() ) ;   

   
   size_t nb_dims = cFE->nb_space_dimensions() ; 
   
   double result = 0.0 ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->start_IP_iterator( QRP_KE ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         double weight = cFE->weight_of_IP() ;
         if( FE::geometry() == FE::axisymmetrical )
         {
            weight *= 2.0*PEL::pi()*cFE->coordinates_of_IP()->coordinate( 0 ) ;
         }
         double dens = DENS->cell_value_at_IP( t_it, cFE ) ;
         double vv = 0.0 ;
         for( size_t d=0; d<nb_dims ; ++d )
         {
            double xx = cFE->value_at_IP( UU, L_UU, d ) ;
            vv += xx * xx ;
         }
         result += weight * vv * dens ;       
      }
   } 
   result *= 0.5 ;
   
   result = communicator()->sum( result ) ;

   if( verbose_level() >= 2 )
   {
      PEL::out() << indent() << "   Kinetic energy = " << result << endl ;
   }
   
   return( result ) ;
}

//----------------------------------------------------------------------------
bool 
CH_AveragesSaver:: center_volume_is_handled( void ) const
//----------------------------------------------------------------------------
{
  return( CENTER ) ;
}

//---------------------------------------------------------------------------
void
CH_AveragesSaver:: center_volume( doubleVector& volume,
                                  doubleArray2D& center,
                                  doubleArray2D& velocity,
                                  size_t_vector& nb_cells ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_AveragesSaver:: center_volume" ) ;
   PEL_CHECK_PRE( center_volume_is_handled() ) ;   

   size_t nb_dims = cFE->nb_space_dimensions() ;
   
   center.re_initialize( CC_MIN.size(), nb_dims ) ;
   velocity.re_initialize( CC_MIN.size(), nb_dims ) ;
   volume.re_initialize( CC_MIN.size() ) ;
   nb_cells.re_initialize( CC_MIN.size() ) ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      double cell_volume = 0.0 ;
      doubleVector cell_vv( nb_dims ) ;
      double cell_cc = 0.0 ;
      cFE->start_IP_iterator( QRP_CV ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         double cc = cFE->value_at_IP( CC, L_CC ) ;
         
         double weight = cFE->weight_of_IP() ;
         if( FE::geometry() == FE::axisymmetrical )
         {
            weight *= 2.0*PEL::pi()*cFE->coordinates_of_IP()->coordinate(0) ;
         }
         cell_volume += weight ;
         for( size_t d=0 ; d<nb_dims ; ++d )
         {
            cell_vv( d ) += weight * cFE->value_at_IP( VV, L_VV, d ) ;
         }
         
         cell_cc += weight * cc ;
      }

      for( size_t i=0 ; i<CC_MIN.size() ; ++i )
      {
         if( cell_cc >= CC_MIN( i ) * cell_volume )
         {
            ++nb_cells( i ) ;
            volume( i ) += cell_volume ;
            GE_Point const* pt = cFE->polyhedron()->center() ;
            for( size_t d=0; d<nb_dims ; ++d )
            {
               center( i, d ) += cell_volume * pt->coordinate( d ) ;
               velocity( i, d ) += cell_vv( d ) ; 
            }    
         }
      }
   }
   
   PEL_Communicator const* com = communicator() ;
   for( size_t i=0 ; i<CC_MIN.size() ; ++i )
   {
      volume( i ) = com->sum( volume( i ) ) ;
      for( size_t d=0; d<nb_dims ; ++d )
      {
         center( i, d ) = com->sum( center( i, d ) ) / volume( i )  ;
         velocity( i, d ) = com->sum( velocity( i, d ) ) / volume( i )  ;
      }      
   }
}

//----------------------------------------------------------------------------
bool 
CH_AveragesSaver:: perimeter_is_handled( void ) const
//----------------------------------------------------------------------------
{
  return( PERIMETER ) ;
}

//---------------------------------------------------------------------------
double
CH_AveragesSaver:: perimeter( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "CH_AveragesSaver:: perimeter" ) ;
   PEL_CHECK_PRE( perimeter_is_handled() ) ;   
   
   double result = 0 ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      double cell_per = 0.0 ;
      cFE->start_IP_iterator( QRP_PER ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         double weight = cFE->weight_of_IP() ;
         if( FE::geometry() == FE::axisymmetrical )
         {
            weight *= 2.0*PEL::pi()*cFE->coordinates_of_IP()->coordinate(0) ;
         }
         double grad_sqr = 0 ;
         for( size_t d = 0 ; d < cFE->nb_space_dimensions() ; ++d )
         {
            grad_sqr +=  PEL::sqr( cFE->gradient_at_IP( CP, L_CP, d ) ) ;
         }
         cell_per += weight * PEL::sqrt( grad_sqr ) ;
      }
      result += cell_per ;
   }
   
   result = communicator()->sum( result )  ;
   
   return( result ) ;
}

//internal--------------------------------------------------------------
void 
CH_AveragesSaver_ERROR:: n0( std::string const& fname  )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << endl << "*** CH_AveragesSaver_ERROR:" << endl << endl ;
   msg << "    unable to open file \"" << fname << "\"" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void 
CH_AveragesSaver_ERROR:: n1( std::string const& module_name, 
                             size_t invalid_dim,
                             size_t_vector valid_dims )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << endl << "*** CH_AveragesSaver_ERROR:" << endl << endl ;
   msg << "    Module " << module_name << " is not allowed in "
       << invalid_dim <<"D." << endl ;
   msg << "    allowed number of dimension : " ;
   for( size_t i = 0 ; i < valid_dims.size() ; ++ i )
   {
      msg << valid_dims( i ) << "D " ;
   }
   msg << endl ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}


