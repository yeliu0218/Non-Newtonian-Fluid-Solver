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

#include <FE_ComparatorWithAnalytic.hh>

#include <FE.hh>
#include <FE_TimeIterator.hh>
#include <doubleVector.hh>

#include <GE_Point.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_QRprovider.hh>
#include <GE_ReferencePolyhedron.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Communicator.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Root.hh>
#include <PEL_Variable.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfBCs.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using std::istringstream ;
using std::string ;
using std::endl ;
using std::ios_base ;
using std::setprecision ; using std::setw ;

struct FE_ComparatorWithAnalytic_ERROR
{
   static void n0( std::string const& field_name,
                   std::string const& error ) ;
   static void n1( std::string const& save_name ) ;
   static void n2( GE_Mpolyhedron const* poly ) ;
   static void n3( std::string const& fname ) ;
   static void n4( std::string const& norm_name, double pp ) ;
} ;

//----------------------------------------------------------------------
FE_ComparatorWithAnalytic const*
FE_ComparatorWithAnalytic:: PROTOTYPE = new FE_ComparatorWithAnalytic() ;
//----------------------------------------------------------------------

//----------------------------------------------------------------------
FE_ComparatorWithAnalytic:: FE_ComparatorWithAnalytic( void )
//----------------------------------------------------------------------
   : FE_OneStepIteration( "FE_ComparatorWithAnalytic" )
   , FIELD( 0 )
   , LEVEL( PEL::bad_index() )
   , cFE( 0 )
   , sFE( 0 )
   , bFE( 0 )
   , QRP( 0 )
   , NORMS( intVector( 0 ) )
   , EXPOS( doubleVector( 0 ) )
   , IDX_EXPOS( doubleVector( 0 ) )
   , NORM_NAMES( stringVector( 0 ) )
   , NORM_SAVE_NAMES( stringVector( 0 ) )
   , COORDS( 0 )
   , TT( 0 )
   , SOL( 0 )
   , DSOL( 0 )
   , NULLIFY_INTEGRAL( false )
   , DIST_REF( PEL::bad_double() )
   , COM( 0 )
   , SAVING_NUMBER( PEL::bad_index() )
   , NEXT_SAVING_TIME( PEL::bad_double() )
   , SAVING_TIMES( doubleVector( 0 ) )
{
}

//----------------------------------------------------------------------
FE_ComparatorWithAnalytic:: ~FE_ComparatorWithAnalytic( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
FE_ComparatorWithAnalytic*
FE_ComparatorWithAnalytic:: create_replica( PEL_Object* a_owner,
                                         PDE_DomainAndFields const* dom,
                                         FE_SetOfParameters const* prms,
                                         PEL_ModuleExplorer* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_ComparatorWithAnalytic* result =
      new FE_ComparatorWithAnalytic( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_ComparatorWithAnalytic:: FE_ComparatorWithAnalytic(
                                         PEL_Object* a_owner,
                                         PDE_DomainAndFields const* dom,
                                         FE_SetOfParameters const* prms,
                                         PEL_ModuleExplorer* exp )
//----------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , FIELD( dom->set_of_discrete_fields()->item(
                                         exp->string_data( "field" ) ) )
   , LEVEL( (size_t) exp->int_data( "level" ) )
   , cFE( dom->create_LocalFEcell( this ) )
   , sFE( dom->create_CursorFEside( this ) )
   , bFE( dom->create_LocalFEbound( this ) )
   , QRP( GE_QRprovider::object(
           exp->string_data( "quadrature_rule_provider" ) ) )
   , BCs( dom->set_of_boundary_conditions() )
   , NORMS( intVector( 0 ) )
   , EXPOS( doubleVector( 0 ) )
   , IDX_EXPOS( doubleVector( 0 ) )
   , NORM_NAMES( exp->stringVector_data( "norms" ) )
   , NORM_SAVE_NAMES( exp->stringVector_data( "norm_saving_names" ) )
   , COORDS( 0 )
   , TT( 0 )
   , SOL( 0 )
   , DSOL( 0 )
   , NULLIFY_INTEGRAL(
            exp->has_entry( "nullify_integral" ) ?
               exp->bool_data( "nullify_integral" ) : false )
   , DIST_REF(
      exp->has_entry( "reference_distance" ) ?
      exp->double_data( "reference_distance" ) : 0.01 )
   , COM( PEL_Exec::communicator() )
   , FILE_OUT( exp->has_entry( "output_file" ) ? 
               exp->string_data( "output_file" ) : "" )
   , SAVING_NUMBER( PEL::bad_index() )
   , NEXT_SAVING_TIME( PEL::bad_double() )
   , SAVING_TIMES( exp->has_entry( "saving_times" ) ? 
                      exp->doubleVector_data( "saving_times" ) : 
                   doubleVector( 0 ) )
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: FE_ComparatorWithAnalytic" ) ;
   PEL_CHECK_INV( invariant() ) ;

   check_field_storage_depth( FIELD, LEVEL ) ;

   // Preprocessing :
   if( NULLIFY_INTEGRAL && FIELD->nb_components() != 1 )
   {
      std::string const msg =
         "    nullify integral for 1 component\n"
         "    field only." ;
      FE_ComparatorWithAnalytic_ERROR::n0( FIELD->name(), msg ) ;
   }

   // Context :
   PEL_ContextSimple* ct = PEL_ContextSimple::create( this ) ;
   COORDS = PEL_DoubleVector::create(
      ct, doubleVector( dom->nb_space_dimensions() ) ) ;
   ct->extend( PEL_Variable::object( "DV_X" ), COORDS ) ;
   TT = PEL_Double::create( ct, 0.0 ) ;
   ct->extend( PEL_Variable::object( "DS_T" ), TT ) ;
   
   // Norms :
   if( NORM_NAMES.size()!=NORM_SAVE_NAMES.size() )
   {
      std::string const msg =
         "    \"norms\" and \"norm_saving_names\"\n"
         "    need to have the same dimension." ;
      FE_ComparatorWithAnalytic_ERROR::n0( FIELD->name(), msg ) ;
   }
   for( size_t i=0 ; i<NORM_NAMES.size() ; ++i )
   {
      for( size_t j=i+1 ; j<NORM_NAMES.size() ; ++j )
      {
         if( NORM_NAMES(j) == NORM_NAMES(i) )
         {
            std::string const msg =
               "    norm \""+NORM_NAMES(i)+"\" is defined several times." ;
            FE_ComparatorWithAnalytic_ERROR::n0( FIELD->name(), msg ) ;
         }
         if( NORM_SAVE_NAMES(j) == NORM_SAVE_NAMES(i) )
         {
            std::string const msg =
               "    norm saving names \""+NORM_SAVE_NAMES(i)
                                       +"\" is defined several times." ;
            FE_ComparatorWithAnalytic_ERROR::n0( FIELD->name(), msg ) ;
         }
      }
   }
   
   size_t nb_norms = NORM_NAMES.size() ;
   NORMS.re_initialize( nb_norms ) ;
   EXPOS.re_initialize( nb_norms ) ;
   
   MAX_LENN = 0 ;
   for( size_t i=0 ; i<nb_norms ; ++i )
   {
      NormType nt = Invalid ;
      read_norm( NORM_NAMES( i ), nt, EXPOS( i ), exp ) ;
      add_exponent( EXPOS( i ) ) ;
      NORMS( i ) = nt ;
      if( NORM_NAMES( i ).length() > MAX_LENN )
      {
         MAX_LENN = NORM_NAMES( i ).length() ;
      }
   }
            
   // Solutions :
   if( exp->has_entry( "solution" ) )
   {
      SOL = exp->abstract_data( this, "solution", ct ) ;
      cFE->require_field_calculation( FIELD, PDE_LocalFE::N ) ;
      sFE->require_field_calculation( FIELD, PDE_LocalFE::node ) ;
      bFE->require_field_calculation( FIELD, PDE_LocalFE::node ) ;
   }
   if( exp->has_entry( "d_solution" ) )
   {
      DSOL = exp->abstract_data( this, "d_solution", ct ) ;
      cFE->require_field_calculation( FIELD, PDE_LocalFE::dN ) ;
   }
   for( size_t i=0 ; i<NORMS.size() ; ++i )
   {
      NormType nt = (NormType) NORMS(i) ;
      if( nt == err_Linf || nt == sol_Linf || nt == err_Lp || nt == sol_Lp ||
          nt == eip_Lp || nt == err_LpD || nt == err_W1p || nt == sol_W1p || 
          nt == eip_W1p || nt == err_H1D || nt == err_H1DirD )
      {
         if( SOL==0 )
         {
            std::string const msg =
               "    \"solution\" expected to compute the requested norms" ;
            FE_ComparatorWithAnalytic_ERROR::n0( FIELD->name(), msg ) ;
         }
      }
      else if( nt == err_sW1p || nt == sol_sW1p || nt == eip_sW1p || 
               nt == err_W1p || nt == sol_W1p || nt == eip_W1p )
      {
         if( DSOL==0 )
         {
            std::string const msg =
               "    \"d_solution\" expected  to compute the requested norms" ;
            FE_ComparatorWithAnalytic_ERROR::n0( FIELD->name(), msg ) ;
         }
      }
   }
   
   // Check solution :
   if( SOL!=0 )
   {
      if( !SOL->value_can_be_evaluated() )
      {
         PEL_Error::object()->raise_not_evaluable(
                        exp, "solution", SOL->undefined_variables() ) ;
      }
      if( SOL->data_type()!=PEL_Data::DoubleVector )
      {
         PEL_Error::object()->raise_bad_data_type(
            exp, "solution", PEL_Data::DoubleVector ) ;
         
      }
   }
   if( DSOL!=0 )
   {
      if( !DSOL->value_can_be_evaluated() )
      {
         PEL_Error::object()->raise_not_evaluable(
                     exp, "d_solution", DSOL->undefined_variables() ) ;
      }
      if( DSOL->data_type()!=PEL_Data::DoubleArray2D )
      {
         PEL_Error::object()->raise_bad_data_type(
            exp, "d_solution", PEL_Data::DoubleArray2D ) ;
      }
   }

   if( exp->has_entry( "reference_distance" ) )
   {
      exp->test_data( "reference_distance", "reference_distance>0." ) ;
   }

   // Post-processing initialization:
   if( COM->rank()==0 && !FILE_OUT.empty() )
   {
      static stringVector FILE_NAMES( 0 ) ;
      if( !FILE_NAMES.has( FILE_OUT ) )
      {
         FILE_NAMES.append( FILE_OUT ) ;
         std::ofstream ofs ;
         ofs.open( FILE_OUT.c_str(), std::ios::out | std::ios::trunc ) ;
         if( !ofs ) FE_ComparatorWithAnalytic_ERROR::n3( FILE_OUT ) ;
         ofs << "#" << endl ;
         ofs << "# FE_ComparatorWithAnalytic generated file" << endl ;
         ofs << "# meaning of the columns:" << endl ;
         ofs << "#    TIME : current time" << endl ;
         ofs << "#    TIST : time step" << endl ;
         ofs << "#    XH   : maximum of the cell diameters" << endl ;
         ofs << "#    XD   : maximum of the cell equivalent ball diameters" 
             << endl ;
         ofs << "#    other columns : given by the data of keyword" << endl ;
         ofs << "#                    \"norm_saving_names\" in the data deck"
             << endl << "#" << endl ;
         ofs << "#       TIME         TIST           XH           XD " ;
         for( size_t i=0 ; i<NORMS.size() ; ++i )
         {
            ofs << setw( 12 ) << NORM_SAVE_NAMES(i) << " " ;
         }
         ofs << endl ;
         ofs.close() ;
      }
   }
}

//----------------------------------------------------------------------
void
FE_ComparatorWithAnalytic:: do_before_time_stepping(
                                           FE_TimeIterator const* t_it ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: do_before_time_stepping" ) ;
   PEL_CHECK_PRE( do_before_time_stepping_PRE( t_it ) ) ;

   if( SAVING_TIMES.size() != 0 )
   {
      if( !t_it->table_of_times_is_valid( SAVING_TIMES ) )
      {
         t_it->raise_invalid_table_of_times( "FE_ComparatorWithAnalytic", 
                                             "saving_times", 
                                             SAVING_TIMES ) ;
      }
      SAVING_NUMBER = 0 ;
      if( SAVING_TIMES(0) == t_it->time() )
      {
         save_norms( t_it, 0 ) ;
      }
      NEXT_SAVING_TIME = t_it->next_time_in_table( SAVING_TIMES ) ;
   }

   // Check solution dimensions :
   TT->set( t_it->time() ) ;
   cFE->start() ;
   COORDS->set( cFE->polyhedron()->center()->coordinate_vector() ) ;
   if( SOL!=0 )
   {
      doubleVector const& theo = SOL->to_double_vector() ;
      if( theo.size() != FIELD->nb_components() )
      {
         std::string const msg =
            "    Bad size for table \"solution\"." ;
         FE_ComparatorWithAnalytic_ERROR::n0( FIELD->name(), msg ) ;
      }
   }
   if( DSOL != 0 )
   {
      doubleArray2D d_theo = DSOL->to_double_array2D() ;
      if( d_theo.index_bound(0) != FIELD->nb_components() ||
          d_theo.index_bound(1) != cFE->nb_space_dimensions() )
      {
         std::string const msg =
            "    Bad size for table \"d_solution\"." ;
         FE_ComparatorWithAnalytic_ERROR::n0( FIELD->name(), msg ) ;
      }
   }
}

//----------------------------------------------------------------------
void
FE_ComparatorWithAnalytic:: do_one_inner_iteration(
                                           FE_TimeIterator const* t_it ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;
}

//----------------------------------------------------------------------
void
FE_ComparatorWithAnalytic:: do_after_time_adaptation(
                                           FE_TimeIterator const* t_it ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: do_after_time_adaptation" ) ;
   PEL_CHECK_PRE( do_after_time_adaptation_PRE( t_it ) ) ;
   
   if( SAVING_TIMES.size() != 0 )
   {
      if( FE_TimeIterator::greater_or_equal( t_it->time(),
                                             NEXT_SAVING_TIME ) )
      {
         start_total_timer( 
               "FE_ComparatorWithAnalytic:: do_after_time_adaptation" ) ;
         
         save_norms( t_it, 0 ) ;
         NEXT_SAVING_TIME = t_it->next_time_in_table( SAVING_TIMES ) ;
         
         stop_total_timer() ;      
      }
   }
}

//----------------------------------------------------------------------
void
FE_ComparatorWithAnalytic:: save_other_than_time_and_fields(
                                                 FE_TimeIterator const* t_it,  
                                                 PDE_ResultSaver* rs )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: save_other_than_time_and_fields" ) ;
   PEL_CHECK_PRE( save_other_than_time_and_fields_PRE( t_it, rs ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( SAVING_TIMES.size() == 0 )
   {
      start_total_timer( 
            "FE_ComparatorWithAnalytic:: save_other_than_time_and_fields" ) ;
      
      save_norms( t_it, rs ) ;
      
      stop_total_timer() ;
   }      
}

//----------------------------------------------------------------------
void
FE_ComparatorWithAnalytic:: print( std::ostream& os,
                                   size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;

   FE_OneStepIteration::print( os, indent_width ) ;
   
   std::string const s( indent_width+3, ' ' ) ;
   os << s << "field : \"" << FIELD->name() << "\"" << std::endl ;
   os << s << "time step saved in \"TIST\"" << std::endl ;
   os << s << "h (inter vertices distance)  saved in \"XH\"" << std::endl ;
   os << s << "h (equivalent ball diameter) saved in \"XD\"" << std::endl ;
   if( NULLIFY_INTEGRAL )
   {
      os << s << "nullify integral" << std::endl ;
   }
   
   size_t lll = 0 ;
   for( size_t i=0 ; i<NORM_SAVE_NAMES.size() ; ++i )
   {
      size_t ll = NORM_SAVE_NAMES( i ).length() ;
      if( ll > lll ) lll = ll ;
   }
   for( size_t i=0 ; i<NORMS.size() ; ++i )
   {
      os << s << setw( MAX_LENN+2) << ( "\"" + NORM_NAMES( i ) +"\"") 
         << "  saved in  " 
         << setw( lll+2 ) << ( "\"" + NORM_SAVE_NAMES(i) + "\"" ) ; 
      if( EXPOS( i ) != PEL::bad_double() )
      {
         os << " (exponent: " << EXPOS( i ) << ")" ;
      }
      os << std::endl ;
   }
}

//----------------------------------------------------------------------
void
FE_ComparatorWithAnalytic:: save_norms( FE_TimeIterator const* t_it,  
                                        PDE_ResultSaver* rs )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: save_norms" ) ;

   bool has_ofs = false ;
   std::ofstream ofs ;
   if( COM->rank()==0 && !FILE_OUT.empty() )
   {
      ofs.open( FILE_OUT.c_str(), std::ios::out | std::ios::app ) ;
      ofs.setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
      ofs << std::setprecision( 5 ) ;
      has_ofs = true ;
   }

   double const h = max_inter_vertices_distance_of_cells() ;
   double const d = max_equivalent_ball_diameter_of_cells() ;
   if( rs != 0 )
   {
      if( !rs->has_variable( "TIST" ) )
      {
         rs->save_variable( t_it->time_step(), "TIST" ) ;
      }
      if( !rs->has_variable( "XH" ) )
      {
         rs->save_variable( h, "XH" ) ;
      }
      if( !rs->has_variable( "XD" ) )
      {
         rs->save_variable( d, "XD" ) ;
      }
   }
   if( has_ofs )
   {
      ofs << setw( 12 ) << t_it->time() << " " ;
      ofs << setw( 12 ) << t_it->time_step() << " " ;
      ofs << setw( 12 ) << h << " " ;
      ofs << setw( 12 ) << d << " " ;
   }

   if( verbose_level() != 0 )
   {
      increase_indent() ;
      if( t_it->is_started() )
      {
         PEL::out() << indent() << "time step : " << t_it->time_step()
                    << std::endl ;
      }
      PEL::out() << indent() << "h (inter vertices distance)  : "
                 << h << std::endl ;
      PEL::out() << indent() << "h (equivalent ball diameter) : "
                 << d << std::endl ;
      decrease_indent() ;
   }

   doubleVector norm_values( NORMS.size() ) ;
   compute_norms( t_it, norm_values ) ;
   for( size_t i=0 ; i<NORMS.size() ; ++i )
   {
      if( rs != 0 )
      {
         if( rs->has_variable( NORM_SAVE_NAMES(i) ) )
            FE_ComparatorWithAnalytic_ERROR::n1( NORM_SAVE_NAMES(i) ) ;
         rs->save_variable( norm_values(i), NORM_SAVE_NAMES(i) ) ;
      }
      if( has_ofs ) ofs << setw( 12 ) << norm_values(i) << " " ;
   }

   if( has_ofs ) 
   {
      ofs << endl ;
      ofs.close() ;
   }
}

//----------------------------------------------------------------------
void
FE_ComparatorWithAnalytic:: compute_norms( FE_TimeIterator const* t_it,
                                           doubleVector& norm_values )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: compute_norms" ) ;
   PEL_CHECK( t_it != 0 ) ;
   PEL_CHECK( norm_values.size() == NORMS.size() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   TT->set( t_it->time() ) ;

   double sol_shift = 0.0 ;
   double app_shift = 0.0 ;
   compute_shift( sol_shift, app_shift ) ;
   
   double x_sol_Linf = PEL::bad_double() ;
   double x_err_Linf = PEL::bad_double() ;
   
   size_t nbn = nb_exponents() ;
   doubleVector x_sol_Lp( nbn ) ; x_sol_Lp.set( PEL::bad_double() ) ;
   doubleVector x_eip_Lp( nbn ) ; x_eip_Lp.set( PEL::bad_double() ) ;
   doubleVector x_err_Lp( nbn ) ; x_err_Lp.set( PEL::bad_double() ) ;

   doubleVector x_err_LpD( nbn ) ; x_err_LpD.set( PEL::bad_double() ) ;

   doubleVector x_sol_sW1p( nbn ) ; x_sol_sW1p.set( PEL::bad_double() ) ;
   doubleVector x_eip_sW1p( nbn ) ; x_eip_sW1p.set( PEL::bad_double() ) ;
   doubleVector x_err_sW1p( nbn ) ; x_err_sW1p.set( PEL::bad_double() ) ;
   
   double x_err_H1D = PEL::bad_double() ;
   double x_err_H1DirD = PEL::bad_double() ;
   
   for( size_t i=0 ; i<NORMS.size() ; ++i )
   {
      NormType const norm = (NormType) NORMS( i ) ;
      double pp = EXPOS( i ) ;
      size_t p = idx_of_exponent( pp ) ;
      if( ( norm == sol_Linf || norm == err_Linf ) &&  
          x_sol_Linf == PEL::bad_double() )
      {
         calc_Linf( sol_shift, app_shift, x_sol_Linf, x_err_Linf ) ;
      }
      else if( ( norm == sol_Lp || norm == eip_Lp || norm == err_Lp ) &&
               x_sol_Lp( p ) == PEL::bad_double() )
      {
         calc_Lp( sol_shift, app_shift, pp, 
                  x_sol_Lp( p ), x_eip_Lp( p ), x_err_Lp( p ) ) ;
      }
      else if( ( norm == err_LpD ) && 
               x_err_LpD( p ) == PEL::bad_double() )
      {
         calc_LpD( sol_shift, app_shift, pp, x_err_LpD( p ) ) ;
      }
      else if( ( norm == sol_sW1p || norm == eip_sW1p || norm == err_sW1p ) &&
               x_sol_sW1p( p ) == PEL::bad_double() )
      {
         calc_sW1p( pp, x_sol_sW1p( p ), x_eip_sW1p( p ), x_err_sW1p( p ) ) ;
      }
      else if( norm == sol_W1p || norm == eip_W1p || norm == err_W1p )
      {
         if( x_sol_Lp( p ) == PEL::bad_double() )
         {
            calc_Lp( sol_shift, app_shift, pp, 
                     x_sol_Lp( p ), x_eip_Lp( p ), x_err_Lp( p ) ) ;
         }
         if( x_sol_sW1p( p ) == PEL::bad_double() )
         {
            calc_sW1p( pp, 
                       x_sol_sW1p( p ), x_eip_sW1p( p ), x_err_sW1p( p ) ) ;
         }
      }
   }
   
   for( size_t i=0 ; i<NORMS.size() ; ++i )
   {
      NormType const norm = (NormType) NORMS( i ) ;
      double pp = EXPOS( i ) ;
      size_t p = idx_of_exponent( pp ) ;
      norm_values( i ) = -PEL::max_double() ;
      switch( norm )
      {
         case err_Linf :
            norm_values( i ) = x_err_Linf ;
            break ;
         case sol_Linf :
            norm_values( i ) = x_sol_Linf ;
            break ;
         case err_Lp :
            norm_values( i ) = PEL::pow( x_err_Lp( p ), 1.0/pp ) ;
            break ;
         case err_LpD :
            norm_values( i ) = PEL::pow( x_err_LpD( p ), 1.0/pp ) ;
            break ;
         case sol_Lp :
            norm_values( i ) = PEL::pow( x_sol_Lp( p ), 1.0/pp ) ;
            break ;
         case eip_Lp :
            norm_values( i ) = PEL::pow( x_eip_Lp( p ), 1.0/pp ) ;
            break ;
         case err_sW1p :
            norm_values( i ) = PEL::pow( x_err_sW1p( p ), 1.0/pp ) ;
            break ;
         case err_H1D :
            if( x_err_H1D == PEL::bad_double() )
            {
               calc_H1D( x_err_H1D ) ;
            }
            norm_values( i ) = PEL::sqrt( x_err_H1D ) ;
            break ;
         case err_H1DirD :
            if( x_err_H1DirD == PEL::bad_double() )
            {
               calc_H1DirD( x_err_H1DirD ) ;
            }
            norm_values( i ) = PEL::sqrt( x_err_H1DirD ) ;
            break ;
         case sol_sW1p :
            norm_values( i ) = PEL::pow( x_sol_sW1p( p ), 1.0/pp ) ;
            break ;
         case eip_sW1p :
            norm_values( i ) = PEL::pow( x_eip_sW1p( p ), 1.0/pp  ) ;
            break ;
         case err_W1p :
            norm_values( i ) = PEL::pow( x_err_Lp( p ) + x_err_sW1p( p ), 
                                         1.0/pp ) ;
            break ;
         case sol_W1p :
            norm_values( i ) = PEL::pow( x_sol_Lp( p ) + x_sol_sW1p( p ), 
                                         1.0/pp ) ;
            break ;
         case eip_W1p :
            norm_values( i ) = PEL::pow( x_eip_Lp( p ) + x_eip_sW1p( p ), 
                                         1.0/pp ) ;
            break ;
         default :
            PEL_Error::object()->raise_internal(
               "*** FE_ComparatorWithAnalytic error :\n"
               "    try to use a norm which is not implemented." ) ;
      }
   }

   if( verbose_level() != 0 )
   {
      ios_base::fmtflags original_flags = PEL::out().flags() ;
      PEL::out().setf( ios_base::uppercase | ios_base::scientific ) ;
      increase_indent() ;
      PEL::out() << indent() << "field : \"" << FIELD->name() << "\""
                 << std::endl ;
      for( size_t i=0 ; i<NORM_NAMES.size() ; ++i )
      {
         PEL::out() << indent() << setw( MAX_LENN+1) << NORM_NAMES( i )
                    << " : " << setprecision( 8 ) << setw( 15 ) 
                    << norm_values( i ) << std::endl ;
      }
      decrease_indent() ;
      PEL::out().flags( original_flags ) ;
   }
}

//----------------------------------------------------------------------
void
FE_ComparatorWithAnalytic:: calc_Linf( double sol_shift, double app_shift,
                                       double& sol_norm, double& err_norm )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: calc_Linf" ) ;
   PEL_CHECK_INV( invariant() ) ;

   err_norm = 0.0 ;
   sol_norm = 0.0 ;

   size_t const nbc = FIELD->nb_components() ;
   
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      for( size_t loc_n=0 ; loc_n < cFE->nb_local_nodes( FIELD ) ; ++loc_n )
      {
         if( cFE->local_node_is_in_mesh( FIELD, loc_n ) )
         {
            cFE->set_calculation_point( 
                              cFE->local_node_location( FIELD, loc_n ) ) ;
            doubleVector const& dv = 
                              cFE->calculation_point()->coordinate_vector() ;
            COORDS->set( dv ) ;
            doubleVector const& theo = SOL->to_double_vector() ;
            for( size_t ic=0 ; ic<nbc ; ++ic   )
            {
               double app = cFE->value_at_pt( FIELD, LEVEL, ic )-app_shift ;
               double sol = theo( ic ) - sol_shift ;
               err_norm = PEL::max( err_norm, PEL::abs( sol - app ) ) ;
               sol_norm = PEL::max( sol_norm, PEL::abs( sol       ) ) ;
            }
         }
      }
   }
   err_norm = COM->max( err_norm ) ;
   sol_norm = COM->max( sol_norm   ) ;
}

//----------------------------------------------------------------------
void
FE_ComparatorWithAnalytic:: calc_Lp( 
                 double sol_shift, double app_shift, double exponent,                    
                 double& sol_norm, double& eip_norm, double& err_norm )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: calc_Lp" ) ;
   PEL_CHECK_INV( invariant() ) ;

   sol_norm = 0.0 ;
   eip_norm = 0.0 ; 
   err_norm = 0.0 ;

   size_t const nbc = FIELD->nb_components() ;
   doubleVector inter( nbc ) ;
 
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->set_row_and_col_fields( FIELD, FIELD ) ;
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         inter.set( 0.0 ) ;
         doubleVector const& N = cFE->Ns_at_IP( PDE_LocalFE::row ) ;
         for( size_t i=0 ; i<N.size() ; ++i )
         {
            GE_Point const* node = cFE->local_node_location( FIELD, i ) ;
            COORDS->set( node->coordinate_vector() ) ;
            doubleVector const& sol_at_node = SOL->to_double_vector() ;
            for( size_t ic=0 ; ic<nbc ; ++ic   )
            {
               //????? a sortir de la boucle ????????????????????????????,
               inter( ic ) += sol_at_node( ic ) * N( i ) ;
            }
         }

         doubleVector const& dv =
                             cFE->coordinates_of_IP()->coordinate_vector() ;
         COORDS->set( dv ) ;
         double weight = cFE->weight_of_IP() ;
         if( FE::geometry() == FE::axisymmetrical )
         {
            weight *= dv( 0 ) ;
         }
         doubleVector const& theo = SOL->to_double_vector() ;
         for( size_t ic=0 ; ic<nbc ; ++ic   )
         {
            double app = cFE->value_at_IP( FIELD, LEVEL, ic ) - app_shift ;
            double sol = theo( ic ) - sol_shift ;
            double eip = sol - inter( ic ) ;
            double err = sol - app ;
            sol_norm += weight * PEL::pow( PEL::abs( sol ), exponent ) ;
            eip_norm += weight * PEL::pow( PEL::abs( eip ), exponent ) ;
            err_norm += weight * PEL::pow( PEL::abs( err ), exponent ) ;
         }
      }
   }
   sol_norm = COM->sum( sol_norm ) ;
   eip_norm = COM->sum( eip_norm ) ;
   err_norm = COM->sum( err_norm ) ;
}

//----------------------------------------------------------------------
void
FE_ComparatorWithAnalytic:: calc_sW1p( double exponent,
                                       double& sol_snorm,
                                       double& eip_snorm,
                                       double& err_snorm )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: calc_sW1p" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   sol_snorm = 0.0 ;
   eip_snorm = 0.0 ;
   err_snorm = 0.0 ;

   size_t const nbc = FIELD->nb_components() ;
   size_t const nb_dims = cFE->nb_space_dimensions() ;
   doubleArray2D inter( nbc, nb_dims ) ;
  
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      cFE->set_row_and_col_fields( FIELD, FIELD ) ;
      cFE->start_IP_iterator( QRP ) ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      { 
         inter.set( 0.0 ) ;
         doubleArray2D const& dN = cFE->dNs_at_IP( PDE_LocalFE::row ) ;
         for( size_t i=0 ; i<dN.index_bound(0) ; ++i )
         {
            GE_Point const* node = cFE->local_node_location( FIELD, i ) ;
            COORDS->set( node->coordinate_vector() ) ;
            doubleVector const& sol_at_node = SOL->to_double_vector() ;
            for( size_t ic=0 ; ic<nbc ; ++ic   )
            {
               for( size_t d=0 ; d<nb_dims ; ++d )
               {
                  inter( ic, d ) += sol_at_node( ic )*dN( i, d ) ;
               }
            }
         }

         doubleVector const& dv =
                             cFE->coordinates_of_IP()->coordinate_vector() ;
         COORDS->set( dv ) ;
         double weight = cFE->weight_of_IP() ;
         if( FE::geometry() == FE::axisymmetrical )
         {
            weight *= dv( 0 ) ;
         }
         doubleArray2D theo = DSOL->to_double_array2D() ;
         for( size_t ic=0 ; ic<nbc ; ++ic   )
         {
            for( size_t d=0 ; d<nb_dims ; ++d )
            {
               double app = cFE->gradient_at_IP( FIELD, LEVEL, d, ic ) ;
               double sol = theo( ic, d ) ;
               double eip = sol - inter( ic, d ) ;
               double err = sol - app ;
               sol_snorm += weight * PEL::pow( PEL::abs( sol ), exponent ) ;
               eip_snorm += weight * PEL::pow( PEL::abs( eip ), exponent ) ;
               err_snorm += weight * PEL::pow( PEL::abs( err ), exponent ) ;
            }
         }
      }
   }
   sol_snorm = COM->sum( sol_snorm   ) ;
   eip_snorm = COM->sum( eip_snorm ) ;
   err_snorm = COM->sum( err_snorm ) ;
}

//----------------------------------------------------------------------
void
FE_ComparatorWithAnalytic:: calc_LpD( double sol_shift, double app_shift,
                                      double exponent,
                                      double& err_norm )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: calc_LpD" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   err_norm = 0.0 ;

   size_t const nbc = FIELD->nb_components() ;
 
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      GE_Point const* pt = D_norm_center( cFE->polyhedron() ) ;
      
      double vol = FE::cell_measure( cFE ) ;

      COORDS->set( pt->coordinate_vector() ) ;
      doubleVector const& theo = SOL->to_double_vector() ;

      check_nb_local_nodes( FIELD, cFE, 1 ) ;
      size_t n = cFE->global_node( FIELD, 0 ) ;

      for( size_t ic=0 ; ic<nbc ; ++ic   )
      {
         double app = FIELD->DOF_value( LEVEL, n, ic ) - app_shift ;
         double sol = theo( ic ) - sol_shift ;
         double err = sol - app ;
         err_norm += vol * PEL::pow( PEL::abs( err ), exponent ) ;
      }
   }
   err_norm = COM->sum( err_norm ) ;
}

//----------------------------------------------------------------------
void
FE_ComparatorWithAnalytic:: calc_H1D( double& err_norm )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: calc_H1D" ) ;

   size_t const nbc = FIELD->nb_components() ;

   err_norm = 0.0 ;
   for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
   {
      PDE_LocalFEcell const* cFE0 = sFE->adjacent_localFEcell( 0 ) ;
      PDE_LocalFEcell const* cFE1 = sFE->adjacent_localFEcell( 1 ) ;

      GE_Point const* pt0 = cFE0->polyhedron()->finite_volume_center() ;
      if( pt0 == 0 ) FE_ComparatorWithAnalytic_ERROR::n2( cFE0->polyhedron() ) ;
      COORDS->set( pt0->coordinate_vector() ) ;
      doubleVector th0 = SOL->to_double_vector() ;

      GE_Point const* pt1 = cFE1->polyhedron()->finite_volume_center() ;
      if( pt1 == 0 ) FE_ComparatorWithAnalytic_ERROR::n2( cFE1->polyhedron() ) ;
      COORDS->set( pt1->coordinate_vector() ) ;
      doubleVector th1 = SOL->to_double_vector() ;

      check_nb_local_nodes( FIELD, cFE0, 1 ) ;
      check_nb_local_nodes( FIELD, cFE1, 1 ) ;
      size_t n_0 = cFE0->global_node( FIELD, 0 ) ;
      size_t n_1 = cFE1->global_node( FIELD, 0 ) ;

      double h = sFE->distance_to_adjacent_finite_volume_center( 0 ) +
                 sFE->distance_to_adjacent_finite_volume_center( 1 ) ;
      double area = FE::side_measure( sFE ) ;

      for( size_t ic=0 ; ic<nbc ; ++ic )
      {
         double app0 = FIELD->DOF_value( LEVEL, n_0, ic ) ;
         double app1 = FIELD->DOF_value( LEVEL, n_1, ic ) ;
         double err = ( app1 - th1( ic ) ) - ( app0 - th0( ic ) ) ;
         err_norm += area/h * PEL::sqr( err ) ;
      }
   }
   err_norm = COM->sum( err_norm ) ;
}

//----------------------------------------------------------------------
void
FE_ComparatorWithAnalytic:: calc_H1DirD( double& err_norm )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: calc_H1DirD" ) ;

   calc_H1D( err_norm ) ;

   size_t const nbc = FIELD->nb_components() ;
   
   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      check_nb_local_nodes( FIELD, bFE, 1 ) ;
      size_t n = bFE->global_node( FIELD, 0 ) ;

      GE_Color const* color = bFE->color() ;
      double area = FE::bound_measure( bFE ) ;

      cFE->go_i_th( bFE->adjacent_cell_id() ) ;
      GE_Point const* pt = cFE->polyhedron()->finite_volume_center() ;
      if( pt == 0 ) FE_ComparatorWithAnalytic_ERROR::n2( cFE->polyhedron() ) ;
      COORDS->set( pt->coordinate_vector() ) ;
      doubleVector th = SOL->to_double_vector() ;

      for( size_t ic=0 ; ic<nbc ; ++ic )
      {
         if( BCs->has_BC( color, FIELD, ic ) )
         {
            PEL_ModuleExplorer const* ee = 
                                      BCs->BC_explorer( color, FIELD, ic ) ;
            string const& type = ee->string_data( "type" ) ;
            if( type == "Dirichlet" )
            {
               double h = bFE->distance_to_adjacent_finite_volume_center() ;
               double err = FIELD->DOF_value( LEVEL, n, ic ) - th( ic ) ;
               err_norm += area/h * PEL::sqr( err ) ;
            }
         }
      }
   }
   err_norm = COM->sum( err_norm ) ;
}

//----------------------------------------------------------------------
void
FE_ComparatorWithAnalytic:: compute_shift( double& sol_shift,
                                           double& app_shift )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: compute_shift" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   sol_shift = 0. ;
   app_shift = 0. ;

   if( NULLIFY_INTEGRAL )
   {
      double w_tot = 0. ;
      for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
      {
         cFE->start_IP_iterator( QRP ) ;
         for( ; cFE->valid_IP() ; cFE->go_next_IP() )
         {
            doubleVector const& dv =
               cFE->coordinates_of_IP()->coordinate_vector() ;
            COORDS->set( dv ) ;
            double weight = cFE->weight_of_IP() ;
            if( FE::geometry() == FE::axisymmetrical )
            {
               weight *= dv(0) ;
            }
            w_tot += weight ;
            double const val = cFE->value_at_IP( FIELD, LEVEL ) ;
            doubleVector const& theo = SOL->to_double_vector() ;
            double const sol = theo(0) ;
            sol_shift += weight*sol ;
            app_shift += weight*val ;
         }
      }
      sol_shift = COM->sum( sol_shift ) ;
      app_shift = COM->sum( app_shift ) ;
      w_tot = COM->sum( w_tot ) ;
      sol_shift /= w_tot ;
      app_shift /= w_tot ;
   }
}

//----------------------------------------------------------------------
double
FE_ComparatorWithAnalytic:: max_inter_vertices_distance_of_cells( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: max_inter_vertices_distance_of_cells" ) ;
   double result = 0. ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      double h = cFE->polyhedron()->inter_vertices_maximum_distance() ;
      if( h > result ) result = h ;      
   }
   result = COM->max( result ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
double
FE_ComparatorWithAnalytic:: max_equivalent_ball_diameter_of_cells( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: max_equivalent_ball_diameter_of_cells" ) ;
   double result = 0. ;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      double h = cFE->polyhedron()->equivalent_ball_diameter() ;
      if( h > result ) result = h ;      
   }
   result = COM->max( result ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_Point const*
FE_ComparatorWithAnalytic:: D_norm_center(
                                      GE_Mpolyhedron const* poly ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: D_norm_center" ) ;
   PEL_CHECK( poly != 0 ) ;
   PEL_CHECK( poly->dimension() == poly->nb_space_dimensions() ) ;
   
   static GE_Point* pt_ref =
      GE_Point::create( PEL_Root::object(), poly->nb_space_dimensions() ) ;
   GE_Point const* result = poly->finite_volume_center() ;
   if( result != 0 )
   {
      poly->apply_inverse_mapping( result, pt_ref ) ;
      if( !poly->reference_polyhedron()->contains( pt_ref, -DIST_REF ) )
      {
         result = 0 ;
      }
   }
   if( result == 0 ) result = poly->center() ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->nb_coordinates() == poly->nb_space_dimensions() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
FE_ComparatorWithAnalytic:: read_norm( std::string const& norm_name,
                                       NormType& nt, double& exponent,
                                       PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: read_norm" ) ;
   
   nt = Invalid ;
   exponent = PEL::bad_double() ;
   
   if( norm_name == "Infinity_error_norm" )
   {
      nt = err_Linf ;
   }
   else if( norm_name == "Infinity_solution_norm" )
   {
      nt = sol_Linf ;
   }
   else if( norm_name == "H1_0_error_norm" )
   {
      nt = err_sW1p ;
      exponent = 2.0 ;
   }
   else if( norm_name == "H1_D_error_D_norm" )
   {
      nt = err_H1D ;
   }
   else if( norm_name == "H1_Dirichlet_D_error_D_norm" )
   {
      nt = err_H1DirD ;
   }
   else if( norm_name == "H1_0_solution_norm" )
   {
      nt = sol_sW1p ;
      exponent = 2.0 ;
   }
   else if( norm_name == "H1_0_interpolation_error_norm" )
   {
      nt = eip_sW1p ;
      exponent = 2.0 ;
   }
   else if( norm_name == "H1_error_norm" )
   {
      nt = err_W1p ;
      exponent = 2.0 ;
   }
   else if( norm_name == "H1_solution_norm" )
   {
      nt = sol_W1p ;
      exponent = 2.0 ;
   }
   else if( norm_name == "H1_interpolation_error_norm" )
   {
      nt = eip_W1p ;
      exponent = 2.0 ;
   }
   else if( ( norm_name.length() > 5 ) && 
            ( norm_name[0] == 'L' ) )
   {
      string nn ;
      parse_norm_name( norm_name, "L", exponent, nn ) ;
      if( nn == "error_norm" )
      {
         nt = err_Lp ;
      }
      else if( nn == "error_D_norm" )
      {
         nt = err_LpD ;
      }
      else if( nn == "solution_norm" )
      {
         nt = sol_Lp ;
      }
      else if( nn == "interpolation_error_norm" )
      {
         nt = eip_Lp ;
      }
   }
   else if( ( norm_name.length() > 5 ) && 
            ( norm_name.substr( 0, 3 ) == "sW1" ) )
   {
      string nn ;
      parse_norm_name( norm_name, "sW1", exponent, nn ) ;
      if( nn == "error_norm" )
      {
         nt = err_sW1p ;
      }
      else if( nn == "solution_norm" )
      {
         nt = sol_sW1p ;
      }
      else if( nn == "interpolation_error_norm" )
      {
         nt = eip_sW1p ;
      }
   }
   else if( ( norm_name.length() > 5 ) && 
            ( norm_name.substr( 0, 2 ) == "W1" ) )
   {
      string nn ;
      parse_norm_name( norm_name, "W1", exponent, nn ) ;
      if( nn == "error_norm" )
      {
         nt = err_W1p ;
      }
      else if( nn == "solution_norm" )
      {
         nt = sol_W1p ;
      }
      else if( nn == "interpolation_error_norm" )
      {
         nt = eip_W1p ;
      }
   }
   
   if( nt == Invalid )
   {
      PEL_ASSERT( exponent == PEL::bad_index() ) ;
      PEL_Error::object()->raise_bad_data_value(
         exp, "norms",
         "   - \"Infinity_solution_norm\"\n"
         "   - \"Infinity_error_norm\"\n"
         "   - \"Lp_solution_norm\" (with a decent value for p)\n"
         "   - \"Lp_interpolation_error_norm\" (with a decent value for p)\n"
         "   - \"Lp_error_norm\" (with a decent value for p)\n"
         "   - \"sW1p_solution_norm\" (with a decent value for p)\n"
         "   - \"sW1p_interpolation_error_norm\" (with a decent value for p)\n"
         "   - \"sW1p_error_norm\" (with a decent value for p)\n"
         "   - \"W1p_solution_norm\" (with a decent value for p)\n"
         "   - \"W1p_interpolation_error_norm\" (with a decent value for p)\n"
         "   - \"W1p_error_norm\" (with a decent value for p)\n"
         "   - \"H1_0_error_norm\"\n"
         "   - \"H1_0_solution_norm\"\n"
         "   - \"H1_0_interpolation_error_norm\"\n"
         "   - \"H1_error_norm\"\n"
         "   - \"H1_solution_norm\"\n"
         "   - \"H1_interpolation_error_norm\""
         "   - \"Lp_error_D_norm\" (with a decent value for p)\n"
         "   - \"H1_D_error_D_norm\"\n"
         "   - \"H1_Dirichlet_D_error_D_norm\"\n"
         ) ;
   }
}

//----------------------------------------------------------------------
void
FE_ComparatorWithAnalytic:: parse_norm_name( std::string const& norm_name,
                                             std::string const& start,
                                             double& pp, std::string& end )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_ComparatorWithAnalytic:: parse_norm_name" ) ;
   
   size_t nl = norm_name.length() ;
   size_t sl = start.length() ;
   
   bool ok = ( nl > sl+1 ) ;
   if( ok )
   {
      size_t jj = sl + 1 ;
      do
      {
         if( norm_name[jj] == '_' ) break ;
         ++jj ;
      }
      while( jj <= nl ) ;
      ok = ( jj != nl ) ;
      if( ok )
      {
         istringstream iss( norm_name.substr( sl, jj-1 ) ) ;
         ok = ( iss >> pp ) ;
         if( ok )
         {
            end = norm_name.substr( jj+1, nl ) ;
         }
      }
   }
   
   if( pp < 1.0 ) FE_ComparatorWithAnalytic_ERROR:: n4( norm_name, pp ) ; 
}

//----------------------------------------------------------------------
void
FE_ComparatorWithAnalytic:: add_exponent( double pp )
//----------------------------------------------------------------------
{
   if( pp != PEL::bad_double() )
   {
      size_t idx = idx_of_exponent( pp ) ;
      if( idx == PEL::bad_index() )
      {
         IDX_EXPOS.append( pp ) ;
      }
   }
}

//----------------------------------------------------------------------
size_t
FE_ComparatorWithAnalytic:: nb_exponents( void ) const
//----------------------------------------------------------------------
{
   return( IDX_EXPOS.size() ) ;
}

//----------------------------------------------------------------------
size_t
FE_ComparatorWithAnalytic:: idx_of_exponent( double pp ) const
//----------------------------------------------------------------------
{
   size_t result = IDX_EXPOS.index_of( pp, 1.e-8, 1.e-14 ) ;
   return( result ) ;
}

//internal--------------------------------------------------------------
void 
FE_ComparatorWithAnalytic_ERROR:: n0( std::string const& field_name,
                                      std::string const& error )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << endl << "*** FE_ComparatorWithAnalytic:" << endl << endl ;
   msg << "    field: \"" << field_name << "\"" << endl ;
   msg << error ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void 
FE_ComparatorWithAnalytic_ERROR:: n1( std::string const& save_name )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << endl << "*** FE_ComparatorWithAnalytic:" << endl << endl ;
   msg << "    the saving name \"" << save_name << "\"" << endl ;
   msg << "    is multiply used" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void 
FE_ComparatorWithAnalytic_ERROR:: n2( GE_Mpolyhedron const* poly )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << endl << "*** FE_ComparatorWithAnalytic:" << endl << endl ;
   msg << "    the \"finite volume center\" has not been defined" << endl ;
   msg << "    for the following polyhedron:" << endl << endl ;
   poly->print( msg, 4 ) ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void 
FE_ComparatorWithAnalytic_ERROR:: n3( std::string const& fname  )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << endl << "*** FE_ComparatorWithAnalytic:" << endl << endl ;
   msg << "    unable to open file \"" << fname << "\"" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}

//internal--------------------------------------------------------------
void 
FE_ComparatorWithAnalytic_ERROR:: n4( std::string const& norm_name, double pp )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << endl << "*** FE_ComparatorWithAnalytic:" << endl << endl ;
   msg << "    the norm of name \"" << norm_name << "\"" << endl ;
   msg << "    is associated to L^{" << pp << "} which is" << endl ;
   msg << "    certainly an error" ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
