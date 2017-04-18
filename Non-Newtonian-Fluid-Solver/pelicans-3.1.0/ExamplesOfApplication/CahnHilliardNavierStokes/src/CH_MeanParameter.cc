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

#include <CH_MeanParameter.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_CursorFEside.hh>

#include <PDE_SetOfDiscreteFields.hh>

#include <FE.hh>
#include <FE_SetOfParameters.hh>
#include <iostream>

using std::string ;

CH_MeanParameter const* 
CH_MeanParameter:: PROTOTYPE = new CH_MeanParameter() ;

//-------------------------------------------------------------------------
CH_MeanParameter:: CH_MeanParameter( void )
//-------------------------------------------------------------------------
   : FE_Parameter( "CH_MeanParameter" )
   , PHs( 0 )
   , dPHs( 0 )
{
}

//-------------------------------------------------------------------------
CH_MeanParameter*
CH_MeanParameter:: create_replica( PEL_Object* a_owner,
                                   PDE_DomainAndFields const* dom,
                                   PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "CH_MeanParameter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   CH_MeanParameter* result = new CH_MeanParameter( a_owner, dom, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
CH_MeanParameter:: CH_MeanParameter( PEL_Object* a_owner,
                                     PDE_DomainAndFields const* dom,
                                     PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : FE_Parameter( a_owner, exp->string_data( "name" ) )
   , NB_PHASES( 0 )
   , PHs( 0 )
   , dPHs( 0 )
   , type( -1  )
{
   PDE_SetOfDiscreteFields const* dfs = dom->set_of_discrete_fields() ;
   PEL_ModuleExplorer* se = exp->create_subexplorer( 0, "phase_fields" ) ;
   se->start_module_iterator() ;
   for( ; se->is_valid_module() ; se->go_next_module() )
   {
      PEL_ModuleExplorer const* sse = se->create_subexplorer( 0 ) ;
      PDE_DiscreteField const* ff = dfs->item( sse->string_data( "field" ) ) ;
      CCs.push_back( ff ) ;
      L_CCs.push_back( sse->int_data( "level_of_phase_field" ) ) ;
      P_CCs.push_back( sse->double_data( "parameter_value" ) ) ;
      sse->destroy() ;
      
      ++NB_PHASES ;
   }
   se->destroy() ;
   
   ++NB_PHASES ;
   P_CCs.push_back( exp->double_data( "parameter_value_in_last_phase" ) ) ;
   
   PHs.re_initialize( NB_PHASES ) ;
   dPHs.re_initialize( NB_PHASES ) ;

   string average = exp->string_data( "average" ) ;
   if( average == "arithmetic" )
   {  
      type = 1 ;
   }
   else if( average == "harmonic" )
   {  
      type = 2 ;
   }
   else if( average == "smoothed_Heavyside" )
   {  
      type = 3 ;
   }
   else
   {
      PEL_Error::object()->raise_bad_data_value( exp, "average", 
                                          "    \"harmonic\"\n"
                                          "    \"arithmetic\"\n"
                                          "    \"smoothed_Heavyside\"\n" ) ;
   }
}

//-------------------------------------------------------------------------
CH_MeanParameter:: ~CH_MeanParameter( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
size_t
CH_MeanParameter:: nb_components( void ) const
//-------------------------------------------------------------------------
{
   size_t result = 1 ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
CH_MeanParameter:: cell_value_at_IP( FE_TimeIterator const* t_it,
                                     PDE_LocalFEcell const* fe,
                                     size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "CH_MeanParameter:: cell_value_at_IP" ) ;
   PEL_CHECK_PRE( cell_value_at_IP_PRE( t_it, fe, ic ) ) ;

   size_t idx_max = NB_PHASES-1 ;
   double sum_c = 0.0 ;
   for( size_t i=0 ; i<idx_max ; ++i )
   {
      double xx = fe->value_at_IP( CCs[i], L_CCs[i] ) ;
      PHs( i ) = xx ;
      sum_c += xx ;
   }
   PHs( idx_max ) = 1.0 - sum_c ;

   return( mean_value() ) ;
}

//-------------------------------------------------------------------------
double
CH_MeanParameter:: cell_value_at_pt( FE_TimeIterator const* t_it,
                                     PDE_LocalFEcell const* fe,
                                     size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "CH_MeanParameter:: cell_value_at_pt" ) ;
   PEL_CHECK_PRE( cell_value_at_pt_PRE( t_it, fe, ic ) ) ;

   size_t idx_max = NB_PHASES-1 ;
   double sum_c = 0.0 ;
   for( size_t i=0 ; i<idx_max ; ++i )
   {
      double xx = fe->value_at_pt( CCs[i], L_CCs[i] ) ;
      PHs( i ) = xx ;
      sum_c += xx ;
   }
   PHs( idx_max ) = 1.0 - sum_c ;

   return( mean_value() ) ;
}

//-------------------------------------------------------------------------
double
CH_MeanParameter:: bound_value_at_IP( FE_TimeIterator const* t_it,
                                      PDE_LocalFEbound const* fe,
                                      size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "CH_MeanParameter:: bound_value_at_IP" ) ;
   PEL_CHECK_PRE( bound_value_at_IP_PRE( t_it, fe, ic ) ) ;

   size_t idx_max = NB_PHASES-1 ;
   double sum_c = 0.0 ;
   for( size_t i=0 ; i<idx_max ; ++i )
   {
      double xx = fe->value_at_IP( CCs[i], L_CCs[i] ) ;
      PHs( i ) = xx ;
      sum_c += xx ;
   }
   PHs( idx_max ) = 1.0 - sum_c ;

   return( mean_value() ) ;
}

//-------------------------------------------------------------------------
double
CH_MeanParameter:: cell_gradient_at_IP( FE_TimeIterator const* t_it,
                                        PDE_LocalFEcell const* fe,
                                        size_t a,
                                        size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "CH_MeanParameter:: cell_gradient_at_IP" ) ;
   PEL_CHECK_PRE( cell_gradient_at_IP_PRE( t_it, fe, a, ic ) ) ;
   
   size_t idx_max = NB_PHASES-1 ;
   double sum_xx = 0.0 ;
   double sum_yy = 0.0 ;
   for( size_t i=0 ; i<idx_max ; ++i )
   {
      double xx = fe->value_at_IP( CCs[i], L_CCs[i] ) ;
      PHs( i ) = xx ;
      sum_xx += xx ;
      
      double yy = fe->gradient_at_IP( CCs[i], L_CCs[i], a, 0 ) ;
      dPHs( i ) = yy ;
      sum_yy += yy ;
   }
   PHs( idx_max ) = 1.0 - sum_xx ;
   dPHs( idx_max ) = - sum_yy ;

   return( grad_mean_value() ) ;
}

//-------------------------------------------------------------------------
void
CH_MeanParameter:: prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{ 
   PEL_LABEL( "CH_MeanParameter:: prepare_for_value_on_cells" ) ;
   PEL_CHECK( prepare_for_value_on_cells_PRE( fe ) ) ; 

   for( size_t i=0 ; i<CCs.size() ; ++i ) 
   {
      fe->require_field_calculation( CCs[i], PDE_LocalFE::N )  ;   
   }
}

//-------------------------------------------------------------------------
void
CH_MeanParameter:: prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const 
//-------------------------------------------------------------------------
{ 
   PEL_LABEL( "CH_MeanParameter:: prepare_for_value_on_bounds" ) ;
   PEL_CHECK( prepare_for_value_on_bounds_PRE( fe ) ) ; 

   for( size_t i=0 ; i<CCs.size() ; ++i ) 
   {
      fe->require_field_calculation( CCs[i], PDE_LocalFE::N )  ;   
   }
}

//-------------------------------------------------------------------------
void
CH_MeanParameter:: prepare_for_value_on_sides( PDE_CursorFEside* fe ) const 
//-------------------------------------------------------------------------
{ 
   PEL_LABEL( "CH_MeanParameter:: prepare_for_value_on_sides" ) ;
   PEL_CHECK( prepare_for_value_on_sides_PRE( fe ) ) ; 

   for( size_t i=0 ; i<CCs.size() ; ++i ) 
   {
      fe->require_field_calculation( CCs[i], PDE_LocalFE::N )  ;   
   }
}

//-------------------------------------------------------------------------
void
CH_MeanParameter:: prepare_for_gradient_on_cells( PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{   
   PEL_LABEL( "CH_MeanParameter:: prepare_for_gradient_on_cells" ) ;
   PEL_CHECK( prepare_for_gradient_on_cells_PRE( fe ) ) ;

   for( size_t i=0 ; i<CCs.size() ; ++i ) 
   {
      fe->require_field_calculation( CCs[i], PDE_LocalFE::N )  ;   
      fe->require_field_calculation( CCs[i], PDE_LocalFE::dN )  ;   
   }
}

//-------------------------------------------------------------------------
double
CH_MeanParameter:: mean_value( void )  const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "CH_MeanParameter:: mean_value" ) ;

   double result = 0.0 ;
   
   //????????? qui est eps
   double eps = 1.e-5 ;
  
   if( type == 1 )
   {
      for( size_t i=0 ; i<NB_PHASES ; ++i )
      {
         double c = PHs( i ) ;
         
         double xx = PEL::bad_double() ;
         if( c < eps )
         {
            xx = 0.0 ;
         }
         else if ( c > (1.0 - eps) )
         {
            xx = P_CCs[ i ] ;
         }
         else
         {
            xx = P_CCs[ i ] * c ;
         }
         result += xx ; 
      }
   }
   else if( type == 2 )
   {
      for( size_t i=0 ; i<NB_PHASES ; ++i )
      {
         double c = PHs( i ) ;
         
         double xx = PEL::bad_double() ;
         if( c < eps )
         {
            xx = 0.0 ;
         }
         else if ( c > (1.0-eps) )
         {
            xx = 1.0 / P_CCs[ i ] ;
         }
         else
         {
            xx = c / P_CCs[ i ] ;
         }
         result += xx ;
      }
      result = 1.0 / result ;
   }
   else if( type == 3 )
   {
      double tot_ww = 0. ;
      for( size_t i=0 ; i<NB_PHASES ; ++i )
      {
         //??????????????? qui est eps3
         double eps_3 = 0.5 ;
         double ww = smoothed_Heavyside( PHs( i )-0.5, eps_3 ) ;
         result += P_CCs[ i ] * ww ;
         tot_ww += ww  ;
      }
      result /= tot_ww ;
   }

   return( result ) ;
}

//-------------------------------------------------------------------------
double
CH_MeanParameter:: grad_mean_value( void )  const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "CH_MeanParameter:: grad_mean_value" ) ;
   
   //???
   //??? NOT TESTED
   //???

   double result = 0.0 ;
   //????????? qui est eps
   double eps = 1.e-5 ;

   if( type == 1 )
   {     
      for( size_t i=0 ; i<NB_PHASES ; ++i )
      {
         double c = PHs( i ) ;
         double dc = dPHs( i ) ;
         
         double xx = PEL::bad_double() ;
         if( c < eps )
         {
            xx = 0.0 ;
         }
         else if ( c > (1.0 - eps) )
         {
            xx = P_CCs[ i ] ;
         }
         else
         {
            xx = P_CCs[ i ] * dc ;
         }
         result += xx ; 
      }
   }
   else if( type == 2 )
   {
      PEL_Error::object()->raise_not_implemented( this, "grad_mean_value") ;
   }
   else if( type == 3 )
   {
      double tot_ww   = 0.0 ;
      double d_tot_ww = 0.0 ;
      double xx = 0.0 ;
      double yy = 0.0 ;
      for( size_t i=0 ; i<NB_PHASES ; ++i )
      {
         //??????????????? qui est eps3
         double eps_3 = 0.5 ;
         
         double ww  = smoothed_Heavyside( PHs( i )-0.5, eps_3 )  ;
         tot_ww += ww  ;
         xx += P_CCs[ i ] * ww ;                           
         
         double d_ww = dPHs( i ) * smoothed_delta( PHs( i )-0.5, eps_3 ) ;
         d_tot_ww += d_ww ;
         yy += P_CCs[ i ] * d_ww ;
      }   
      result = ( yy - xx*d_tot_ww/tot_ww ) / tot_ww ; 
   }
   return( result ) ;
}

//---------------------------------------------------------------------------
double
CH_MeanParameter:: smoothed_delta( double x, double eps )
//---------------------------------------------------------------------------
{
   double result = PEL::max_double() ;
   if( PEL::abs( x )< eps )
   {
      result = 0.5*( 1.0 + PEL::cos( PEL::pi() * x / eps ) )/eps ; 
   }
   else 
   {
      result = 0.0 ;
   }
   return( result ) ;
}

//---------------------------------------------------------------------------
double
CH_MeanParameter:: smoothed_Heavyside( double x, double eps )
//---------------------------------------------------------------------------
{
   double result = PEL::max_double() ;
   if( x < -eps )
   {
      result = 0.0 ; 
   }
   else if( x > eps )
   {
      result = 1.0 ;
   }
   else
   {
      result = 0.5 * ( 1. + x/eps + PEL::sin( PEL::pi()*x/eps)/PEL::pi() ) ;
   }       
   return( result ) ;
}
