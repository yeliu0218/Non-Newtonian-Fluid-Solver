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

#include <FE_ProductOfParameters.hh>

#include <FE_SetOfParameters.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <iostream>

using std::string ;

FE_ProductOfParameters const* 
FE_ProductOfParameters:: PROTOTYPE = new FE_ProductOfParameters() ;

//-------------------------------------------------------------------------
FE_ProductOfParameters:: FE_ProductOfParameters( void )
//-------------------------------------------------------------------------
   : FE_Parameter( "FE_ProductOfParameters" )
{
}

//-------------------------------------------------------------------------
FE_ProductOfParameters*
FE_ProductOfParameters:: create_replica( PEL_Object* a_owner,
                                         PDE_DomainAndFields const* dom,
                                         PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ProductOfParameters:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   FE_ProductOfParameters* result = new FE_ProductOfParameters( a_owner, dom, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_ProductOfParameters:: FE_ProductOfParameters( PEL_Object* a_owner,
                                                 PDE_DomainAndFields const* dom,
                                                 PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : FE_Parameter( a_owner, exp->string_data( "name" ) )
   , NBCS( PEL::bad_index() )
{
   PEL_ModuleExplorer* e = exp->create_subexplorer( 0, "list_of_parameters" ) ;
   e->start_module_iterator() ;
   for( ; e->is_valid_module() ; e->go_next_module() )
   {
      PEL_ModuleExplorer* ee = e->create_subexplorer( e ) ;

      FE_Parameter* prm = 0 ;
      string const& pmtype = ee->string_data( "type" ) ;
      if( pmtype == "already_defined" )
      {
         NAMES.push_back( ee->string_data( "name" ) ) ;
      }
      else if( pmtype == "to_be_defined" )
      {
         PEL_ModuleExplorer* eee = ee->create_subexplorer( ee, 
                                                           "FE_Parameter" ) ;
         prm = FE_Parameter::make( this, dom, eee ) ;
         NAMES.push_back( prm->name() ) ;
      }
      else
      {
         PEL_Error::object()->raise_bad_data_value( ee, "type",
                             "\"already_defined\" and \"to_be_defined\"" ) ;
      }

      if( prm != 0 ) infer_nb_components( prm ) ;

      PRMS.push_back( prm ) ;
   }
   e->destroy() ;

   if( PRMS.size() == 0 )
   {
      std::string mesg =
         "*** FE_ProductOfParameters error:\n"
         "    empty list_of_parameters" ;
      PEL_Error::object()->raise_plain( mesg ) ;
   }
}

//-------------------------------------------------------------------------
FE_ProductOfParameters:: ~FE_ProductOfParameters( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
FE_ProductOfParameters:: do_the_links( FE_SetOfParameters const* prms )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ProductOfParameters:: do_the_links" ) ;
   PEL_CHECK_PRE( do_the_links_PRE( prms ) ) ;

   for( size_t i=0 ; i<PRMS.size() ; ++i )
   {
      FE_Parameter* prm = PRMS[i] ;
      if( prm == 0 )
      {
         prm = prms->item( NAMES[i] ) ;
         PRMS[i] = prm ;
         infer_nb_components( prm ) ;
      }
      prm->do_the_links( prms ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_ProductOfParameters:: reset( PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ProductOfParameters:: reset" ) ;

   for( size_t i=0 ; i<PRMS.size() ; ++i )
   {
      PEL_ModuleExplorer* e = 0 ;
      if( exp != 0 )
      {
         e = exp->create_subexplorer( 0, NAMES[i] ) ;
      }
      PRMS[i]->reset( e ) ;
      if( e != 0 )
      {
         e->destroy() ; e = 0 ;
      }
   }
}

//-------------------------------------------------------------------------
size_t
FE_ProductOfParameters:: nb_components( void ) const
//-------------------------------------------------------------------------
{
   return( NBCS ) ;
}

//-------------------------------------------------------------------------
double
FE_ProductOfParameters:: cell_value( FE_TimeIterator const* t_it,
                                     PDE_LocalFEcell const* fe,
                                     size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ProductOfParameters:: cell_value" ) ;
   PEL_CHECK_PRE( cell_value_PRE( t_it, fe, ic ) ) ;

   PEL::out() << "dans cell_value" << std::endl ;
   size_t iic = ( PRMS[0]->nb_components() == 1 ? 0 : ic ) ;
   PEL::out() << "iic=" << iic << std::endl ;
   double result = PRMS[0]->cell_value( t_it, fe, iic ) ;
   PEL::out() << "result0=" << result << std::endl ;
   for( size_t i=1 ; i<PRMS.size() ; ++i )
   {
      iic = ( PRMS[i]->nb_components() == 1 ? 0 : ic ) ;
      result *= PRMS[i]->cell_value( t_it, fe, iic ) ;
      PEL::out() << "result" << i << "=" << result << std::endl ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ProductOfParameters:: cell_value_at_pt( FE_TimeIterator const* t_it,
                                           PDE_LocalFEcell const* fe,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ProductOfParameters:: cell_value_at_pt" ) ;
   PEL_CHECK_PRE( cell_value_at_pt_PRE( t_it, fe, ic ) ) ;

   size_t iic = ( PRMS[0]->nb_components() == 1 ? 0 : ic ) ;
   double result = PRMS[0]->cell_value_at_pt( t_it, fe, iic ) ;
   for( size_t i=1 ; i<PRMS.size() ; ++i )
   {
      iic = ( PRMS[i]->nb_components() == 1 ? 0 : ic ) ;
      result *= PRMS[i]->cell_value_at_pt( t_it, fe, iic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ProductOfParameters:: cell_value_at_IP( FE_TimeIterator const* t_it,
                                           PDE_LocalFEcell const* fe,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ProductOfParameters:: cell_value_at_IP" ) ;
   PEL_CHECK_PRE( cell_value_at_IP_PRE( t_it, fe, ic ) ) ;

   size_t iic = ( PRMS[0]->nb_components() == 1 ? 0 : ic ) ;
   double result = PRMS[0]->cell_value_at_IP( t_it, fe, iic ) ;
   for( size_t i=1 ; i<PRMS.size() ; ++i )
   {
      iic = ( PRMS[i]->nb_components() == 1 ? 0 : ic ) ;
      result *= PRMS[i]->cell_value_at_IP( t_it, fe, iic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ProductOfParameters:: bound_value( FE_TimeIterator const* t_it,
                                      PDE_LocalFEbound const* fe,
                                      size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ProductOfParameters:: bound_value" ) ;
   PEL_CHECK_PRE( bound_value_PRE( t_it, fe, ic ) ) ;

   size_t iic = ( PRMS[0]->nb_components() == 1 ? 0 : ic ) ;
   double result = PRMS[0]->bound_value( t_it, fe, iic ) ;
   for( size_t i=1 ; i<PRMS.size() ; ++i )
   {
      iic = ( PRMS[i]->nb_components() == 1 ? 0 : ic ) ;
      result *= PRMS[i]->bound_value( t_it, fe, iic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ProductOfParameters:: bound_value_at_pt( FE_TimeIterator const* t_it,
                                            PDE_LocalFEbound const* fe,
                                            size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ProductOfParameters:: bound_value_at_pt" ) ;
   PEL_CHECK_PRE( bound_value_at_pt_PRE( t_it, fe, ic ) ) ;

   size_t iic = ( PRMS[0]->nb_components() == 1 ? 0 : ic ) ;
   double result = PRMS[0]->bound_value_at_pt( t_it, fe, iic ) ;
   for( size_t i=1 ; i<PRMS.size() ; ++i )
   {
      iic = ( PRMS[i]->nb_components() == 1 ? 0 : ic ) ;
      result *= PRMS[i]->bound_value_at_pt( t_it, fe, iic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ProductOfParameters:: bound_value_at_IP( FE_TimeIterator const* t_it,
                                            PDE_LocalFEbound const* fe,
                                            size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ProductOfParameters:: bound_value_at_IP" ) ;
   PEL_CHECK_PRE( bound_value_at_IP_PRE( t_it, fe, ic ) ) ;

   size_t iic = ( PRMS[0]->nb_components() == 1 ? 0 : ic ) ;
   double result = PRMS[0]->bound_value_at_IP( t_it, fe, iic ) ;
   for( size_t i=1 ; i<PRMS.size() ; ++i )
   {
      iic = ( PRMS[i]->nb_components() == 1 ? 0 : ic ) ;
      result *= PRMS[i]->bound_value_at_IP( t_it, fe, iic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ProductOfParameters:: side_value( FE_TimeIterator const* t_it,
                                     PDE_CursorFEside const* fe,
                                     size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ProductOfParameters:: side_value" ) ;
   PEL_CHECK_PRE( side_value_PRE( t_it, fe, ic ) ) ;

   size_t iic = ( PRMS[0]->nb_components() == 1 ? 0 : ic ) ;
   double result = PRMS[0]->side_value( t_it, fe, iic ) ;
   for( size_t i=1 ; i<PRMS.size() ; ++i )
   {
      iic = ( PRMS[i]->nb_components() == 1 ? 0 : ic ) ;
      result *= PRMS[i]->side_value( t_it, fe, iic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ProductOfParameters:: side_value_at_pt( FE_TimeIterator const* t_it,
                                           PDE_CursorFEside const* fe,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ProductOfParameters:: side_value_at_pt" ) ;
   PEL_CHECK_PRE( side_value_at_pt_PRE( t_it, fe, ic ) ) ;

   size_t iic = ( PRMS[0]->nb_components() == 1 ? 0 : ic ) ;
   double result = PRMS[0]->side_value_at_pt( t_it, fe, iic ) ;
   for( size_t i=1 ; i<PRMS.size() ; ++i )
   {
      iic = ( PRMS[i]->nb_components() == 1 ? 0 : ic ) ;
      result *= PRMS[i]->side_value_at_pt( t_it, fe, iic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_ProductOfParameters:: side_value_at_IP( FE_TimeIterator const* t_it,
                                           PDE_CursorFEside const* fe,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ProductOfParameters:: side_value_at_IP" ) ;
   PEL_CHECK_PRE( side_value_at_IP_PRE( t_it, fe, ic ) ) ;

   size_t iic = ( PRMS[0]->nb_components() == 1 ? 0 : ic ) ;
   double result = PRMS[0]->side_value_at_IP( t_it, fe, iic ) ;
   for( size_t i=1 ; i<PRMS.size() ; ++i )
   {
      iic = ( PRMS[i]->nb_components() == 1 ? 0 : ic ) ;
      result *= PRMS[i]->side_value_at_IP( t_it, fe, iic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
void
FE_ProductOfParameters:: prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ProductOfParameters:: prepare_for_value_on_cells" ) ;
   PEL_CHECK( prepare_for_value_on_cells_PRE( fe ) ) ;

   for( size_t i=0 ; i<PRMS.size() ; ++i )
   {
      PRMS[i]->transfer_cell_calculation_requirements( fe, 
                                                       FE_Parameter::Val ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_ProductOfParameters:: prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ProductOfParameters:: prepare_for_value_on_bounds" ) ;
   PEL_CHECK( prepare_for_value_on_bounds_PRE( fe ) ) ;

   for( size_t i=0 ; i<PRMS.size() ; ++i )
   {
      PRMS[i]->transfer_bound_calculation_requirements( fe, 
                                                        FE_Parameter::Val ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_ProductOfParameters:: prepare_for_value_on_sides( PDE_CursorFEside* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ProductOfParameters:: prepare_for_value_on_sides" ) ;
   PEL_CHECK( prepare_for_value_on_sides_PRE( fe ) ) ;

   for( size_t i=0 ; i<PRMS.size() ; ++i )
   {
      PRMS[i]->transfer_side_calculation_requirements( fe, 
                                                       FE_Parameter::Val ) ;
   }
}

//--------------------------------------------------------------------------
void
FE_ProductOfParameters:: infer_nb_components( FE_Parameter const* prm )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "FE_ProductOfParameters:: infer_nb_components" ) ;
   
   if( NBCS == PEL::bad_index() )
   {
      NBCS = prm->nb_components() ;
   }
   else
   {
      if( prm->nb_components() != 1 )
      {
         if( NBCS != 1 )
         {
            std::string mesg =
               "*** FE_ProductOfParameters error:\n"
               "    should be defined from FE_Parameter instances\n"
               "    such that at most one of them have more than\n"
               "    one component" ;
            PEL_Error::object()->raise_plain( mesg ) ;

         }
         else
         {
            NBCS = prm->nb_components() ;
         }
      }
   }
}

