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

#include <FE_SumOfParameters.hh>

#include <FE_SetOfParameters.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

using std::string ;

FE_SumOfParameters const* 
FE_SumOfParameters:: PROTOTYPE = new FE_SumOfParameters() ;

//-------------------------------------------------------------------------
FE_SumOfParameters:: FE_SumOfParameters( void )
//-------------------------------------------------------------------------
   : FE_Parameter( "FE_SumOfParameters" )
{
}

//-------------------------------------------------------------------------
FE_SumOfParameters*
FE_SumOfParameters:: create_replica( PEL_Object* a_owner,
                                     PDE_DomainAndFields const* dom,
                                     PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   FE_SumOfParameters* result = new FE_SumOfParameters( a_owner, dom, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_SumOfParameters:: FE_SumOfParameters( PEL_Object* a_owner,
                                         PDE_DomainAndFields const* dom,
                                         PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : FE_Parameter( a_owner, exp->string_data( "name" ) )
   , NBCS( 0 )
{
   PEL_ModuleExplorer* e = exp->create_subexplorer( 0, "list_of_parameters" ) ;
   e->start_module_iterator() ;
   for( ; e->is_valid_module() ; e->go_next_module() )
   {
      PEL_ModuleExplorer* ee = e->create_subexplorer( e ) ;
      COEFS.push_back( ee->double_data( "coefficient" ) ) ;

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

      if( prm != 0 && NBCS == 0 )
      {
         NBCS = prm->nb_components() ;
      }
      else if( prm != 0 && NBCS != prm->nb_components() )
      {
         std::string mesg =
            "*** FE_SumOfParameters error:\n"
            "    should be defined from FE_Parameter instances\n"
            "    with the same number of components" ;
         PEL_Error::object()->raise_plain( mesg ) ;
      }
      PRMS.push_back( prm ) ;
   }
   e->destroy() ;

   if( PRMS.size() == 0 )
   {
      std::string mesg =
         "*** FE_SumOfParameters error:\n"
         "    empty list_of_parameters" ;
      PEL_Error::object()->raise_plain( mesg ) ;
   }
}

//-------------------------------------------------------------------------
FE_SumOfParameters:: ~FE_SumOfParameters( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
FE_SumOfParameters:: do_the_links( FE_SetOfParameters const* prms )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: do_the_links" ) ;
   PEL_CHECK_PRE( do_the_links_PRE( prms ) ) ;

   for( size_t i=0 ; i<PRMS.size() ; ++i )
   {
      FE_Parameter* prm = PRMS[i] ;
      if( prm == 0 )
      {
         prm = prms->item( NAMES[i] ) ;
         PRMS[i] = prm ;
         if( NBCS == 0 )
         {
            NBCS = prm->nb_components() ;
         }
         else if( NBCS != prm->nb_components() )
         {
            std::string mesg = "FE_SumOfParameters should be defined from\n" ;
            mesg += "FE_Parameter instances with the same number of components" ;
            PEL_Error::object()->raise_plain( mesg ) ;
         }
      }
      prm->do_the_links( prms ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_SumOfParameters:: reset( PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: reset" ) ;

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
FE_SumOfParameters:: nb_components( void ) const
//-------------------------------------------------------------------------
{
   return( NBCS ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: cell_value( FE_TimeIterator const* t_it,
                                 PDE_LocalFEcell const* fe,
                                 size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: cell_value" ) ;
   PEL_CHECK_PRE( cell_value_PRE( t_it, fe, ic ) ) ;

   double result = 0.0 ;
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->cell_value( t_it, fe, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: cell_gradient( FE_TimeIterator const* t_it,
                                    PDE_LocalFEcell const* fe,
                                    size_t a,
                                    size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: cell_gradient" ) ;
   PEL_CHECK_PRE( cell_gradient_PRE( t_it, fe, a, ic ) ) ;

   double result = 0.0 ;
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->cell_gradient( t_it, fe, a, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: cell_value_at_pt( FE_TimeIterator const* t_it,
                                       PDE_LocalFEcell const* fe,
                                       size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: cell_value_at_pt" ) ;
   PEL_CHECK_PRE( cell_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double result = 0.0 ;
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->cell_value_at_pt( t_it, fe, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: cell_gradient_at_pt( FE_TimeIterator const* t_it,
                                          PDE_LocalFEcell const* fe,
                                          size_t a,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: cell_gradient_at_pt" ) ;
   PEL_CHECK_PRE( cell_gradient_at_pt_PRE( t_it, fe, a, ic ) ) ;

   double result = 0.0 ;
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->cell_gradient_at_pt( t_it, fe, a, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: cell_value_at_IP( FE_TimeIterator const* t_it,
                                       PDE_LocalFEcell const* fe,
                                       size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: cell_value_at_IP" ) ;
   PEL_CHECK_PRE( cell_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double result = 0.0 ;
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->cell_value_at_IP( t_it, fe, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: cell_gradient_at_IP( FE_TimeIterator const* t_it,
                                          PDE_LocalFEcell const* fe,
                                          size_t a,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: cell_gradient_at_IP" ) ;
   PEL_CHECK_PRE( cell_gradient_at_IP_PRE( t_it, fe, a, ic ) ) ;

   double result = 0.0 ;
   //???? gain de temps
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->cell_gradient_at_IP( t_it, fe, a, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: bound_value( FE_TimeIterator const* t_it,
                                  PDE_LocalFEbound const* fe,
                                  size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: bound_value" ) ;
   PEL_CHECK_PRE( bound_value_PRE( t_it, fe, ic ) ) ;

   double result = 0.0 ;
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->bound_value( t_it, fe, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: bound_gradient( FE_TimeIterator const* t_it,
                                     PDE_LocalFEbound const* fe,
                                     size_t a,
                                     size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: bound_gradient" ) ;
   PEL_CHECK_PRE( bound_gradient_PRE( t_it, fe, a, ic ) ) ;

   double result = 0.0 ;
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->bound_gradient( t_it, fe, a, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: bound_value_at_pt( FE_TimeIterator const* t_it,
                                        PDE_LocalFEbound const* fe,
                                        size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: bound_value_at_pt" ) ;
   PEL_CHECK_PRE( bound_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double result = 0.0 ;
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->bound_value_at_pt( t_it, fe, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: bound_gradient_at_pt( FE_TimeIterator const* t_it,
                                           PDE_LocalFEbound const* fe,
                                           size_t a,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: bound_gradient_at_pt" ) ;
   PEL_CHECK_PRE( bound_gradient_at_pt_PRE( t_it, fe, a, ic ) ) ;

   double result = 0.0 ;
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->bound_gradient_at_pt( t_it, fe, a, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: bound_value_at_IP( FE_TimeIterator const* t_it,
                                        PDE_LocalFEbound const* fe,
                                        size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: bound_value_at_IP" ) ;
   PEL_CHECK_PRE( bound_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double result = 0.0 ;
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->bound_value_at_IP( t_it, fe, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: bound_gradient_at_IP( FE_TimeIterator const* t_it,
                                           PDE_LocalFEbound const* fe,
                                           size_t a,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: bound_gradient_at_IP" ) ;
   PEL_CHECK_PRE( bound_gradient_at_IP_PRE( t_it, fe, a, ic ) ) ;

   double result = 0.0 ;
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->bound_gradient_at_IP( t_it, fe, a, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: side_value( FE_TimeIterator const* t_it,
                                 PDE_CursorFEside const* fe,
                                 size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: side_value" ) ;
   PEL_CHECK_PRE( side_value_PRE( t_it, fe, ic ) ) ;

   double result = 0.0 ;
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->side_value( t_it, fe, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: side_gradient( FE_TimeIterator const* t_it,
                                    PDE_CursorFEside const* fe,
                                    size_t a,
                                    size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: side_gradient" ) ;
   PEL_CHECK_PRE( side_gradient_PRE( t_it, fe, a, ic ) ) ;

   double result = 0.0 ;
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->side_gradient( t_it, fe, a, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: side_value_at_pt( FE_TimeIterator const* t_it,
                                       PDE_CursorFEside const* fe,
                                       size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: side_value_at_pt" ) ;
   PEL_CHECK_PRE( side_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double result = 0.0 ;
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->side_value_at_pt( t_it, fe, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: side_gradient_at_pt( FE_TimeIterator const* t_it,
                                          PDE_CursorFEside const* fe,
                                          size_t a,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: side_gradient_at_pt" ) ;
   PEL_CHECK_PRE( side_gradient_at_pt_PRE( t_it, fe, a, ic ) ) ;

   double result = 0.0 ;
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->side_gradient_at_pt( t_it, fe, a, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: side_value_at_IP( FE_TimeIterator const* t_it,
                                       PDE_CursorFEside const* fe,
                                       size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: side_value_at_IP" ) ;
   PEL_CHECK_PRE( side_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double result = 0.0 ;
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->side_value_at_IP( t_it, fe, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_SumOfParameters:: side_gradient_at_IP( FE_TimeIterator const* t_it,
                                          PDE_CursorFEside  const* fe,
                                          size_t a,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: side_gradient_at_IP" ) ;
   PEL_CHECK_PRE( side_gradient_at_IP_PRE( t_it, fe, a, ic ) ) ;

   double result = 0.0 ;
   for( size_t i=0 ; i<COEFS.size() ; ++i )
   {
      result += COEFS[i] * PRMS[i]->side_gradient_at_IP( t_it, fe, a, ic ) ;
   }
   return( result ) ;
}

//-------------------------------------------------------------------------
void
FE_SumOfParameters:: prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: prepare_for_value_on_cells" ) ;
   PEL_CHECK( prepare_for_value_on_cells_PRE( fe ) ) ;

   for( size_t i=0 ; i<PRMS.size() ; ++i )
   {
      PRMS[i]->transfer_cell_calculation_requirements( fe, 
                                                       FE_Parameter::Val ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_SumOfParameters:: prepare_for_gradient_on_cells( PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: prepare_for_gradient_on_cells" ) ;
   PEL_CHECK( prepare_for_gradient_on_cells_PRE( fe ) ) ;

   for( size_t i=0 ; i<PRMS.size() ; ++i )
   {
      PRMS[i]->transfer_cell_calculation_requirements( fe, 
                                                       FE_Parameter::Grad ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_SumOfParameters:: prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: prepare_for_value_on_bounds" ) ;
   PEL_CHECK( prepare_for_value_on_bounds_PRE( fe ) ) ;

   for( size_t i=0 ; i<PRMS.size() ; ++i )
   {
      PRMS[i]->transfer_bound_calculation_requirements( fe, 
                                                        FE_Parameter::Val ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_SumOfParameters:: prepare_for_gradient_on_bounds( PDE_LocalFEbound* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: prepare_for_gradient_on_bounds" ) ;
   PEL_CHECK( prepare_for_gradient_on_bounds_PRE( fe ) ) ;

   for( size_t i=0 ; i<PRMS.size() ; ++i )
   {
      PRMS[i]->transfer_bound_calculation_requirements( fe, 
                                                        FE_Parameter::Grad ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_SumOfParameters:: prepare_for_value_on_sides( PDE_CursorFEside* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: prepare_for_value_on_sides" ) ;
   PEL_CHECK( prepare_for_value_on_sides_PRE( fe ) ) ;

   for( size_t i=0 ; i<PRMS.size() ; ++i )
   {
      PRMS[i]->transfer_side_calculation_requirements( fe, 
                                                       FE_Parameter::Val ) ;
   }
}

//-------------------------------------------------------------------------
void
FE_SumOfParameters:: prepare_for_gradient_on_sides( PDE_CursorFEside* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SumOfParameters:: prepare_for_gradient_on_sides" ) ;
   PEL_CHECK( prepare_for_gradient_on_sides_PRE( fe ) ) ;

   for( size_t i=0 ; i<PRMS.size() ; ++i )
   {
      PRMS[i]->transfer_side_calculation_requirements( fe, 
                                                       FE_Parameter::Grad ) ;
   }
}
