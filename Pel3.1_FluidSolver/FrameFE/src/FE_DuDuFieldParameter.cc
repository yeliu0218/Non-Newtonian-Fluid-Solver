/*
 *  Copyright : 
 *    "Institut de Radioprotection et de S�ret� Nucl�aire - IRSN" (1995-2008)
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

#include <FE.hh>
#include <FE_DuDuFieldParameter.hh>
#include <FE_SetOfParameters.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include<cmath>

FE_DuDuFieldParameter const* 
FE_DuDuFieldParameter:: PROTOTYPE = new FE_DuDuFieldParameter() ;

//-------------------------------------------------------------------------
FE_DuDuFieldParameter:: FE_DuDuFieldParameter( void )
//-------------------------------------------------------------------------
   : FE_Parameter( "FE_DuDuFieldParameter" )
   , coef(0.5)
   , power(0.5)
{
}

//-------------------------------------------------------------------------
FE_DuDuFieldParameter*
FE_DuDuFieldParameter:: create_replica( PEL_Object* a_owner,
        	                            PDE_DomainAndFields const* dom,
                                        PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DuDuFieldParameter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   FE_DuDuFieldParameter* result = 
                              new FE_DuDuFieldParameter( a_owner, dom, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_DuDuFieldParameter:: FE_DuDuFieldParameter( PEL_Object* a_owner,
        		                               PDE_DomainAndFields const* dom,
                                               PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : FE_Parameter( a_owner, exp->string_data( "name" ) )
   , FF( dom->set_of_discrete_fields()->item( exp->string_data( "field_name" ) ) )
   , L_FF( exp->int_data( "field_level" ) )
   , NBCS( 1 )
   , coef( exp->double_data( "coef" ) )
   , power(exp->double_data( "power" )  )
{
   if( FF->nb_components() < 1 )
   {
      std::string mesg =
         "*** FE_DuDuFieldParameter : the associated field\n"
         "    should have at least one component" ;
      PEL_Error::object()->raise_plain( mesg ) ;
   }
}

//-------------------------------------------------------------------------
FE_DuDuFieldParameter:: ~FE_DuDuFieldParameter( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
size_t
FE_DuDuFieldParameter:: nb_components( void ) const
//-------------------------------------------------------------------------
{
   return( NBCS ) ;
}

//-------------------------------------------------------------------------
double
FE_DuDuFieldParameter:: cell_value_at_pt( FE_TimeIterator const* t_it,
                                          PDE_LocalFEcell const* fe,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DuDuFieldParameter:: cell_value_at_pt" ) ;
   PEL_CHECK_PRE( cell_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double result = 0;
   FE::add_D_u_D_u_at_pt( result, FF, L_FF, fe, coef );

   return( pow(result, power) ) ;
}

//-------------------------------------------------------------------------
double
FE_DuDuFieldParameter:: cell_value_at_IP( FE_TimeIterator const* t_it,
                                          PDE_LocalFEcell const* fe,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DuDuFieldParameter:: cell_value_at_IP" ) ;
   PEL_CHECK_PRE( cell_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double result = 0;
   FE::add_D_u_D_u_at_IP( result, FF, L_FF, fe, coef );

   return( pow(result, power) ) ;
}

//-------------------------------------------------------------------------
double
FE_DuDuFieldParameter:: bound_value_at_pt( FE_TimeIterator const* t_it,
                                           PDE_LocalFEbound const* fe,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DuDuFieldParameter:: bound_value_at_pt" ) ;
   PEL_CHECK_PRE( bound_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double result = 0;
   FE::add_D_u_D_u_at_pt( result, FF, L_FF, fe, coef );

   return( pow(result, power) ) ;
}

//-------------------------------------------------------------------------
double
FE_DuDuFieldParameter:: bound_value_at_IP( FE_TimeIterator const* t_it,
                                           PDE_LocalFEbound const* fe,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DuDuFieldParameter:: bound_value_at_IP" ) ;
   PEL_CHECK_PRE( bound_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double result = 0;
   FE::add_D_u_D_u_at_IP( result, FF, L_FF, fe, coef );

   return( pow(result, power) ) ;
}

//-------------------------------------------------------------------------
void
FE_DuDuFieldParameter:: prepare_for_value_on_cells( PDE_LocalFEcell* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DuDuFieldParameter:: prepare_for_value_on_cells" ) ;
   PEL_CHECK( prepare_for_value_on_cells_PRE( fe ) ) ;
   
   fe->require_field_calculation( FF, PDE_LocalFE::dN )  ;
}

//-------------------------------------------------------------------------
void
FE_DuDuFieldParameter:: prepare_for_value_on_bounds( PDE_LocalFEbound* fe ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_DuDuFieldParameter:: prepare_for_value_on_bounds" ) ;
   PEL_CHECK( prepare_for_value_on_bounds_PRE( fe ) ) ;

   fe->require_field_calculation( FF, PDE_LocalFE::dN )  ;
}
