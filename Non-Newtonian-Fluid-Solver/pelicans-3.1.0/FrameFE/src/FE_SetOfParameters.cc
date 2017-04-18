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

#include <FE_SetOfParameters.hh>

#include <PEL_DataWithContext.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <PDE_DomainAndFields.hh>
#include <PDE_SetOfDomains.hh>

#include <FE_Parameter.hh>
#include <FE_UniformParameter.hh>

#include <iostream>
#include <sstream>

using std::string ;
using std::ostringstream ;
using std::endl ;

//----------------------------------------------------------------------
FE_SetOfParameters*
FE_SetOfParameters:: create( PEL_Object* a_owner,
                             PDE_DomainAndFields const* dom,
                             PEL_ModuleExplorer* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_SetOfParameters:: create" ) ;

   FE_SetOfParameters* result = new FE_SetOfParameters( a_owner, dom, exp ) ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_SetOfParameters:: FE_SetOfParameters( PEL_Object* a_owner,
                                         PDE_DomainAndFields const* dom,
                                         PEL_ModuleExplorer* exp )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
{
   PEL_LABEL( "FE_SetOfParameters:: FE_SetOfParameters" ) ;

   build_uniform_parameters( exp ) ;

   exp->start_module_iterator() ;
   for( ; exp->is_valid_module() ; exp->go_next_module() )
   {
      PEL_ModuleExplorer* ee = exp->create_subexplorer( 0 ) ; 
      FE_Parameter* prm = FE_Parameter::make( this, dom, ee ) ;
      add_one_parameter( prm ) ;
      ee->destroy() ;
   }

   std::map<std::string,FE_Parameter*>::const_iterator it = MAP.begin() ;
   for( ; it != MAP.end() ; ++it )
   {
      it->second->do_the_links( this ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
FE_SetOfParameters*
FE_SetOfParameters:: create( PEL_Object* a_owner,
                             PDE_SetOfDomains const* sdoms,
                             PEL_ModuleExplorer* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_SetOfParameters:: create" ) ;

   FE_SetOfParameters* result = new FE_SetOfParameters( a_owner, sdoms, exp ) ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_SetOfParameters:: FE_SetOfParameters( PEL_Object* a_owner,
                                         PDE_SetOfDomains const* sdoms,
                                         PEL_ModuleExplorer* exp )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
{
   PEL_LABEL( "FE_SetOfParameters:: FE_SetOfParameters" ) ;

   build_uniform_parameters( exp ) ;

   exp->start_module_iterator() ;
   for( ; exp->is_valid_module() ; exp->go_next_module() )
   {
      PEL_ModuleExplorer* ee = exp->create_subexplorer( 0 ) ; 
      FE_Parameter* prm = FE_Parameter::make( this, sdoms, ee ) ;
      add_one_parameter( prm ) ;
      ee->destroy() ;
   }

   std::map<std::string,FE_Parameter*>::const_iterator it = MAP.begin() ;
   for( ; it != MAP.end() ; ++it )
   {
      it->second->do_the_links( this ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
FE_SetOfParameters:: ~FE_SetOfParameters( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_SetOfParameters:: ~FE_SetOfParameters" ) ;
}

//----------------------------------------------------------------------
size_t
FE_SetOfParameters:: nb_parameters( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_SetOfParameters:: nb_parameters" ) ;

   return( MAP.size() ) ;
}

//----------------------------------------------------------------------
bool
FE_SetOfParameters:: has( std::string const& name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_SetOfParameters:: has" ) ;

   bool result = ( MAP.count( name ) != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_Parameter*
FE_SetOfParameters:: item( std::string const& name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_SetOfParameters:: item(variable_name)" ) ;
   
   FE_Parameter* result = 0 ;

   std::map<std::string,FE_Parameter*>::const_iterator it = MAP.find( name ) ;
   if( it != MAP.end() )
   {
      result = (*it).second ;
   }
   else
   {
      ostringstream mesg ;
      mesg << "FE_SetOfParameters :" << endl ;
      mesg << "   request for an unknown FE_Parameter object" << endl ;
      mesg << "   of name : \"" << name << "\"" << endl ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   PEL_CHECK_POST( result->name() == name ) ;
   return( result ) ;
}

//-----------------------------------------------------------------------
void
FE_SetOfParameters:: build_uniform_parameters( PEL_ModuleExplorer* exp )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "FE_SetOfParameters:: build_uniform_parameters" ) ;

   exp->start_entry_iterator() ;
   for( ; exp->is_valid_entry() ; exp->go_next_entry() )
   {
      PEL_DataWithContext const* data = exp->data( 0 ) ;
      if( !data->value_can_be_evaluated() )
         PEL_Error::object()->raise_not_evaluable(
                         exp, exp->keyword(), data->undefined_variables() ) ;
      if( data->data_type() != PEL_Data::DoubleVector )
         PEL_Error::object()->raise_bad_data_type( exp, 
                                                   exp->keyword(), 
                                                   PEL_Data::DoubleVector ) ;
      FE_Parameter* prm = FE_UniformParameter::create( this, 
                                                   exp->keyword(),
						   data->to_double_vector() ) ;
      add_one_parameter( prm ) ;
      data->destroy() ;
   }
}

//-----------------------------------------------------------------------
void
FE_SetOfParameters:: add_one_parameter( FE_Parameter* prm )
//-----------------------------------------------------------------------
{
   PEL_LABEL( "FE_SetOfParameters:: add_one_parameter" ) ;

   std::map<std::string,FE_Parameter*>::const_iterator it = 
                                                MAP.find( prm->name() ) ;
   if( it == MAP.end() )
   {
      MAP[ prm->name() ] = prm ;
   }
   else
   {
      ostringstream mesg ;
      mesg << "FE_SetOfParameters :" << endl ;
      mesg << "   there should be only one FE_Parameter object" << endl ;
      mesg << "   of name : \"" << prm->name() << "\"" << endl ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }
}
