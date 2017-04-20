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

#include <PDE_InterfaceAndFields.hh>

#include <PEL_Vector.hh>
#include <PEL_assertions.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectReader.hh>
#include <PEL_ObjectWriter.hh>

#include <GE_SetOfPoints.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainBuilder.hh>
#include <PDE_DOFsInitializer.hh>
#include <PDE_InterfaceBuilder.hh>
#include <PDE_LocalFEmortarSide.hh>
#include <PDE_SetOfDiscreteFields.hh>

//-------------------------------------------------------------------------
PDE_InterfaceAndFields*
PDE_InterfaceAndFields:: create( PEL_Object* a_owner,
                                 PEL_ModuleExplorer* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_InterfaceAndFields:: create" ) ;

   PDE_InterfaceAndFields* result = 
                    new PDE_InterfaceAndFields( a_owner, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_InterfaceAndFields:: PDE_InterfaceAndFields( PEL_Object* a_owner,
                                                 PEL_ModuleExplorer* exp )
//-------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , NAME( exp->string_data( "name" ) )
   , DOM0( 0 )
   , DOM1( 0 )
   , BUILDER( 0 )
{
   PDE_DomainBuilder* db0 = PDE_DomainBuilder::object( 
                                exp->string_data( "adjacent_domain_0" ) ) ;
   PDE_DomainBuilder* db1 = PDE_DomainBuilder::object( 
                                exp->string_data( "adjacent_domain_1" ) ) ;
   
   BUILDER = PDE_InterfaceBuilder::create( this, exp, db0, db1 ) ;
   
   PDE_DOFsInitializer::object()->initialize_discrete_fields( this, exp ) ;

   DOM0 = db0->facade() ;
   DOM1 = db1->facade() ;
}

//-------------------------------------------------------------------------
PDE_InterfaceAndFields:: ~PDE_InterfaceAndFields( void )
//-------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
std::string const& 
PDE_InterfaceAndFields:: name( void ) const
//----------------------------------------------------------------------
{
   return( NAME ) ;
}

//----------------------------------------------------------------------
PDE_DomainAndFields const*
PDE_InterfaceAndFields:: adjacent_domain( size_t i_adj ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_InterfaceAndFields:: name_of_adjacent_domain" ) ;
   PEL_CHECK_PRE( i_adj==0 || i_adj==1 ) ;

   PDE_DomainAndFields const* result = ( i_adj==0 ? DOM0 : DOM1 ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_SetOfDiscreteFields*
PDE_InterfaceAndFields:: set_of_discrete_fields( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_InterfaceAndFields:: set_of_discrete_fields" ) ;

   PDE_SetOfDiscreteFields* result = BUILDER->set_of_discrete_fields() ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   PEL_CHECK_POST(
      FORALL( ( result->start() ; result->is_valid() ; result->go_next() ),
              result->item()->owner()==result ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
GE_SetOfPoints const*
PDE_InterfaceAndFields:: set_of_vertices( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_InterfaceAndFields:: set_of_vertices" ) ;
   PEL_CHECK_INV( invariant() ) ;

   GE_SetOfPoints* result = BUILDER->set_of_vertices() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_LocalFEmortarSide*
PDE_InterfaceAndFields:: create_LocalFEmortarSide( PEL_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_InterfaceAndFields:: create_LocalFEmortarSide" ) ;

   PDE_LocalFEmortarSide* result = 
                   new PDE_LocalFEmortarSide( a_owner, 
                                              BUILDER->nb_space_dimensions(), 
                                              BUILDER->mortar_sides() ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   return( result ) ; 
}

//-------------------------------------------------------------------------
void
PDE_InterfaceAndFields:: save_state( PEL_ObjectWriter* writer ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_InterfaceAndFields:: save_state" ) ;
   PEL_CHECK_PRE( save_state_PRE( writer ) ) ;

   writer->start_new_object( "PDE_InterfaceAndFields" ) ;
  
   // Fields storing :
   set_of_discrete_fields()->save_state( writer ) ;

   writer->finalize_object() ;

   PEL_CHECK_POST( save_state_POST( writer ) ) ;
}

//-------------------------------------------------------------------------
void
PDE_InterfaceAndFields:: restore_state( PEL_ObjectReader* reader )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_InterfaceAndFields:: restore_state" ) ;
   PEL_CHECK_PRE( restore_state_PRE( reader ) ) ;

   reader->start_object_retrieval( "PDE_InterfaceAndFields" ) ;
              
   // Retrieving fields :
   set_of_discrete_fields()->restore_state( reader ) ;

   reader->end_object_retrieval() ;

   PEL_CHECK_POST( restore_state_POST( reader ) ) ;
}
