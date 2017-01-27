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

#include <PDE_PlainFieldProjector.hh>

#include <PEL_assertions.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_LocalFEcell.hh>

//----------------------------------------------------------------------
PDE_PlainFieldProjector*
PDE_PlainFieldProjector:: create(
                      PEL_Object* a_owner,
                      PDE_DiscreteField* a_field,
                      size_t a_field_level,
                      PDE_DiscreteField const* a_source_field,
                      size_t a_source_field_level,
                      PDE_DomainAndFields const* a_dom,
                      PEL_ModuleExplorer* a_exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_PlainFieldProjector:: create" ) ;
   PEL_CHECK_PRE( a_field!=0 ) ;
   PEL_CHECK_PRE( a_field_level<a_field->storage_depth() ) ;
   PEL_CHECK_PRE( a_source_field!=0 ) ;
   PEL_CHECK_PRE( a_field->nb_components()==a_source_field->nb_components() ) ;
   PEL_CHECK_PRE( a_source_field_level<a_source_field->storage_depth() ) ;
   PEL_CHECK_PRE( a_dom!=0 ) ;
   PEL_CHECK_PRE( a_exp!=0 ) ;

   PDE_PlainFieldProjector* result =
      new PDE_PlainFieldProjector( a_owner,
                                   a_field, a_field_level,
                                   a_source_field, a_source_field_level,
                                   a_dom,
                                   a_exp ) ;

   PEL_CHECK_POST( result!=0 && result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->field()==a_field ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
PDE_PlainFieldProjector:: PDE_PlainFieldProjector(
                      PEL_Object* a_owner,
                      PDE_DiscreteField* a_field,
                      size_t a_field_level,
                      PDE_DiscreteField const* a_source_field,
                      size_t a_source_field_level,
                      PDE_DomainAndFields const* a_dom,
                      PEL_ModuleExplorer* a_exp )
//----------------------------------------------------------------------
   : PDE_ProjectorForDOFsSetting( a_owner, a_field, a_field_level, a_dom, a_exp )
   , SRC_FIELD( a_source_field )
   , SRC_LEVEL( a_source_field_level )
{
   PEL_LABEL( "PDE_PlainFieldProjector:: PDE_PlainFieldProjector" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   add_field_requirement_on_cells( SRC_FIELD, PDE_LocalFE::N ) ;
}

//----------------------------------------------------------------------
PDE_PlainFieldProjector:: ~PDE_PlainFieldProjector( void )
//----------------------------------------------------------------------

{
   PEL_LABEL( "PDE_PlainFieldProjector:: ~PDE_PlainFieldProjector" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
PDE_PlainFieldProjector:: compute_value_at_IP( PDE_LocalFEcell const* fe,
                                               doubleVector& result ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_PlainFieldProjector:: compute_value_at_IP" ) ;
   PEL_CHECK( compute_value_at_IP_PRE( fe, result ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   for( size_t i=0 ; i<SRC_FIELD->nb_components() ; ++i )
   {
      result(i) = fe->value_at_IP( SRC_FIELD, SRC_LEVEL, i ) ;
   }
}

//----------------------------------------------------------------------
bool
PDE_PlainFieldProjector:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PDE_ProjectorForDOFsSetting::invariant() ) ;
   return( true ) ;
}
