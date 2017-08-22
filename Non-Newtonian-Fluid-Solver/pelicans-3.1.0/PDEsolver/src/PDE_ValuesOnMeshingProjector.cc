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

#include <PDE_ValuesOnMeshingProjector.hh>

#include <GE_Color.hh>
#include <GE_Point.hh>
#include <GE_Mpolyhedron.hh>

#include <PEL_assertions.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ValuesOnMeshing.hh>

//----------------------------------------------------------------------
PDE_ValuesOnMeshingProjector*
PDE_ValuesOnMeshingProjector:: create(
                      PEL_Object* a_owner,
                      PDE_DiscreteField* a_field,
                      size_t a_field_level,
                      PDE_ValuesOnMeshing* a_field_initializer,
                      PDE_DomainAndFields const* a_dom,
                      PEL_ModuleExplorer* a_exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ValuesOnMeshingProjector:: create" ) ;
   PEL_CHECK_PRE( a_field!=0 ) ;
   PEL_CHECK_PRE( a_field_level<a_field->storage_depth() ) ;
   PEL_CHECK_PRE( a_field_initializer!=0 ) ;
   PEL_CHECK_PRE( a_dom!=0 ) ;
   PEL_CHECK_PRE( a_exp!=0 ) ;

   PDE_ValuesOnMeshingProjector* result =
      new PDE_ValuesOnMeshingProjector( a_owner,
                                               a_field, a_field_level,
                                               a_field_initializer,
                                               a_dom,
                                               a_exp ) ;

   PEL_CHECK_POST( result!=0 && result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->field()==a_field ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
PDE_ValuesOnMeshingProjector:: PDE_ValuesOnMeshingProjector(
                      PEL_Object* a_owner,
                      PDE_DiscreteField* a_field,
                      size_t a_field_level,
                      PDE_ValuesOnMeshing* a_field_initializer,
                      PDE_DomainAndFields const* a_dom,
                      PEL_ModuleExplorer* a_exp )
//----------------------------------------------------------------------
   : PDE_ProjectorForDOFsSetting( a_owner, a_field, a_field_level, a_dom, a_exp )
   , INITIALIZER( a_field_initializer )
   , POLY( 0 )
{
   PEL_LABEL( "PDE_ValuesOnMeshingProjector:: PDE_ValuesOnMeshingProjector" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PDE_ValuesOnMeshingProjector:: ~PDE_ValuesOnMeshingProjector( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ValuesOnMeshingProjector:: ~PDE_ValuesOnMeshingProjector" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
PDE_ValuesOnMeshingProjector:: compute_value_at_IP( PDE_LocalFEcell const* fe,
                                                    doubleVector& result ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ValuesOnMeshingProjector:: compute_value_at_IP" ) ;
   PEL_CHECK( compute_value_at_IP_PRE( fe, result ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( POLY==0 || POLY!=fe->polyhedron() )
   {
      POLY = fe->polyhedron() ;
      INITIALIZER->set_mesh( POLY, fe->color() ) ;
   }
   INITIALIZER->compute_value( fe->coordinates_of_IP(), result ) ;
}

//----------------------------------------------------------------------
bool
PDE_ValuesOnMeshingProjector:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PDE_ProjectorForDOFsSetting::invariant() ) ;
   return( true ) ;
}
