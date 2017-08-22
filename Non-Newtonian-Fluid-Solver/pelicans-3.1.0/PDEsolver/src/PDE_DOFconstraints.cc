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

#include <PDE_DOFconstraints.hh>

#include <PEL.hh>

#include <PDE_DiscreteField.hh>

#include <iomanip>
#include <iostream>

using std::endl ;
using std::string ;
using std::vector ;

//----------------------------------------------------------------------
PDE_DOFconstraints:: PDE_DOFconstraints( PEL_Object* a_owner,
                                         size_t a_nb_nodes,
                                         size_t a_nb_comps )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , NB_NODES( a_nb_nodes )
   , NB_COMPS( a_nb_comps )
   , CONSTRAINT_IDX( NB_NODES, NB_COMPS, PEL::bad_index() ) 
{
}

//----------------------------------------------------------------------
PDE_DOFconstraints:: ~PDE_DOFconstraints( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PDE_DOFconstraints:: raise_nb_nodes( size_t a_nb_nodes )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFconstraints:: raise_nb_nodes" ) ;
   PEL_CHECK_PRE( a_nb_nodes > nb_nodes() ) ;
   
   NB_NODES = a_nb_nodes ;
   CONSTRAINT_IDX.raise_first_index_bound( NB_NODES, PEL::bad_index() ) ;
   
   PEL_CHECK_POST( nb_nodes() == a_nb_nodes ) ;
}

//----------------------------------------------------------------------
size_t
PDE_DOFconstraints:: nb_nodes( void ) const
//----------------------------------------------------------------------
{
   return( NB_NODES ) ;
}

//----------------------------------------------------------------------
size_t
PDE_DOFconstraints:: nb_components( void ) const
//----------------------------------------------------------------------
{
   return( NB_COMPS ) ;
}

//----------------------------------------------------------------------
void
PDE_DOFconstraints:: remove( size_t slave_n, size_t slave_ic )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFconstraints:: remove_constraint_for_DOF" ) ;
   PEL_CHECK_PRE( slave_n  < nb_nodes() ) ;
   PEL_CHECK_PRE( slave_ic < nb_components() ) ;
   PEL_CHECK_PRE( has( slave_n, slave_ic ) ) ;
   
   CONSTRAINT_IDX( slave_n, slave_ic ) = PEL::bad_index() ;
   
   PEL_CHECK_POST( !has( slave_n, slave_ic ) ) ; 
}

//----------------------------------------------------------------------
void
PDE_DOFconstraints:: add( size_t slave_n, size_t slave_ic,
                                            size_t master_n, size_t master_ic, 
                                            double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFconstraints:: add_constraint_for_DOF" ) ;
   PEL_CHECK_PRE( slave_n < nb_nodes() ) ;
   PEL_CHECK_PRE( slave_ic < nb_components() ) ;
   PEL_CHECK_PRE( master_n < nb_nodes() ) ;
   PEL_CHECK_PRE( master_ic < nb_components() ) ;
   PEL_CHECK_PRE( ! has( master_n, master_ic ) ) ; //???
   
   size_t ict = CONSTRAINT_IDX( slave_n, slave_ic ) ;
   if( ict == PEL::bad_index() )
   {
      CONSTRAINT_IDX( slave_n, slave_ic ) = CONSTRAINTS.size() ;
      CONSTRAINTS.push_back( 
         vector< ConstraintElement >( 1, 
               ConstraintElement( master_n, master_ic, coef ) ) ) ;
   }
   else
   {
      vector< ConstraintElement >& constraint = CONSTRAINTS[ ict ] ;
      for( size_t i=0 ; i<constraint.size() ; ++i )
      {
         //???? peut ne pas etre vrai si on raffine en haut, puis en bas etc
         //???? si ca avait déjà été raffiné avant ?????????????????
         PEL_ASSERT( ( constraint[i].n != master_n ) ||
                     ( ( constraint[i].n == master_n ) && 
                       ( constraint[i].ic != master_ic ) ) ) ;
      }
      constraint.push_back( ConstraintElement( master_n, master_ic, coef ) ) ;
   }
   
   PEL_CHECK_POST( has( slave_n, slave_ic ) ) ; 
}

//----------------------------------------------------------------------
bool
PDE_DOFconstraints:: has( size_t n, size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFconstraints:: DOF_is_constrained" ) ;
   PEL_CHECK_PRE( n < nb_nodes() ) ;
   PEL_CHECK_PRE( ic < nb_components() ) ;
   
   bool result = ( ( CONSTRAINT_IDX.index_bound( 0 ) != 0 ) &&
                   ( CONSTRAINT_IDX( n, ic ) != PEL::bad_index() ) ) ;
   return( result ) ;
}

