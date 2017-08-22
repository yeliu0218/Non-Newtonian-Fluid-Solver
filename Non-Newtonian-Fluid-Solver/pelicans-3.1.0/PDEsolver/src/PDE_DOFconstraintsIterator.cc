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

#include <PDE_DOFconstraintsIterator.hh>

#include <PDE_DOFconstraints.hh>

#include <iomanip>
#include <iostream>

using std::endl ;
using std::string ;
using std::vector ;

//----------------------------------------------------------------------
PDE_DOFconstraintsIterator:: PDE_DOFconstraintsIterator( PEL_Object* a_owner,
                                       PDE_DiscreteField const* a_ff,
                                       PDE_DOFconstraints const* a_cstr )
//----------------------------------------------------------------------
    : PEL_Object( a_owner )
    , FF( a_ff )
    , CSTR( a_cstr )
    , VEC_OF_IT( 0 )
    , OK_IT( false )
{
}

//----------------------------------------------------------------------
PDE_DOFconstraintsIterator:: ~PDE_DOFconstraintsIterator( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
PDE_DiscreteField const*
PDE_DOFconstraintsIterator:: field( void ) const
//----------------------------------------------------------------------
{
   return( FF ) ;
}

//----------------------------------------------------------------------
void
PDE_DOFconstraintsIterator:: start( size_t n, size_t ic )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFconstraintsIterator:: start_constraint_iterator" ) ;
   PEL_CHECK_PRE( n < field()->nb_nodes() ) ;
   PEL_CHECK_PRE( ic < field()->nb_components() ) ;
   PEL_CHECK_PRE( field()->DOF_is_constrained( n, ic ) ) ;
   
   VEC_OF_IT = & ( CSTR->CONSTRAINTS[ CSTR->CONSTRAINT_IDX( n, ic ) ] ) ;
   IT = VEC_OF_IT->begin() ;
   OK_IT = ( IT != VEC_OF_IT->end() ) ;
   
   PEL_CHECK_POST( is_valid() ) ;
}

//----------------------------------------------------------------------
bool
PDE_DOFconstraintsIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   return( OK_IT ) ;
}

//----------------------------------------------------------------------
void
PDE_DOFconstraintsIterator:: go_next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFconstraintsIterator:: go_next_constraint_element" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   
   ++IT ;
   OK_IT = ( IT != VEC_OF_IT->end() ) ;
}

//----------------------------------------------------------------------
size_t
PDE_DOFconstraintsIterator:: node_of_constraining_DOF( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFconstraintsIterator:: node_of_constraining_DOF" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   
   return( IT->n ) ;
}

//----------------------------------------------------------------------
size_t
PDE_DOFconstraintsIterator:: component_of_constraining_DOF( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFconstraintsIterator:: component_of_constraining_DOF" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   
   return( IT->ic ) ;
}

//----------------------------------------------------------------------
double
PDE_DOFconstraintsIterator:: constraint_coefficient( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DOFconstraintsIterator:: constraint_coefficient" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   
   return( IT->coef ) ;
}
