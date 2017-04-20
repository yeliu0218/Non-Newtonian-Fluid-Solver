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

#include <PEL_DoubleVector.hh>

#include <iostream>

#include <PEL_assertions.hh>
#include <PEL_Container.hh>
#include <PEL_Iterator.hh>
#include <doubleVector.hh>

//----------------------------------------------------------------------------
PEL_DoubleVector*
PEL_DoubleVector:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DoubleVector:: create_clone" ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_DoubleVector* result = 0 ;
   result = new PEL_DoubleVector( a_owner, dv ) ;
   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_DoubleVector*
PEL_DoubleVector:: create( PEL_Object* a_owner,
                           doubleVector const& aDoubleVector )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DoubleVector:: create(doubleVector)" ) ;
   PEL_DoubleVector* result = new PEL_DoubleVector( a_owner, aDoubleVector ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->to_double_vector().size() == aDoubleVector.size() ) ;
   PEL_CHECK_POST(
      FORALL( ( size_t i=0 ; i<aDoubleVector.size() ; ++i ),
              result->to_double_vector()(i) == aDoubleVector(i) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_DoubleVector:: PEL_DoubleVector( PEL_Object* a_owner,
                                     doubleVector const& aDoubleVector )
//----------------------------------------------------------------------------
   : PEL_Data( a_owner )
   , dv( aDoubleVector )
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
PEL_DoubleVector:: PEL_DoubleVector( PEL_Object* a_owner,
                                     PEL_Container const* aDataList,
                                     bool& error )
//----------------------------------------------------------------------------
   : PEL_Data( a_owner )
   , dv( aDataList->count() )
{
   PEL_Iterator* list_iterator = aDataList->create_iterator( 0 ) ;
   size_t cpt=0 ;
   for( list_iterator->start() ;
        list_iterator->is_valid() && !error ;
        list_iterator->go_next() )
   {
      PEL_Data const* val =
         dynamic_cast<PEL_Data const*>( list_iterator->item() ) ;
      error |= ( val==0  ||
                 val->data_type()!=Double ||
                 !val->value_can_be_evaluated( 0 ) ) ;
      if( !error ) dv(cpt++) = val->to_double() ;
   }
   
   list_iterator->destroy() ; list_iterator=0 ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
PEL_DoubleVector*
PEL_DoubleVector:: create( PEL_Object* a_owner,
                           PEL_Container const* aDataList ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DoubleVector:: create(PEL_Container)" ) ;
   PEL_CHECK_PRE( aDataList!=0 ) ;
   PEL_CHECK_PRE( aDataList->count()!=0 ) ;

   bool error = false ;
   PEL_DoubleVector* result = new PEL_DoubleVector( 0, aDataList, error ) ;
   if( error )
   {
      result->destroy() ;
      result = 0 ;
   }
   else if( a_owner != 0 )
   {
      result->set_owner( a_owner ) ;
   }

   PEL_CHECK_POST( IMPLIES( result != 0, result->owner() == a_owner ) ) ;
   PEL_CHECK_POST( IMPLIES( result != 0, result->to_double_vector().size() == aDataList->count() ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_DoubleVector:: ~PEL_DoubleVector( void )
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
void
PEL_DoubleVector:: set( doubleVector const& other )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DoubleVector:: set" ) ;
   dv=other ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<other.size() ; i++ ),
                           to_double_vector()(i)==other(i) ) ) ;
}

//----------------------------------------------------------------------------
PEL_Data::Type
PEL_DoubleVector:: data_type( void ) const
//----------------------------------------------------------------------------
{
   return( DoubleVector ) ;
}

//----------------------------------------------------------------------------
doubleVector const& 
PEL_DoubleVector:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DoubleVector:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE(ct) ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( dv ) ;
}

//----------------------------------------------------------------------------
void
PEL_DoubleVector:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DoubleVector:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   std::ios::fmtflags oldoptions = os.flags( std::ios::scientific ) ;
   PEL_CHECK_INV( invariant() ) ;
   os << space << "< " ;
   for( size_t i=0 ; i<dv.size() ; i++ )
   {
      PEL::print_double( os, dv(i) ) ;
      os << " " ;
   }
   os << "> " ;
   os.flags( oldoptions ) ;
}

//----------------------------------------------------------------------
PEL_Data*
PEL_DoubleVector:: create_derivative( PEL_Object* a_owner,
                                      PEL_Variable const* var,
                                      PEL_Context const* ct ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_DoubleVector:: create_derivative" ) ;
   PEL_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;

   doubleVector adv(dv.size()) ;
   adv.set( 0.0 ) ;
   
   PEL_DoubleVector* result = create( a_owner, adv ) ;
   
   PEL_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
bool
PEL_DoubleVector:: invariant( void ) const 
//-------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Data::invariant() ) ;
   return( true ) ;
}
