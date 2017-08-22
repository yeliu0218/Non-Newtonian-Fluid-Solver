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

#include <PEL_Double.hh>

#include <iostream>

#include <PEL_assertions.hh>
#include <PEL.hh>
#include <PEL_DoubleComparatorExact.hh>

//----------------------------------------------------------------------------
PEL_Double*
PEL_Double:: create( PEL_Object* a_owner, double val )
//----------------------------------------------------------------------------
{
   return( new PEL_Double( a_owner, val ) ) ;
}

//----------------------------------------------------------------------------
PEL_Double*
PEL_Double:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   return( new PEL_Double( a_owner, VALUE ) ) ;
}

//----------------------------------------------------------------------------
PEL_Double:: PEL_Double( PEL_Object* a_owner, double val )
//----------------------------------------------------------------------------
   : PEL_Data( a_owner )
   , VALUE( val )
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
PEL_Double:: ~PEL_Double( void )
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Double*
PEL_Double:: create_derivative( PEL_Object* a_owner,
                                PEL_Variable const* var,
                                PEL_Context const* ct  ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Double:: create_derivative" ) ;
   PEL_CHECK_PRE( create_derivative_PRE( a_owner, var, ct ) ) ;
   
   PEL_Double* result = create( a_owner, 0.0 ) ;
   
   PEL_CHECK_POST( create_derivative_POST( a_owner, var, result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_DoubleComparator const*
PEL_Double:: double_comparator( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Double:: double_comparator" ) ;
   static PEL_DoubleComparator const* result =
                                   PEL_DoubleComparatorExact::object() ;
   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_Double:: is_equal( PEL_Object const* other ) const 
//----------------------------------------------------------------------
{
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;

   bool result = ( three_way_comparison( other ) == 0 ) ;

   PEL_CHECK_POST( is_equal_POST( result, other ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
int
PEL_Double:: three_way_comparison( PEL_Object const* other ) const
//----------------------------------------------------------------------
{
   PEL_CHECK_PRE( three_way_comparison_PRE( other ) ) ;

   static PEL_DoubleComparator const* DBL_COMP = double_comparator() ;
   PEL_Double const* otherDouble = static_cast<PEL_Double const*>( other ) ;
   int result = DBL_COMP->three_way_comparison( VALUE, otherDouble->VALUE )  ;

   PEL_CHECK_POST( three_way_comparison_POST( result, other ) ) ;      
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PEL_Double:: hash_code( void ) const 
//----------------------------------------------------------------------
{
   return( (size_t) VALUE ) ;
}

//----------------------------------------------------------------------------
PEL_Data::Type
PEL_Double:: data_type( void ) const
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
   return( PEL_Data::Double ) ;
}

//----------------------------------------------------------------------------
double
PEL_Double:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------------
{
   PEL_CHECK_PRE( to_double_PRE( ct ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( VALUE ) ;
}

//----------------------------------------------------------------------------
void
PEL_Double:: set( double val )
//----------------------------------------------------------------------------
{
   VALUE = val ;
}

//----------------------------------------------------------------------------
void
PEL_Double:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
   std::string space( indent_width, ' ' ) ;
   os << space ;
   PEL::print_double( os, VALUE ) ;
}
