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

#include <PEL_Int.hh>

#include <numeric>
#include <iostream>

#include <PEL_assertions.hh>

//----------------------------------------------------------------------------
PEL_Int*
PEL_Int:: create( PEL_Object* a_owner, int val )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Int:: create" ) ;
   return( new PEL_Int( a_owner, val ) ) ;
}



//----------------------------------------------------------------------------
PEL_Int*
PEL_Int:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Int:: create_clone" ) ;
   return( new PEL_Int( a_owner, myValue ) ) ;
}



//----------------------------------------------------------------------------
PEL_Int:: PEL_Int( PEL_Object* a_owner, int val )
//----------------------------------------------------------------------------
   : PEL_Data( a_owner ),
     myValue( val )
{
   PEL_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------------
PEL_Int:: ~PEL_Int( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Int:: ~PEL_Int" ) ;
   PEL_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------------
PEL_Data::Type
PEL_Int:: data_type( void ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Int:: data_type" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( PEL_Data::Int ) ;
}



//----------------------------------------------------------------------------
int
PEL_Int:: to_int( PEL_Context const* ct ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Int:: to_int" ) ;
   PEL_CHECK_PRE( to_int_PRE(ct) ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( myValue ) ;
}



//----------------------------------------------------------------------------
void
PEL_Int:: set( int val )
//----------------------------------------------------------------------------
{
   myValue = val ;
}



//----------------------------------------------------------------------------
void
PEL_Int:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Int:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;
   std::string space( indent_width, ' ' ) ;
   os << space << myValue ;
}
