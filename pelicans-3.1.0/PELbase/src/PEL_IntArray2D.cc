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

#include <PEL_IntArray2D.hh>

#include <iostream>
#include <numeric>

#include <PEL_assertions.hh>
#include <intVector.hh>
#include <PEL.hh>
#include <PEL_List.hh>

//----------------------------------------------------------------------------
PEL_IntArray2D*
PEL_IntArray2D:: create( PEL_Object* a_owner,
                         intArray2D const& val )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_IntArray2D:: create(PEL_Object*,intArray2D const&)" ) ;

   PEL_IntArray2D* result = new PEL_IntArray2D( a_owner, val ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_IntArray2D*
PEL_IntArray2D:: create( PEL_Object* a_owner,
                         size_t_array2D const& val )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_IntArray2D:: create(PEL_Object*,intArray2D const&)" ) ;

   PEL_IntArray2D* result = new PEL_IntArray2D( a_owner, val ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_IntArray2D*
PEL_IntArray2D:: create( PEL_Object* a_owner,
                            PEL_List const* list )
//----------------------------------------------------------------------------
{
   PEL_CHECK_PRE( list!=0 ) ;
   PEL_CHECK_PRE( list->count()!=0 ) ;
   PEL_CHECK_PRE( FORALL( (size_t i=0 ; i<list->count() ; i++ ),
                          dynamic_cast<PEL_Data const*>( list->at(i) )!=0 ) ) ;

   size_t nb_rows = list->count() ;
   size_t nb_cols = PEL::bad_index() ;
   bool error = nb_rows==0 ;
   intArray2D array(0,0);
   if( !error )
   {   
      for( size_t i=0 ; !error && i<nb_rows ; i++ )
      {
         PEL_Data const* val =
            static_cast<PEL_Data const*>( list->at(i) ) ;
         error = val==0 || !val->value_can_be_evaluated(0) || val->data_type()!=IntVector ;
         if( !error )
         {
            intVector const& vec = val->to_int_vector() ;
            if( i==0 )
            {
               nb_cols = vec.size() ;
               error = nb_cols==0 ;
               if( !error )
               {
                  array.re_initialize(nb_rows,nb_cols);
               }
            }
            error = error || vec.size()!=nb_cols ;
            for( size_t j=0 ; !error && j<nb_cols ; j++ )
            {
               array(i,j) = vec(j) ;
            }
         }
      }
   }

   return ( error ? 0 : new PEL_IntArray2D( a_owner, array ) ) ;
   
}

//----------------------------------------------------------------------------
PEL_IntArray2D*
PEL_IntArray2D:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   return( new PEL_IntArray2D( a_owner, myValue ) ) ;
}

//----------------------------------------------------------------------------
PEL_IntArray2D:: PEL_IntArray2D( PEL_Object* a_owner,
                                 intArray2D const& val )
//----------------------------------------------------------------------------
   : PEL_Data( a_owner )
   , myValue( val )
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
PEL_IntArray2D:: PEL_IntArray2D( PEL_Object* a_owner,
                                 size_t_array2D const& val )
//----------------------------------------------------------------------------
   : PEL_Data( a_owner )
   , myValue( 0, 0 )
{
   myValue.set( val ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
PEL_IntArray2D:: ~PEL_IntArray2D( void )
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}



//----------------------------------------------------------------------------
PEL_Data::Type
PEL_IntArray2D:: data_type( void ) const
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
   return( PEL_Data::IntArray2D ) ;
}



//----------------------------------------------------------------------------
intArray2D const&
PEL_IntArray2D:: to_int_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------------
{
   PEL_CHECK_PRE( to_int_array2D_PRE(ct) ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( myValue ) ;
}



//----------------------------------------------------------------------------
void
PEL_IntArray2D:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
   std::string space( indent_width, ' ' ) ;
   PEL_CHECK_INV( invariant() ) ;
   os << space << "[ " << std::endl << space ;
   for( size_t i=0 ; i<myValue.index_bound(0) ; i++ )
   {
      os << "< " ;
      for( size_t j=0 ; j<myValue.index_bound(1) ; j++ )
      {
         if( (j+1) % 10 == 0 )
         {
	   os << std::endl << space ;
         }
         os << myValue(i,j) << " " ;
      }
      os << "  >" ;
      if( i!=myValue.index_bound(0)-1 )
      {
         os << "," ;
      }
      os << std::endl << space ;
   }
   os << " ]" ;
}
