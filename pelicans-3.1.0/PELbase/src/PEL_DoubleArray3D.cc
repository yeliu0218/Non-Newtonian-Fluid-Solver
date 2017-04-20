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

#include <PEL_DoubleArray3D.hh>

#include <iostream>
#include <numeric>

#include <PEL_assertions.hh>

//----------------------------------------------------------------------------
PEL_DoubleArray3D*
PEL_DoubleArray3D:: create( PEL_Object* a_owner,
                          doubleArray3D const& val )
//----------------------------------------------------------------------------
{
   return( new PEL_DoubleArray3D( a_owner, val ) ) ;
}

//----------------------------------------------------------------------------
PEL_DoubleArray3D*
PEL_DoubleArray3D:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   return( new PEL_DoubleArray3D( a_owner, myValue ) ) ;
}

//----------------------------------------------------------------------------
PEL_DoubleArray3D:: PEL_DoubleArray3D( PEL_Object* a_owner,
                                   doubleArray3D const& val )
//----------------------------------------------------------------------------
   : PEL_Data( a_owner ),
     myValue( val )
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
PEL_DoubleArray3D:: ~PEL_DoubleArray3D( void )
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------------
PEL_Data::Type
PEL_DoubleArray3D:: data_type( void ) const
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
   return( PEL_Data::DoubleArray3D ) ;
}

//----------------------------------------------------------------------------
doubleArray3D const&
PEL_DoubleArray3D:: to_double_array3D( PEL_Context const* ct ) const
//----------------------------------------------------------------------------
{
   PEL_CHECK_PRE( to_double_array3D_PRE( ct ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( myValue ) ;
}

//----------------------------------------------------------------------------
void
PEL_DoubleArray3D:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
   std::string space( indent_width, ' ' ) ;
   std::ios::fmtflags oldoptions = os.flags( std::ios::scientific ) ;
   PEL_CHECK_INV( invariant() ) ;
   os << space << "array( " << std::endl << space ;
   for( size_t i=0 ; i<myValue.index_bound(0) ; i++ )
   {
      os << space << "array( " << std::endl << space ;
      for( size_t j=0 ; j<myValue.index_bound(1) ; j++ )
      {
         os << "  < " ;
         for( size_t k=0 ; k<myValue.index_bound(2) ; k++ )
         {
            if( (k+1) % 10 == 0 )
            {
               os << std::endl << space ;
            }
            PEL::print_double( os, myValue(i,j,k) ) ;
            os << " " ;
         }
         os << ">" ;
         if( j!=myValue.index_bound(1)-1 )
         {
            os << "," ;
         }
      }
      os << "  ) " ;
      if( i!=myValue.index_bound(0)-1 )
      {
         os << "," ;
      }
      os << std::endl << space ;
   }
   os << " )" ;
   os.flags( oldoptions ) ;
}
