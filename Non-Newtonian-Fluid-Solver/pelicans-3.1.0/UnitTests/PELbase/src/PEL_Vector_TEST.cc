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

#include <PEL_Vector_TEST.hh>

#include <PEL_assertions.hh>
#include <PEL_Iterator.hh>
#include <PEL_Vector.hh>

#include <PEL_Double.hh>


#include <PEL_Root.hh>

//-------------------------------------------------------------------------
PEL_Vector_TEST * PEL_Vector_TEST::registered_test = new PEL_Vector_TEST() ;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
PEL_Vector *
PEL_Vector_TEST:: emptyList( PEL_Object * a_owner ) const
//-------------------------------------------------------------------------
{
   return PEL_Vector::create( a_owner, 0 ) ;
}



//-------------------------------------------------------------------------
void
PEL_Vector_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   // Base ancestor tests
   PEL_Sequence_TEST::process_all_tests() ;

   // Test of PEL_Vector only fonctionalities
   const size_t n = 10 ;
   
   PEL_Vector * vector = PEL_Vector::create( this, n ) ;
   for( size_t i = 0 ; i<n ; i++ )
   {
      vector->set_at( i, PEL_Double::create( this, i ) ) ;
   }

   size_t m = n/2 ;
   vector->resize( m ) ;
   bool ok = true ;
   for( size_t i = 0 ; ok && i<m ; i++ )
   {
      PEL_Object* o = vector->at( i ) ;
      PEL_Double * d = dynamic_cast<PEL_Double * >(o) ;
      ok = ok && d->to_double()==(double)i ;
   }   
   notify_one_test_result( "resize", ok && vector->index_limit()==m ) ;

   vector->nullify() ;
   ok = vector->count()==0 && vector->index_limit()==m ;
   for( size_t i=0 ; ok && i<m ; i++ )
   {
      ok = ok && vector->at( i )==0 ;
   }   
   notify_one_test_result( "nullify", ok ) ;

   
}



//-------------------------------------------------------------------------
PEL_Vector_TEST:: PEL_Vector_TEST( void )
//-------------------------------------------------------------------------
   :   PEL_Sequence_TEST( "PEL_Vector", "PEL_Vector_TEST" )
{
}



//-------------------------------------------------------------------------
PEL_Vector_TEST:: ~PEL_Vector_TEST( void )
//-------------------------------------------------------------------------
{
}



