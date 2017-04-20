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

#include <PEL_Int_TEST.hh>

#include <PEL.hh>
#include <PEL_Int.hh>
#include <PEL_Map.hh>

//-------------------------------------------------------------------------
PEL_Int_TEST* PEL_Int_TEST::registered_test = new PEL_Int_TEST() ;
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void
PEL_Int_TEST:: process_all_tests( void )
//-------------------------------------------------------------------------
{
   int d = 2 ;
   PEL_Int const* p_d = PEL_Int::create( this, d ) ;
   notify_one_test_result("  to_int", p_d->to_int() == d ) ;
   PEL_Int* another_d = p_d->create_clone( this ) ;
   notify_one_test_result("  clone", another_d->to_int() == d ) ;
//    notify_one_test_result("  is_equal", another_d->is_equal( p_d ) ) ;
//    PEL_Int const* small_d = PEL_Int::create( this, d-1.0 ) ;
//    notify_one_test_result("  < = >", small_d->three_way_comparison( p_d ) == -1 ) ;
//    notify_one_test_result("  hash", another_d->hash_code()==p_d->hash_code() ) ;
}



//-------------------------------------------------------------------------
PEL_Int_TEST:: PEL_Int_TEST( void )
//-------------------------------------------------------------------------
   :   PEL_ObjectTest( "PEL_Int", "PEL_Int_TEST" )
{
}



//-------------------------------------------------------------------------
PEL_Int_TEST:: ~PEL_Int_TEST( void )
//-------------------------------------------------------------------------
{
}



