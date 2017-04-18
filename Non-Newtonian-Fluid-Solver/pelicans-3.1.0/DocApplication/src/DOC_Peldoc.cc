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

#include <DOC_Peldoc.hh>

#include <iostream>

#include <PEL_assertions.hh>
#include <PEL_ModuleExplorer.hh>
#include <DOC_Category.hh>
#include <DOC_Class.hh>
#include <DOC_Tools.hh>
#include <DOC_Package.hh>
#include <DOC_Writer.hh>

#include <DOC_SyntaxTest.hh>

using std::cout ;
using std::endl ;

DOC_Peldoc const* DOC_Peldoc::PROTOTYPE = new DOC_Peldoc() ;

//----------------------------------------------------------------------------
DOC_Peldoc:: DOC_Peldoc( void )
//----------------------------------------------------------------------------
   : PEL_Application( "peldoc" )
{
}

//----------------------------------------------------------------------------
DOC_Peldoc*
DOC_Peldoc:: create_replica( PEL_Object* a_owner, 
                      PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Peldoc:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   DOC_Peldoc* result = new DOC_Peldoc( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
DOC_Peldoc:: run( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Peldoc:: run" ) ;
   
   /* Les files standards du C++ ne sont pas lus */
   DOC_Tools::declare_read( "cstddef" ) ;
   DOC_Tools::declare_read( "cstdio" ) ;
   DOC_Tools::declare_read( "cstdlib" ) ;
   DOC_Tools::declare_read( "ctime" ) ;
   DOC_Tools::declare_read( "exception" ) ;
   DOC_Tools::declare_read( "fstream" ) ;
   DOC_Tools::declare_read( "iostream" ) ;
   DOC_Tools::declare_read( "iosfwd" ) ;
   DOC_Tools::declare_read( "iomanip" ) ;
   DOC_Tools::declare_read( "new" ) ;
   DOC_Tools::declare_read( "numeric" ) ;
   DOC_Tools::declare_read( "string" ) ;
   DOC_Tools::declare_read( "typeinfo" ) ;
   
   DOC_Tools::read_files() ;
   
   DOC_Class::display_all() ;
   
   if( DOC_Tools::message() )
   {
      DOC_Category::display_list( std::cout, (size_t)0 ) ;
   }
   
}

//----------------------------------------------------------------------------
DOC_Peldoc:: DOC_Peldoc( PEL_Object* a_owner, PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
{
   if( exp->has_entry( "packaging" ) )
   {
      DOC_Package::read( exp->string_data( "packaging" ) ) ;
   }
   if( exp->has_entry( "upper" ) )
   {
      DOC_Class::declare_upper( exp->stringVector_data( "upper" ) ) ;
   }
   DOC_Writer::parse_explorer( exp ) ;
   DOC_Tools::parse_explorer( exp ) ;
}

//----------------------------------------------------------------------------
DOC_Peldoc:: ~DOC_Peldoc( void )
//----------------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}



