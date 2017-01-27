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

#include <PDE_SystemNumbering_TEST.hh>

#include <doubleVector.hh>
#include <stringVector.hh>

#include <PEL_assertions.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <PDE_CrossProcessNodeNumbering.hh>
#include <PDE_CrossProcessUnknownNumbering.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <fstream>
#include <ios>
#include <iostream>
#include <iomanip>
#include <sstream>

using std::endl ;
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;
using std::ostringstream ;

PDE_SystemNumbering_TEST const*
PDE_SystemNumbering_TEST:: REGISTRATOR = new PDE_SystemNumbering_TEST() ;

//---------------------------------------------------------------------------
PDE_SystemNumbering_TEST:: PDE_SystemNumbering_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "PDE_SystemNumbering", 
                     "PDE_SystemNumbering_TEST" )
{
}

//----------------------------------------------------------------------------
PDE_SystemNumbering_TEST:: ~PDE_SystemNumbering_TEST( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
PDE_SystemNumbering_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering_TEST:: process_one_test" ) ;

   TEST_NAME = exp->name() ;
   out() << "| ... "  <<  TEST_NAME << endl ;
   
   std::ofstream ofs( exp->string_data( "output_file" ).c_str() ) ;
   
   bool distributed = true ;
   if( exp->has_entry( "distributed" ) )
   {
      distributed = exp->bool_data( "distributed" ) ;
   }
   PEL_ModuleExplorer* ee = 
                       exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   PDE_DomainAndFields* dom = ( distributed ?
                PDE_DomainAndFields::create( 0, ee, PEL_Exec::communicator() ) :
                PDE_DomainAndFields::create( 0, ee ) ) ;
   ee->destroy() ; ee = 0 ;
   
   stringVector const& fnames = exp->stringVector_data( "fields_name" ) ;
   stringVector const& orders = exp->stringVector_data( "orderings" ) ;
   PEL_ASSERT( fnames.size() == orders.size() ) ;
   size_t nn = fnames.size() ;
   PEL_Vector* fields = PEL_Vector::create( 0, nn ) ;
   PEL_Vector* links  = PEL_Vector::create( 0, nn ) ;
   
   for( size_t i=0 ; i<nn ; ++i )
   {
      PDE_DiscreteField* ff = dom->set_of_discrete_fields()->item( fnames(i) ) ;
      fields->set_at( i, ff ) ;
   }
   
   // *** 
   print_disc( dom, fields, ofs ) ;
   
   // *** One instance of PDE_SystemNumbering for each field
   for( size_t i=0 ; i<nn ; ++i )
   {
      PDE_DiscreteField const* ff = 
                        static_cast< PDE_DiscreteField* >( fields->at( i ) ) ;
      PDE_LinkDOF2Unknown* ll = PDE_LinkDOF2Unknown::create( 0, ff,
                                                             orders( i ),
                                                             true ) ;
      ofs << "****************************************************" << endl ;
      ofs << "PDE_SystemNumbering::create( PEL_Object*, "           << endl ;
      ofs << "                             PDE_LinkDOF2Unknown* ) " << endl ;
      ofs << "****************************************************" << endl ;
      PDE_SystemNumbering const* nmb = PDE_SystemNumbering::create( 0, ll ) ;
      nmb->print( ofs, 0 ) ;
      nmb->destroy() ; nmb = 0 ;
   }
   
   // *** An instance of PDE_SystemNumbering mixing all fields
   for( size_t i=0 ; i<nn ; ++i )
   {
      PDE_DiscreteField const* ff = 
                        static_cast< PDE_DiscreteField* >( fields->at( i ) ) ;
      PDE_LinkDOF2Unknown* ll = PDE_LinkDOF2Unknown::create( 0, ff,
                                                             orders( i ),
                                                             true ) ;
      links->set_at( i, ll ) ;
   }
   std::string const& system_ordering = exp->string_data( "system_ordering" ) ;
   ofs << "*******************************************************" << endl ;
   ofs << "PDE_SystemNumbering::create( PEL_Object*, "              << endl ;
   ofs << "                             PEL_Vector const*, "        << endl ;
   ofs << "                             std::string const& ) "      << endl ;
   ofs << "*******************************************************" << endl ;
   PDE_SystemNumbering const* nmb = 
             PDE_SystemNumbering::create( 0, links, system_ordering ) ;
   nmb->print( ofs, 0 ) ;
   nmb->destroy() ; nmb = 0 ;
   
   // ***
   fields->destroy() ;
   links->destroy() ;
   dom->destroy() ;
}

//-----------------------------------------------------------------------------
void
PDE_SystemNumbering_TEST:: print_disc( PDE_DomainAndFields const* dom,
                                       PEL_Vector const* fields,
                                       std::ostream& os )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SystemNumbering_TEST:: print_disc" ) ;

   PDE_LocalFEcell* cfe = dom->create_LocalFEcell( 0 ) ;
   
   for( size_t i=0 ; i<fields->count() ; ++i )
   {
      PDE_DiscreteField const* ff = 
                        static_cast< PDE_DiscreteField* >( fields->at( i ) ) ;
      cfe->require_field_calculation( ff, PDE_LocalFE::N ) ;
   }
   
   for( cfe->start() ; cfe->is_valid() ; cfe->go_next() )
   {
      os << "------------------------------" << endl ;
      cfe->print( os, 0 ) ;
   }
   os << "------------------------------" << endl ;
   
   cfe->destroy() ; cfe = 0 ;
}
