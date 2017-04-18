#include <PEL_RunTest_DEMO1.hh>

#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>


#include <fstream>
#include <iostream>

using std::endl ;
using std::ofstream ;

PEL_RunTest_DEMO1 const* 
PEL_RunTest_DEMO1::PROTOTYPE = new PEL_RunTest_DEMO1() ;

//----------------------------------------------------------------------------
PEL_RunTest_DEMO1:: PEL_RunTest_DEMO1( void )
//----------------------------------------------------------------------------
   : PEL_Application( "PEL_RunTest_DEMO1" )
{
}

//----------------------------------------------------------------------------
PEL_RunTest_DEMO1*
PEL_RunTest_DEMO1:: create_replica( PEL_Object* a_owner, 
                                 PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest_DEMO1:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PEL_RunTest_DEMO1* result = new PEL_RunTest_DEMO1( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_RunTest_DEMO1:: PEL_RunTest_DEMO1( PEL_Object* a_owner, 
                                 PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , OUT_NAME( exp->string_data( "output_file" ) )
{
}

//----------------------------------------------------------------------------
PEL_RunTest_DEMO1:: ~PEL_RunTest_DEMO1( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
PEL_RunTest_DEMO1:: run( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_RunTest_DEMO1:: run" ) ;
   
   ofstream ofs( OUT_NAME.c_str() ) ;

   for( size_t i=0 ; i<4 ; ++i )
   {
      ofs << i << "," << 2.0*i << "," << 3.0*i << endl ;
   }

   ofs.close() ;
}

