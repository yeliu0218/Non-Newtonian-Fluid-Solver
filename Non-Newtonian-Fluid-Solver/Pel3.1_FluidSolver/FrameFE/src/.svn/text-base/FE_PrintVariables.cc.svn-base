#include <FE_PrintVariables.hh>

#include <PEL.hh>
#include <PEL_Bool.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Int.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectReader.hh>
#include <PEL_ObjectWriter.hh>
#include <PEL_Variable.hh>
#include <PEL_assertions.hh>

#include <stringVector.hh>
#include <iostream>
#include <sstream>
using std::endl ;

//----------------------------------------------------------------------
FE_PrintVariables*
FE_PrintVariables:: create( PEL_Object* a_owner,
                          PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_PrintVariables:: create" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   
   FE_PrintVariables* result =  new FE_PrintVariables( a_owner,  exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ; 

   return( result ) ;
}

//----------------------------------------------------------------------
FE_PrintVariables:: FE_PrintVariables( PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
{
   PEL_LABEL( "FE_PrintVariables:: FE_PrintVariables" ) ;

     PEL::out()<< "****************************" << std::endl;
     PEL::out()<< "*** Print Variables      ***" << std::endl;
     PEL::out()<< "****************************" << std::endl;
     if( exp->has_entry( "names" ) )
     {
       stringVector name= exp->stringVector_data( "names" );
       doubleVector value = exp->doubleVector_data( "values" );
       
       for(int i=0;i<value.size();i++)
         PEL::out()<< name(i) << " = " << value(i) << std::endl;
       
      }
     PEL::out()<< "****************************" << std::endl;
}

//----------------------------------------------------------------------
FE_PrintVariables:: ~FE_PrintVariables( void )
//----------------------------------------------------------------------
{
}


