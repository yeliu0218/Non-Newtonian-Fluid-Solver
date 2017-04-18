#include <GE_Color_DEMO1.hh>

#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <GE_Color.hh>

#include <fstream>
#include <iostream>

using std::endl ;
using std::ofstream ;

GE_Color_DEMO1 const* GE_Color_DEMO1::PROTOTYPE = new GE_Color_DEMO1() ;

//----------------------------------------------------------------------------
GE_Color_DEMO1:: GE_Color_DEMO1( void )
//----------------------------------------------------------------------------
   : PEL_Application( "GE_Color_DEMO1" )
{
}

//----------------------------------------------------------------------------
GE_Color_DEMO1*
GE_Color_DEMO1:: create_replica( PEL_Object* a_owner, 
                                 PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Color_DEMO1:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   GE_Color_DEMO1* result = new GE_Color_DEMO1( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
GE_Color_DEMO1:: GE_Color_DEMO1( PEL_Object* a_owner, 
                                 PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , OUT_NAME( exp->string_data( "trace_file" ) )
{
}

//----------------------------------------------------------------------------
GE_Color_DEMO1:: ~GE_Color_DEMO1( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
GE_Color_DEMO1:: run( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Color_DEMO1:: run" ) ;
   
   ofstream ofs( OUT_NAME.c_str() ) ;

   GE_Color::extend( "blue" ) ;
   GE_Color::extend( "red" ) ;
   GE_Color::extend( "green" ) ;

   stringVector nl( 2 ) ;
   nl( 0 ) = "blue" ; nl( 1 ) = "red" ;
   GE_Color::extend( "bluered", nl ) ;

   nl( 0 ) = "green" ;
   GE_Color::extend( "greenred", nl ) ;

   GE_Color const* b  = GE_Color::object( "blue" ) ;
   GE_Color const* g  = GE_Color::object( "green" ) ;
   GE_Color const* br = GE_Color::object( "bluered" ) ;
   GE_Color const* gr = GE_Color::object( "greenred" ) ;

   ofs << "exists \"blue\" : " << GE_Color::exist( "blue" ) << endl ;
   ofs << "exists \"modry\" : " << GE_Color::exist( "modry" ) << endl ;
   ofs << "b->is_composite() : " << b->is_composite() << endl ;
   ofs << "br->is_composite() : " << br->is_composite() << endl ;

   ofs << "br->has( \"blue\" ) : " << br->has( "blue" ) << endl ;
   ofs << "br->has( \"modry\" ) : " << br->has( "modry" ) << endl ;

   ofs << "br->is_matching( b ) : " << br->is_matching( b ) << endl ;
   ofs << "br->is_matching( g ) : " << br->is_matching( g ) << endl ;
   ofs << "br->is_overlapping( b ) : " << br->is_overlapping( b ) << endl ;
   ofs << "br->is_overlapping( g ) : " << br->is_overlapping( g ) << endl ;
   ofs << "br->is_overlapping( gr ) : " << br->is_overlapping( gr ) << endl ;

   ofs.close() ;
}

