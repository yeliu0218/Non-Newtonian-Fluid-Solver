#include <AP_LawTest.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <doubleVector.hh>
#include <doubleArray2D.hh>
#include <doubleArray4D.hh>

#include <AP_KinematicState.hh>
#include <AP_ConstitutiveLaw.hh>

#include <fstream>
#include <iostream>

using std::cout ;
using std::endl ;

AP_LawTest const* AP_LawTest::PROTOTYPE = new AP_LawTest() ;

//----------------------------------------------------------------------------
AP_LawTest:: AP_LawTest( void )
//----------------------------------------------------------------------------
   : PEL_Application( "AP_LawTest" )
{
}

//----------------------------------------------------------------------------
AP_LawTest*
AP_LawTest:: create_replica( PEL_Object* a_owner, 
                             PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LawTest:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   AP_LawTest* result = new AP_LawTest( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
AP_LawTest:: AP_LawTest( PEL_Object* a_owner, 
                         PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , LAW( 0 )
   , OUTPUT_FILE( exp->string_data( "output_file_name" ) )
{
   PEL_ModuleExplorer* se =
                       exp->create_subexplorer( 0, "AP_ConstitutiveLaw" ) ;
   LAW = AP_ConstitutiveLaw::create( this, se ) ;
   se->destroy() ;
}

//----------------------------------------------------------------------------
AP_LawTest:: ~AP_LawTest( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
AP_LawTest:: run( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LawTest:: run" ) ;
   
   size_t nb_dims = 3 ;
   size_t nbc = 3 ;

   doubleArray2D grad_disp( nb_dims, nb_dims ) ;
   doubleArray2D S( nb_dims, nb_dims ) ;
   doubleArray4D C( nb_dims, nb_dims, nb_dims, nb_dims ) ;

   std::ofstream os( OUTPUT_FILE.c_str() ) ;
   if( !os ) PEL_Error::object()->raise_file_handling( OUTPUT_FILE, "open" ) ;

   for( size_t a=0 ; a<nbc ; ++a )
   {
      for( size_t b=0 ; b<nb_dims ; ++b )
      {
         if( a == b ) grad_disp( a, b ) = 0.5 + a * b ;
         else grad_disp( a, b ) = 1.0 ;
      }
   }
 
   AP_KinematicState* st = AP_KinematicState::create( 0, nb_dims ) ;
   st->set_state( grad_disp ) ;
   LAW->update_S_dSdE( st, S, C ) ;
   st->destroy() ; st = 0 ;
   
   os << endl ;
   os << " 2 PK stress tensor S( i, j ) " << "and" 
      << " Rigidity matrix C( i, j, k, l ) for material : " << endl ;
   os << endl ;

   for( size_t i=0 ; i<nb_dims ; ++i )
   {
      for( size_t j=0 ; j<nb_dims ; ++j )
      {
         os << "  S( " << i << ", " << j << " )           " 
            << S( i, j ) << endl ;
      }
   }
   os << endl ;

   for( size_t i=0 ; i<nb_dims ; ++i )
   {
      for( size_t j=0 ; j<nb_dims ; ++j )
      {
         for( size_t k=0 ; k<nb_dims ; ++k )
         {
            for( size_t l=0 ; l<nb_dims ; ++l )
            {
               os << "  C( " << i << ", " << j << ", " 
                  << k << ", " << l << " )    " 
                  << C( i, j, k, l ) << endl ;
            }
	 }
      }
   }
   os.close() ;
}
