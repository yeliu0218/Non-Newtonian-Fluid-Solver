#include <AP_LawL2_TEST.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>

#include <doubleArray2D.hh>
#include <doubleArray4D.hh>

#include <AP_KinematicState.hh>

#include <iostream>

using std::cout ; using std::endl ;

//----------------------------------------------------------------------------
AP_LawL2_TEST*
AP_LawL2_TEST:: create( PEL_Object* a_owner,
                        PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LawL2_TEST:: create" ) ;

   AP_LawL2_TEST* result = new AP_LawL2_TEST( a_owner, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result  ) ;
}

//----------------------------------------------------------------------------
AP_LawL2_TEST:: AP_LawL2_TEST( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
   : AP_ConstitutiveLaw( a_owner )
{
}

//----------------------------------------------------------------------------
AP_LawL2_TEST:: ~AP_LawL2_TEST( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
AP_LawL2_TEST:: update_S_dSdE( AP_KinematicState const* st,
                               doubleArray2D& Pio2,
                               doubleArray4D& dPio2dGL ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LawL2_TEST:: update_S_dSdE" ) ;
   PEL_CHECK_PRE( update_S_dSdE_PRE( st, Pio2, dPio2dGL ) ) ;

   size_t dim = 3 ;

   doubleArray2D const& grad_disp = st->grad_disp() ;

   for( size_t i=0 ; i<dim ; ++i )
   { 
      for( size_t j=0 ; j<dim ; ++j )
      {
         Pio2( i, j ) = grad_disp( i, j ) ; 
      }
   }
}
