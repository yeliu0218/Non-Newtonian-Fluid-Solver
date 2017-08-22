#include <AP_ConstitutiveLaw.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>
#include <doubleArray2D.hh>
#include <doubleArray4D.hh>

#include <AP_LinearElasticity.hh>
#include <AP_StVenantKirchhoff.hh>
#include <AP_NeoHooke.hh>
#include <AP_MooneyRivlin.hh>
#include <AP_LawL2_TEST.hh>
#include <AP_KinematicState.hh>

#include <string>

//--------------------------------------------------------------------------
AP_ConstitutiveLaw*
AP_ConstitutiveLaw:: create( PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "AP_ConstitutiveLaw:: create" ) ;

   AP_ConstitutiveLaw* result = 0 ;
   std::string const& nn = exp->string_data( "concrete_name" ) ;
   if( nn == "AP_LinearElasticity" )
   {
      result = AP_LinearElasticity::create( a_owner, exp ) ;
   }
   else if( nn == "AP_StVenantKirchhoff" )
   {
      result = AP_StVenantKirchhoff::create( a_owner, exp ) ;
   }
   else if( nn == "AP_NeoHooke" )
   {
      result = AP_NeoHooke::create( a_owner, exp ) ;
   }
   else if( nn == "AP_MooneyRivlin" )
   {
      result = AP_MooneyRivlin::create( a_owner, exp ) ;
   }
   else if( nn == "AP_LawL2_TEST" )
   {
      result = AP_LawL2_TEST::create( a_owner, exp ) ;
   }
   else 
   {
      PEL_Error::object()->raise_plain( "unknown : " + nn ) ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result  ) ;
}

//--------------------------------------------------------------------------
AP_ConstitutiveLaw:: ~AP_ConstitutiveLaw( void )
//--------------------------------------------------------------------------
{
}

//--------------------------------------------------------------------------
AP_ConstitutiveLaw:: AP_ConstitutiveLaw( PEL_Object* a_owner )
//--------------------------------------------------------------------------
   : PEL_Object( a_owner )
{
}

//--------------------------------------------------------------------------
bool
AP_ConstitutiveLaw:: update_S_dSdE_PRE( AP_KinematicState const* st,
                                        doubleArray2D& Pio2,
                                        doubleArray4D& dPio2dGL ) const
//--------------------------------------------------------------------------
{
   PEL_ASSERT( st != 0 ) ;
   PEL_ASSERT( Pio2.index_bound( 0 ) == 3 ) ;
   PEL_ASSERT( Pio2.index_bound( 1 ) == 3 ) ;
   PEL_ASSERT( dPio2dGL.index_bound( 0 ) == 3 ) ;
   PEL_ASSERT( dPio2dGL.index_bound( 1 ) == 3 ) ;
   PEL_ASSERT( dPio2dGL.index_bound( 2 ) == 3 ) ;
   PEL_ASSERT( dPio2dGL.index_bound( 3 ) == 3 ) ;
   return( true ) ;
}

