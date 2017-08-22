#include <AP_LinearElasticity.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>
#include <doubleArray2D.hh>
#include <doubleArray4D.hh>

#include <AP_KinematicState.hh>

//---------------------------------------------------------------------------
AP_LinearElasticity*
AP_LinearElasticity:: create( PEL_Object* a_owner,
                              PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_LinearElasticity:: create" ) ;

   AP_LinearElasticity* result = new AP_LinearElasticity( a_owner, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result  ) ;
}

//---------------------------------------------------------------------------
AP_LinearElasticity:: AP_LinearElasticity( PEL_Object* a_owner, 
                                           PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : AP_ConstitutiveLaw( a_owner )
   , YOUNG( exp->double_data( "Young_modulus" ) )
   , POISSON( exp->double_data( "Poisson_coefficient" ) ) 
{
}

//---------------------------------------------------------------------------
AP_LinearElasticity:: ~AP_LinearElasticity( void )
//---------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
AP_LinearElasticity:: update_S_dSdE( AP_KinematicState const* st,
                                     doubleArray2D& Pio2,
                                     doubleArray4D& dPio2dGL ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_LinearElasticity:: update_S_dSdE" ) ;
   PEL_CHECK_PRE( update_S_dSdE_PRE( st, Pio2, dPio2dGL ) ) ;

   doubleArray2D const& grad_disp = st->grad_disp() ;
  
   for( size_t i=0 ; i<3 ; ++i )
   {
      for( size_t j=0 ; j<3 ; ++j )
      {
         double xx = 0.0 ;
         for( size_t k=0 ; k<3 ; ++k )
         {
            for( size_t l=0 ; l<3 ; ++l )
            {
               dPio2dGL( i, j, k, l )  = compute_rigidity_matrix( 
                                         i, j, k, l, YOUNG, POISSON ) ;
               xx += dPio2dGL( i, j, k, l ) 
                   * CauchyTensor( grad_disp, k, l ) ;
            }
	 }
         Pio2( i, j ) = xx ; 
      }
   }
}

//----------------------------------------------------------------------
double
AP_LinearElasticity:: kronecker( size_t i, size_t j )
         
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_LinearElasticity:: kronecker" ) ;

   if( i == j ) return 1 ;

   else return 0 ;
}

//----------------------------------------------------------------------
double
AP_LinearElasticity:: I( size_t i, size_t j, size_t k, size_t l )
         
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_LinearElasticity:: I" ) ;

   double result = 1. / 2. * ( kronecker( i, k ) * kronecker( j, l ) + 
                               kronecker( i, l ) * kronecker( j, k ) ) ;

   return result ;
}

//----------------------------------------------------------------------
double
AP_LinearElasticity:: J( size_t i, size_t j, size_t k, size_t l )         
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_LinearElasticity:: J" ) ;

   double result = 1. / 3. * kronecker( i, j ) * kronecker( k, l ) ; ;

   return result ;
}

//----------------------------------------------------------------------
double
AP_LinearElasticity:: K( size_t i, size_t j, size_t k, size_t l )         
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_LinearElasticity:: K" ) ;

   double result = I( i, j, k, l ) - J( i, j, k, l ) ;

   return result ;
}

//----------------------------------------------------------------------
double
AP_LinearElasticity:: CauchyTensor( doubleArray2D const& grad_disp, 
                                    size_t i, size_t j )
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_LinearElasticity:: CauchyTensor" ) ;

   double result = 0.5 * ( grad_disp( i, j ) + grad_disp( j, i ) ) ;

   return result ;
}

//----------------------------------------------------------------------
double
AP_LinearElasticity:: compute_rigidity_matrix( size_t i, size_t j, 
                                               size_t k, size_t l, 
                                               double young,
                                               double poisson ) const 
//----------------------------------------------------------------------
{
   PEL_LABEL( "AP_LinearElasticity:: compute_rigidity_matrix" ) ;
 
   double result = young / ( 1 + poisson ) * K( i, j, k, l ) +
                   young / ( 1 - 2 * poisson ) * J( i, j, k, l ) ;
   return result ;
}
