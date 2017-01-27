#include <AP_KinematicState.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_assertions.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_LocalFE.hh>

#include <iostream>

//---------------------------------------------------------------------------
AP_KinematicState*
AP_KinematicState:: create( PEL_Object* a_owner,
                            size_t nb_disp_dims )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState:: create" ) ;
// ??? RegressionTests/FSI/Onde1D ???
   //PEL_CHECK_PRE( nb_disp_dims==2 || nb_disp_dims==3 ) ;

   AP_KinematicState* result = 
                  new AP_KinematicState( a_owner, nb_disp_dims ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_displacement_dimensions() == nb_disp_dims ) ;
   return( result  ) ;
}
//---------------------------------------------------------------------------
AP_KinematicState:: AP_KinematicState( PEL_Object* a_owner, 
                                       size_t nb_disp_dims )
//---------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , NB_DISP_DIMS( nb_disp_dims )
   , grad_DISP( 3, 3 )
   , CALC_F_J( false )
   , F( 3, 3 )
   , J( PEL::bad_double() )
   , CALC_Fm1( false )
   , Fm1( 3, 3 )
   , CALC_C( false )
   , C( 3, 3 )
   , CALC_Cm1( false )
   , Cm1( 3, 3 )
   , CALC_E( false )
   , E( 3, 3 )
{
}

//---------------------------------------------------------------------------
AP_KinematicState:: ~AP_KinematicState( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
AP_KinematicState:: set_state( PDE_DiscreteField const* uu,
                               size_t level,
                               PDE_LocalFE const* fe )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState:: set_state" ) ;
   PEL_CHECK_PRE( uu != 0 ) ;
   PEL_CHECK_PRE( level < uu->storage_depth() ) ;
   PEL_CHECK_PRE( uu->nb_components() == nb_displacement_dimensions() ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( fe->nb_local_nodes( uu ) != 0 )  ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled( uu, PDE_LocalFE::dN ) ) ;

   for( size_t ic=0 ; ic<NB_DISP_DIMS ; ++ic )
   {
      for( size_t d=0 ; d<NB_DISP_DIMS ; ++d )
      {
         grad_DISP( ic, d ) = fe->gradient_at_IP( uu, level, d, ic ) ;
      }
   }
   if( NB_DISP_DIMS == 2 )
   {
      for( size_t i=0 ; i<3 ; ++i )
      {
         grad_DISP( i, 2 ) = 0.0 ;
         grad_DISP( 2, i ) = 0.0 ;
      }
   }

   CALC_F_J = true ;
   CALC_Fm1 = true ;
   CALC_C   = true ;
   CALC_Cm1 = true ;
   CALC_E   = true ;
}

//---------------------------------------------------------------------------
void
AP_KinematicState:: set_state( doubleArray2D const& g_disp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState:: set_state" ) ;
   PEL_CHECK_PRE( g_disp.index_bound( 0 ) == nb_displacement_dimensions() ) ;
   PEL_CHECK_PRE( g_disp.index_bound( 1 ) == nb_displacement_dimensions() ) ;

   grad_DISP = g_disp ;

   CALC_F_J = true ;
   CALC_Fm1 = true ;
   CALC_C   = true ;
   CALC_Cm1 = true ;
   CALC_E = true ;
}

//---------------------------------------------------------------------------
size_t
AP_KinematicState:: nb_displacement_dimensions( void ) const
//---------------------------------------------------------------------------
{
   return( NB_DISP_DIMS ) ;
}

//---------------------------------------------------------------------------
doubleArray2D const&
AP_KinematicState:: grad_disp( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState:: grad_disp" ) ;

   doubleArray2D const& result = grad_DISP ;

   PEL_CHECK_POST( result.index_bound( 0 ) == 3 ) ;
   PEL_CHECK_POST( result.index_bound( 1 ) == 3 ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
doubleArray2D const&
AP_KinematicState:: DG( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState:: DG" ) ;

   if( CALC_F_J ) compute_F_J() ;
   
   doubleArray2D const& result = F ;

   return( result ) ;
}

//---------------------------------------------------------------------------
double
AP_KinematicState:: detF( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState:: detF" ) ;
   if( CALC_F_J ) compute_F_J() ;

   return( J ) ;
}

//---------------------------------------------------------------------------
doubleArray2D const&
AP_KinematicState:: inv_F( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState:: inv_F" ) ;

   if( CALC_Fm1 )
   {
      if( CALC_F_J ) compute_F_J() ;
      compute_Fm1() ;
   }

   doubleArray2D const& result = Fm1 ;

   return( result ) ;
}

//---------------------------------------------------------------------------
doubleArray2D const&
AP_KinematicState:: rCG( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState:: rCG" ) ;

   if( CALC_C )
   {
      if( CALC_F_J ) compute_F_J() ;
      compute_C() ;
   }

   doubleArray2D const& result = C ;

   PEL_CHECK_POST( result.index_bound( 0 ) == 3 ) ;
   PEL_CHECK_POST( result.index_bound( 1 ) == 3 ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
doubleArray2D const&
AP_KinematicState:: inv_rCG( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState:: inv_rCG" ) ;

   if( CALC_Cm1 )
   {
      if( CALC_C ) compute_C() ;
      compute_Cm1() ;
   }

   doubleArray2D const& result = Cm1 ;

   PEL_CHECK_POST( result.index_bound( 0 ) == 3 ) ;
   PEL_CHECK_POST( result.index_bound( 1 ) == 3 ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
doubleArray2D const&
AP_KinematicState:: GL( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState:: GL" ) ;

   if( CALC_E ) compute_E() ;

   doubleArray2D const& result = E ;

   PEL_CHECK_POST( result.index_bound( 0 ) == 3 ) ;
   PEL_CHECK_POST( result.index_bound( 1 ) == 3 ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
AP_KinematicState:: compute_F_J( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState:: compute_F_J" ) ;

   for( size_t a=0 ; a<NB_DISP_DIMS ; ++a )
   {
      for( size_t b=0 ; b<NB_DISP_DIMS ; ++b )
      {
         F( a, b ) = grad_DISP( a, b ) ;

         if( a == b ) F( a, b ) += 1.0 ;
      }
   }
   if( NB_DISP_DIMS == 2 )
   {
      F( 0, 2 ) = 0.0 ;
      F( 1, 2 ) = 0.0 ;
      F( 2, 2 ) = 1.0 ;
      F( 2, 0 ) = 0.0 ;
      F( 2, 1 ) = 0.0 ;
   }   

   J = determinant( F, NB_DISP_DIMS ) ;

   CALC_F_J = false ;
}

//---------------------------------------------------------------------------
void
AP_KinematicState:: compute_Fm1( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState:: compute_Fm1" ) ;

   if( NB_DISP_DIMS == 3 )
   {
      doubleArray2D cof( 3, 3 ) ;

      cof( 0, 0 ) = F( 1, 1 ) * F( 2, 2 ) - F( 1, 2 ) * F( 2, 1 ) ;
      cof( 0, 1 ) = F( 1, 2 ) * F( 2, 0 ) - F( 1, 0 ) * F( 2, 2 ) ;
      cof( 0, 2 ) = F( 1, 0 ) * F( 2, 1 ) - F( 1, 1 ) * F( 2, 0 ) ;
      cof( 1, 0 ) = F( 2, 1 ) * F( 0, 2 ) - F( 2, 2 ) * F( 0, 1 ) ;
      cof( 1, 1 ) = F( 2, 2 ) * F( 0, 0 ) - F( 2, 0 ) * F( 0, 2 ) ; 
      cof( 1, 2 ) = F( 2, 0 ) * F( 0, 1 ) - F( 2, 1 ) * F( 0, 0 ) ;
      cof( 2, 0 ) = F( 0, 1 ) * F( 1, 2 ) - F( 0, 2 ) * F( 1, 1 ) ;
      cof( 2, 1 ) = F( 0, 2 ) * F( 1, 0 ) - F( 0, 0 ) * F( 1, 2 ) ;
      cof( 2, 2 ) = F( 0, 0 ) * F( 1, 1 ) - F( 0, 1 ) * F( 1, 0 ) ;
  
      for( size_t a=0 ; a<3 ; ++a )
      {
         for( size_t b=0 ; b<3 ; ++b )
         {
            Fm1( a, b ) =  cof( b, a ) / J ;
         }
      }
   }
   else if( NB_DISP_DIMS == 2 )
   {
      Fm1( 0, 0 ) = F( 1, 1 ) / J ;
      Fm1( 1, 1 ) = F( 0, 0 ) / J ;
      Fm1( 0, 1 ) = - F( 0, 1 ) / J ;
      Fm1( 1, 0 ) = - F( 1, 0 ) / J ;
      Fm1( 0, 2 ) = 0.0 ;
      Fm1( 1, 2 ) = 0.0 ;
      Fm1( 2, 2 ) = 1.0 ;
      Fm1( 2, 0 ) = 0.0 ;
      Fm1( 2, 1 ) = 0.0 ;
   }

   CALC_Fm1 = false ;
}

//---------------------------------------------------------------------------
void
AP_KinematicState:: compute_C( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState:: compute_C" ) ;

   for( size_t a=0 ; a<NB_DISP_DIMS ; ++a )
   {
      for( size_t b=0 ; b<NB_DISP_DIMS ; ++b )
      {
         double xx = 0.0 ;
         for( size_t k=0 ; k<3 ; ++k )
         {
            xx += F( k, a ) * F( k, b ) ;
         }
         C( a, b ) = xx ;
      }
   }
   if( NB_DISP_DIMS == 2 )
   {
      C( 0, 2 ) = 0.0 ;
      C( 1, 2 ) = 0.0 ;
      C( 2, 2 ) = 1.0 ;
      C( 2, 0 ) = 0.0 ;
      C( 2, 1 ) = 0.0 ;
   }   

   CALC_C = false ; 
}

//---------------------------------------------------------------------------
void
AP_KinematicState:: compute_Cm1( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState:: compute_Cm1" ) ;

   double det = determinant( C, NB_DISP_DIMS ) ;

   if( NB_DISP_DIMS == 3 )
   {
      doubleArray2D cof( 3, 3 ) ;

      cof( 0, 0 ) = C( 1, 1 ) * C( 2, 2 ) - C( 1, 2 ) * C( 2, 1 ) ;
      cof( 0, 1 ) = C( 1, 2 ) * C( 2, 0 ) - C( 1, 0 ) * C( 2, 2 ) ;
      cof( 0, 2 ) = C( 1, 0 ) * C( 2, 1 ) - C( 1, 1 ) * C( 2, 0 ) ;
      cof( 1, 0 ) = C( 2, 1 ) * C( 0, 2 ) - C( 2, 2 ) * C( 0, 1 ) ;
      cof( 1, 1 ) = C( 2, 2 ) * C( 0, 0 ) - C( 2, 0 ) * C( 0, 2 ) ; 
      cof( 1, 2 ) = C( 2, 0 ) * C( 0, 1 ) - C( 2, 1 ) * C( 0, 0 ) ;
      cof( 2, 0 ) = C( 0, 1 ) * C( 1, 2 ) - C( 0, 2 ) * C( 1, 1 ) ;
      cof( 2, 1 ) = C( 0, 2 ) * C( 1, 0 ) - C( 0, 0 ) * C( 1, 2 ) ;
      cof( 2, 2 ) = C( 0, 0 ) * C( 1, 1 ) - C( 0, 1 ) * C( 1, 0 ) ;
  
      for( size_t a=0 ; a<3 ; ++a )
      {
         for( size_t b=0 ; b<3 ; ++b )
         {
            Cm1( a, b ) =  cof( b, a ) / det ;
         }
      }
   }
   else if( NB_DISP_DIMS == 2 )
   {
      Cm1( 0, 0 ) = C( 1, 1 ) / det ;
      Cm1( 1, 1 ) = C( 0, 0 ) / det ;
      Cm1( 0, 1 ) = - C( 0, 1 ) / det ;
      Cm1( 1, 0 ) = - C( 1, 0 ) / det ;
      Cm1( 0, 2 ) = 0.0 ;
      Cm1( 1, 2 ) = 0.0 ;
      Cm1( 2, 2 ) = 1.0 ;
      Cm1( 2, 0 ) = 0.0 ;
      Cm1( 2, 1 ) = 0.0 ;
   }

   CALC_Cm1 = false ;
}

//---------------------------------------------------------------------------
void
AP_KinematicState:: compute_E( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState:: compute_E" ) ;

   for( size_t a=0 ; a<NB_DISP_DIMS ; ++a )
   {
      for( size_t b=0 ; b<NB_DISP_DIMS ; ++b )
      {
         double xx = grad_DISP( a, b ) + grad_DISP( b, a ) ;
         for( size_t k=0 ; k<NB_DISP_DIMS ; ++k )
         {
            xx += grad_DISP( k, a ) * grad_DISP( k, b ) ;
         }
         E( a, b ) = 0.5 * xx ;
      }
   }
   if( NB_DISP_DIMS == 2 )
   {
      E( 0, 2 ) = 0.0 ;
      E( 1, 2 ) = 0.0 ;
      E( 2, 2 ) = 0.0 ;
      E( 2, 0 ) = 0.0 ;
      E( 2, 1 ) = 0.0 ;
   }   

   CALC_E = false ;
}

//---------------------------------------------------------------------------
double
AP_KinematicState:: determinant( doubleArray2D const& mat,
                                 size_t mat_size )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP_KinematicState:: determinant" ) ;
   PEL_CHECK( mat_size==2 || mat_size==3 ) ;

   double result = PEL::bad_double() ;

   if( mat_size == 3 )
   {
      result = mat( 0, 0 ) * mat( 1, 1 ) * mat( 2, 2 )
             + mat( 1, 0 ) * mat( 2, 1 ) * mat( 0, 2 )
             + mat( 2, 0 ) * mat( 0, 1 ) * mat( 1, 2 )
             - mat( 0, 2 ) * mat( 1, 1 ) * mat( 2, 0 )
             - mat( 0, 0 ) * mat( 1, 2 ) * mat( 2, 1 )
             - mat( 0, 1 ) * mat( 1, 0 ) * mat( 2, 2 ) ;
   }
   else if( mat_size == 2 )
   {
      result = mat( 0, 0 ) * mat( 1, 1 )
             - mat( 0, 1 ) * mat( 1, 0 ) ;
   }

   return result ;
}

