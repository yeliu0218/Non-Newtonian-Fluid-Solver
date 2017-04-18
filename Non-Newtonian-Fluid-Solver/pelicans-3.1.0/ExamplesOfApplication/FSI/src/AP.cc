#include <AP.hh>

#include <doubleVector.hh>
#include <doubleArray2D.hh>
#include <doubleArray4D.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEcell.hh>

#include <iostream>

using std::cout ; using std::endl ;

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

//---------------------------------------------------------------------------
void
AP:: add_C_gradsym_row_gradsym_col( PDE_LocalEquation* leq,
                                    PDE_LocalFEcell const* fe,
                                    doubleArray4D& dPio2dGL,
                                    double coef )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP:: add_C_gradsym_row_gradsym_col" ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( fe->field( row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( row ) == fe->field( col ) ) ;
   PEL_CHECK_PRE( fe->field( row )->nb_components() ==
                  fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                  fe->field( row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                  fe->nb_basis_functions( row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                  fe->nb_basis_functions( row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                  fe->field( row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                  fe->field( row )->nb_components() ) ;
//    PEL_CHECK_PRE( dPio2dGL.index_bound( 0 ) == fe->nb_space_dimensions() ) ;
//    PEL_CHECK_PRE( dPio2dGL.index_bound( 1 ) == fe->nb_space_dimensions() ) ;
//    PEL_CHECK_PRE( dPio2dGL.index_bound( 2 ) == fe->nb_space_dimensions() ) ;
//    PEL_CHECK_PRE( dPio2dGL.index_bound( 3 ) == fe->nb_space_dimensions() ) ;


   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_row_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_col_nodes = fe->nb_basis_functions( col ) ;

   double c_w = fe->weight_of_IP() * coef ;
   
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   doubleArray2D const& dN_col = fe->dNs_at_IP( col ) ;

   for( size_t i=0 ; i<nb_row_nodes ; ++i )
   {
      for( size_t j=0 ; j<nb_col_nodes ; ++j )
      {
         for( size_t ic=0 ; ic<nb_dims ; ++ic )
         {
            for( size_t jc=0 ; jc<nb_dims ; ++jc )
            {
               double xx = 0.0 ;             
               for( size_t a=0 ; a<nb_dims ; ++a )
               {
                  for( size_t b=0 ; b<nb_dims ; ++b )
                  {
                     double xx1 = 0.0 ; 
                     if( a == ic ) xx1 += dN_row( i, b ) ;
                     if( b == ic ) xx1 += dN_row( i, a ) ;
                     
                     for( size_t k=0 ; k<nb_dims ; ++k )
                     {
                        for( size_t l=0 ; l<nb_dims ; ++l )
                        {
                           double xx2 = 0.0 ;   
                           if( k == jc ) xx2 += dN_col( j, l ) ;
                           if( l == jc ) xx2 += dN_col( j, k ) ;
                           
                           xx += dPio2dGL( a, b, k, l ) * xx1 * xx2 ;
                        }
                     }
                  }
               }
               xx *= 0.25 * c_w ;
               leq->add_to_matrix( xx, i, j, ic, jc ) ;
            }
	 }
      }
   }
}

//---------------------------------------------------------------------------
void
AP:: add_S_gradsym_row( PDE_LocalEquation* leq,
                        PDE_LocalFEcell const* fe,
                        doubleArray2D& S,
                        double coef )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP:: add_S_gradsym_row" ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( fe->field( row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( col ) != 0 ) ; //?????? on n'ajoute que dans row
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == fe->nb_basis_functions( row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == fe->nb_basis_functions( col ) ) ;
   PEL_CHECK_PRE( fe->field( row )->nb_components() ==
                  fe->field( col )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                  fe->field( row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == 
                  fe->field( col )->nb_components() ) ;
//    PEL_CHECK_PRE( S.index_bound( 0 ) == fe->nb_space_dimensions() ) ;

   size_t nb_row_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_dims = fe->nb_space_dimensions()     ;

   double c_w = fe->weight_of_IP() * coef ;
 
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;

   for( size_t i=0 ; i<nb_row_nodes ; ++i )
   {
      for( size_t ic=0 ; ic<nb_dims ; ++ic )
      {
         double xx = 0.0 ;
         for( size_t a=0 ; a<nb_dims ; ++a )
         {
            for( size_t b=0 ; b<nb_dims ; ++b )
            {
               if( a == ic ) xx += S( a, b ) * dN_row( i, b ) ;
               if( b == ic ) xx += S( a, b ) * dN_row( i, a ) ;
            }
         }
  	 xx *= 0.5 * c_w ;
	 leq->add_to_vector( xx, i, ic ) ;
      }
   }
}

//---------------------------------------------------------------------------
void
AP:: add_C_gradGL_row_gradGL_col( PDE_LocalEquation* leq,
                                  PDE_LocalFEcell const* fe,
                                  doubleArray2D const& grad_disp,
                                  doubleArray4D const& dPio2dGL,
                                  double coef )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP:: add_C_gradGL_row_gradGL_col" ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( fe->field( row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( row ) == fe->field( col ) ) ;
   PEL_CHECK_PRE( fe->field( row )->nb_components() ==
                  fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                  fe->field( row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                  fe->nb_basis_functions( row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                  fe->nb_basis_functions( row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                  fe->field( row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                  fe->field( row )->nb_components() ) ;

   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;

   double c_w = fe->weight_of_IP() * coef ;
   
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;

   static doubleArray4D d_gl( nb_nodes, nb_dims, nb_dims, nb_dims ) ;
   d_gl.re_initialize( nb_nodes, nb_dims, nb_dims, nb_dims ) ;

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      for( size_t ic=0 ; ic<nb_dims ; ++ic )
      {
         for( size_t a=0 ; a<nb_dims ; ++a )
         {
            for( size_t b=0 ; b<nb_dims ; ++b )
            {
               double xx = grad_disp( ic, a ) * dN_row( i, b )
                         + grad_disp( ic, b ) * dN_row( i, a ) ;
               if( a == ic ) xx += dN_row( i, b ) ;
               if( b == ic ) xx += dN_row( i, a ) ;
               d_gl( i, ic, a, b ) = xx ;
            }
         }
      }
   }

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      for( size_t ic=0 ; ic<nb_dims ; ++ic )
      {
         for( size_t j=0 ; j<nb_nodes ; ++j )
         {
            for( size_t jc=0 ; jc<nb_dims ; ++jc )
            {
               double xx = 0.0 ;
               for( size_t a=0 ; a<nb_dims ; ++a )
               {
                  for( size_t b=0 ; b<nb_dims ; ++b )
                  {
                     double x1 = 0.0 ;
                     for( size_t k=0 ; k<nb_dims ; ++k )
                     {
                        for( size_t l=0 ; l<nb_dims ; ++l )
                        {
                           x1 += dPio2dGL( a, b, k, l ) * d_gl( j, jc, k, l ) ;
                        }
                     }
                     xx += d_gl( i, ic, a, b ) * x1 ;
                  }
               }
               xx *= 0.25 * c_w ;
               leq->add_to_matrix( xx, i, j, ic, jc ) ;
            }
	 }
      }
   }
}

//---------------------------------------------------------------------------
void
AP:: add_S_graddGL_row_col( PDE_LocalEquation* leq,
                            PDE_LocalFEcell const* fe,
                            doubleArray2D const& S,
                            double coef )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP:: add_S_graddGL_row_col" ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( fe->field( row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( col ) != 0 ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == fe->nb_basis_functions( row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == fe->nb_basis_functions( col ) ) ;
   PEL_CHECK_PRE( fe->field( row )->nb_components() ==
                  fe->field( col )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                  fe->field( row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == 
                  fe->field( col )->nb_components() ) ;
//    PEL_CHECK_PRE( S.index_bound( 0 ) == fe->nb_space_dimensions() ) ;
//    PEL_CHECK_PRE( S.index_bound( 1 ) == fe->nb_space_dimensions() ) ;

   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_row_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_col_nodes = fe->nb_basis_functions( col ) ;

   double c_w = fe->weight_of_IP() * coef ;
 
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   doubleArray2D const& dN_col = fe->dNs_at_IP( col ) ;

   for( size_t i=0 ; i<nb_row_nodes ; ++i )
   {
      for( size_t j=0 ; j<nb_col_nodes ; ++j )
      {        
         for( size_t ic=0 ; ic<nb_dims ; ++ic )
         {
            double xx = 0.0 ;
            for( size_t a=0 ; a<nb_dims ; ++a )
            {
               for( size_t b=0 ; b<nb_dims ; ++b )
               {
                  xx += S( a, b ) * ( dN_row( i, b ) * dN_col( j, a ) +
                                      dN_row( i, a ) * dN_col( j, b ) ) ;
               }
            }     
            xx *= 0.5 * c_w ;
            leq->add_to_matrix( xx, i, j, ic, ic ) ;
         }
      }
   } 
}

//---------------------------------------------------------------------------
void
AP:: add_S_gradGL_row( PDE_LocalEquation* leq,
                       PDE_LocalFEcell const* fe,
                       doubleArray2D const& grad_disp,
                       doubleArray2D const& S,
                       double coef )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "AP:: add_S_gradGL_row" ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( fe->field( row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( col ) != 0 ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == fe->nb_basis_functions( row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == fe->nb_basis_functions( col ) ) ;
   PEL_CHECK_PRE( fe->field( row )->nb_components() ==
                  fe->field( col )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                  fe->field( row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == 
                  fe->field( col )->nb_components() ) ;
//    PEL_CHECK_PRE( S.index_bound( 0 ) == fe->nb_space_dimensions() ) ;
//    PEL_CHECK_PRE( S.index_bound( 1 ) == fe->nb_space_dimensions() ) ;

   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_row_nodes = fe->nb_basis_functions( row ) ;

   double c_w = fe->weight_of_IP() * coef ;
 
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;

   for( size_t i=0 ; i<nb_row_nodes ; ++i )
   {
      for( size_t ic=0 ; ic<nb_dims ; ++ic )
      {
         double xx = 0.0 ;
         for( size_t a=0 ; a<nb_dims ; ++a )
         {
            for( size_t b=0 ; b<nb_dims ; ++b )
            {
               xx += S( a, b ) * ( grad_disp( ic, a ) * dN_row( i, b ) +
                                   grad_disp( ic, b ) * dN_row( i, a ) ) ;
               if( a == ic ) xx += S( a, b ) * dN_row( i, b ) ;
               if( b == ic ) xx += S( a, b ) * dN_row( i, a ) ;
            }
         }
  	 xx *= 0.5 * c_w ;
	 leq->add_to_vector( xx, i, ic ) ;
      }
   }
}
