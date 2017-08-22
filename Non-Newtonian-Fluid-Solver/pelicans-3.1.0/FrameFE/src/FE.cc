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

#include <FE.hh>

#include <doubleArray2D.hh>
#include <doubleVector.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>

#include <PEL_assertions.hh>

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;

FE::FE_Geometry FE::GEOM = FE::unspecified ;

//---------------------------------------------------------------------------
FE::FE_Geometry
FE:: geometry( void )
//---------------------------------------------------------------------------
{
   return( GEOM ) ;
}

//---------------------------------------------------------------------------
void
FE:: set_geometry( FE_Geometry geom )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE:: set_geometry" ) ;
   PEL_CHECK_PRE( geom == FE::cartesian || geom == FE::axisymmetrical ) ;

   GEOM = geom ;

   PEL_CHECK_POST( FE::geometry() == geom ) ;
}

//---------------------------------------------------------------------------
double
FE:: bound_measure( PDE_LocalFEbound const* fe )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE:: bound_measure" ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( FE::geometry() != FE::unspecified ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;

   double result = fe->polyhedron()->measure() ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      result *= fe->polyhedron()->center()->coordinate(0) ;
   }

   return( result ) ;
}

//---------------------------------------------------------------------------
double
FE:: side_measure( PDE_CursorFEside const* fe )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE:: side_measure" ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( !fe->is_periodic() ) ;
   PEL_CHECK_PRE( FE::geometry() != FE::unspecified ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;

   GE_Mpolyhedron const* poly = fe->polyhedron() ;
   double result = poly->measure() ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      result *= poly->center()->coordinate( 0 ) ;
   }

   return( result ) ;
}

//---------------------------------------------------------------------------
double
FE:: side_measure( PDE_CursorFEside const* fe, size_t i_adj )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE:: side_measure" ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->is_periodic() ) ;
   PEL_CHECK_PRE( i_adj==0 || i_adj==1 ) ;
   PEL_CHECK_PRE( FE::geometry() != FE::unspecified ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;

   GE_Mpolyhedron const* poly = fe->polyhedron( i_adj ) ;
   double result = poly->measure() ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      result *= poly->center()->coordinate( 0 ) ;
   }

   return( result ) ;
}

//---------------------------------------------------------------------------
double
FE:: cell_measure( PDE_LocalFEcell const* fe )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE:: cell_measure" ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( FE::geometry() != FE::unspecified ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;

   double result = fe->polyhedron()->measure() ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      result *= fe->polyhedron()->center()->coordinate(0) ;
   }

   return( result ) ;
}

//---------------------------------------------------------------------------
void
FE:: add_row_col_S( PDE_LocalEquation* leq,
                    PDE_LocalFE const* fe,
                    double coef )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_row_col_S" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) == 
                     fe->field( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::N ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == 
                     fe->field( PDE_LocalFE::col )->nb_components() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double c_w = fe->weight_of_IP() * coef ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   doubleVector const& N_row = fe->Ns_at_IP( row ) ;
   doubleVector const& N_col = fe->Ns_at_IP( col ) ;
   
   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      for( size_t j=i ; j<nb_nodes ; ++j )
      {
         double xx = N_row( i ) * N_col( j ) ;
         xx *= c_w ;

         for( size_t ic=0 ; ic<nbc ; ++ic )
         {
            leq->add_to_matrix( xx, i, j, ic, ic ) ;
            if( i!=j ) leq->add_to_matrix( xx, j, i, ic, ic ) ;
         }
      }
   }
}

//---------------------------------------------------------------------------
void
FE:: add_row_col_NS( PDE_LocalEquation* leq,
                     PDE_LocalFE const* fe,
                     double coef )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_row_col_NS" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::N ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::col ), PDE_LocalFE::N ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                  fe->nb_basis_functions( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row )->nb_components() ==
                     fe->field( PDE_LocalFE::col )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == 
                     fe->field( PDE_LocalFE::col )->nb_components() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_row_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_col_nodes = fe->nb_basis_functions( col ) ;
   double c_w = fe->weight_of_IP() * coef ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   doubleVector const& N_row = fe->Ns_at_IP( row ) ;
   doubleVector const& N_col = fe->Ns_at_IP( col ) ;

   for( size_t i=0 ; i<nb_row_nodes ; ++i )
   {
      for( size_t j=0 ; j<nb_col_nodes ; ++j )
      {
         double xx = N_row( i ) * N_col( j ) ;
         xx *= c_w ;

         for( size_t ic=0 ; ic<nbc ; ++ic )
         {
            leq->add_to_matrix( xx, i, j, ic, ic ) ;
         }
      }
   }
}

//---------------------------------------------------------------------------
void
FE:: add_lumped_row_col( PDE_LocalEquation* leq,
                         PDE_LocalFE const* fe,
                         double coef )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_lumped_row_col" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) == fe->field( col ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::N ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == 
                     fe->field( PDE_LocalFE::col )->nb_components() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double c_w = fe->weight_of_IP() * coef ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   doubleVector const& N_row = fe->Ns_at_IP( row ) ;
   
   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      double xx = N_row( i ) * c_w ;

      for( size_t ic=0 ; ic<nbc ; ++ic )
      {
         leq->add_to_matrix( xx, i, i, ic, ic ) ;
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_grad_row_grad_col_S( PDE_LocalEquation* leq,
                              PDE_LocalFEcell const* fe,
                              double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_grad_row_grad_col_S" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) ==
                     fe->field( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row )->nb_components() == 1 || 
                  fe->field( PDE_LocalFE::row )->nb_components() == 
                     fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( IMPLIES(
                     FE::geometry() == FE::axisymmetrical &&
                        fe->field( PDE_LocalFE::row )->nb_components() != 1,
                     fe->field_calculation_is_handled(
                        fe->field( PDE_LocalFE::row ), PDE_LocalFE::N ) ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   
   size_t nbc = fe->field(row)->nb_components() ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double c_w = fe->weight_of_IP() * coef ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   doubleArray2D const& dN_col = fe->dNs_at_IP( col ) ;
   
   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      for( size_t j=i ; j<nb_nodes ; ++j )
      {
         double xx = 0.0 ;
         for( size_t d=0 ; d<nb_dims ; ++d )
         {
            xx += dN_row( i, d ) * dN_col( j, d ) ;
         }
         xx *= c_w ;
         for( size_t ic=0 ; ic < nbc ; ++ic )
         {
            leq->add_to_matrix( xx, i, j, ic, ic ) ;
            if( i != j ) leq->add_to_matrix( xx, j, i, ic, ic ) ;
         }
      }
   }

   if( FE::geometry() == FE::axisymmetrical && nbc != 1 )
   {
      double radius = fe->coordinates_of_IP()->coordinate(0) ;
      doubleVector const& N_row = fe->Ns_at_IP( row ) ;
      doubleVector const& N_col = fe->Ns_at_IP( col ) ;
      double yy = c_w / radius / radius ;
      
      for( size_t i=0 ; i<nb_nodes ; ++i )
      {
         for( size_t j=i ; j<nb_nodes ; ++j )
         {
            double xx = yy * N_row( i ) * N_col( j ) ;
            leq->add_to_matrix( xx, i, j, 0, 0 ) ;
            if( i != j ) leq->add_to_matrix( xx, j, i, 0, 0 ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_grad_row_grad_col_NS( PDE_LocalEquation* leq,
                               PDE_LocalFEcell const* fe,
                               double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_grad_row_grad_col_NS" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row )->nb_components() == 1 || 
                  fe->field( PDE_LocalFE::row )->nb_components() == 
                     fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col )->nb_components() ==
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( IMPLIES(
                     FE::geometry() == FE::axisymmetrical &&
                        fe->field( PDE_LocalFE::row )->nb_components() != 1,
                     fe->field_calculation_is_handled(
                        fe->field( PDE_LocalFE::row ), PDE_LocalFE::N ) ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::col ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( IMPLIES(
                     FE::geometry() == FE::axisymmetrical &&
                        fe->field( PDE_LocalFE::col )->nb_components() != 1,
                     fe->field_calculation_is_handled(
                        fe->field( PDE_LocalFE::col ), PDE_LocalFE::N ) ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                     fe->nb_basis_functions( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                     fe->field( PDE_LocalFE::col )->nb_components() ) ;

   size_t nbc = fe->field(row)->nb_components() ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_row_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_col_nodes = fe->nb_basis_functions( col ) ;
   double c_w = fe->weight_of_IP() * coef ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   doubleArray2D const& dN_col = fe->dNs_at_IP( col ) ;

   for( size_t i=0 ; i<nb_row_nodes ; ++i )
   {
      for( size_t j=0 ; j<nb_col_nodes ; ++j )
      {
         double xx = 0.0 ;
         for( size_t d=0 ; d<nb_dims ; ++d )
         {
            xx += dN_row( i, d ) * dN_col( j, d ) ;
         }
         xx *= c_w ;
         for( size_t ic=0 ; ic < nbc ; ++ic )
         {
            leq->add_to_matrix( xx, i, j, ic, ic ) ;
         }
      }
   }
   
   if( FE::geometry() == FE::axisymmetrical && nbc != 1 )
   {
      double radius = fe->coordinates_of_IP()->coordinate(0) ;
      doubleVector const& N_row = fe->Ns_at_IP( row ) ;
      doubleVector const& N_col = fe->Ns_at_IP( col ) ;
      double yy = c_w / radius / radius ;
      
      for( size_t i=0 ; i<nb_row_nodes ; ++i )
      {
         for( size_t j=0 ; j<nb_col_nodes ; ++j )
         {
            double xx = yy * N_row( i ) * N_col( j ) ;
            leq->add_to_matrix( xx, i, j, 0, 0 ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_row( PDE_LocalEquation* leq,
              PDE_LocalFE const* fe,
              double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_row" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row )->nb_components() == 1 ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::N ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;

   size_t nb_nodes = fe->nb_basis_functions( PDE_LocalFE::row ) ;
   double c_w = fe->weight_of_IP() * coef ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   doubleVector const& N_row = fe->Ns_at_IP( PDE_LocalFE::row ) ;

   for( size_t i=0 ; i<nb_nodes ; i++ )
   {
      double xx = c_w * N_row( i ) ;
      leq->add_to_vector( xx, i ) ;
   }
}

//----------------------------------------------------------------------
void
FE:: add_row( PDE_LocalEquation* leq,
              PDE_LocalFE const* fe,
              doubleVector const& coef_1,
              double coef_2 )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_row" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::N ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( coef_1.size() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   doubleVector const& N_row = fe->Ns_at_IP( row ) ;
   double w = coef_2 * fe->weight_of_IP() ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      w *= fe->coordinates_of_IP()->coordinate(0) ;
   }

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      double yy = w * N_row( i ) ;
      for( size_t ic=0 ; ic<nbc ; ++ic )
      {
         double xx = yy * coef_1( ic ) ;
         leq->add_to_vector( xx, i, ic ) ;
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_grad_row( PDE_LocalEquation* leq,
                   PDE_LocalFE const* fe,
                   doubleVector const& coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_grad_row" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical  ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( coef.size() == fe->nb_space_dimensions() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   double w = fe->weight_of_IP() ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      w *= fe->coordinates_of_IP()->coordinate(0) ;
   }

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      double xx = 0.0 ;
      for( size_t d=0 ; d<nb_dims ; ++d )
      {
         xx += dN_row( i, d ) * coef( d ) ;
      }
      xx *= w ;
      for( size_t ic=0 ; ic<nbc ; ++ic )
      {
         leq->add_to_vector( xx, i, ic ) ;
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_grad_row( PDE_LocalEquation* leq,
                   PDE_LocalFE const* fe,
                   doubleArray2D const& coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_grad_row" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical  ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( coef.index_bound( 0 ) == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( coef.index_bound( 1 ) == fe->nb_space_dimensions() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   double w = fe->weight_of_IP() ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
 
   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      for( size_t ic=0 ; ic<nbc ; ++ic )
      {
         double xx = 0.0 ;
         for( size_t d=0 ; d<nb_dims ; ++d )
         {
            xx += dN_row( i, d ) * coef( ic, d ) ;
         }
         xx *= w ;
         leq->add_to_vector( xx, i, ic ) ;
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_div_row( PDE_LocalEquation* leq,
                  PDE_LocalFE const* fe,
                  double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_div_row" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row )->nb_components() ==
                     fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->field_calculation_is_handled(
                              fe->field( PDE_LocalFE::row ),
                              PDE_LocalFE::N ) ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                  fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                  fe->field( PDE_LocalFE::row )->nb_components() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double c_w = fe->weight_of_IP() * coef ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      for( size_t ic=0 ; ic<nbc ; ++ic )
      {
         double xx = c_w*dN_row( i, ic ) ;
         leq->add_to_vector( xx, i, ic ) ;
      }
   }
   if( FE::geometry() == FE::axisymmetrical )
   {
      // Axisymmetrical : div(f) = df/dr+df/dz+f/r
      double radius = fe->coordinates_of_IP()->coordinate(0) ;
      double yy = c_w / radius ;
      doubleVector const& N_row = fe->Ns_at_IP( row ) ;
      for( size_t i=0 ; i<nb_nodes ; ++i )
      {
         double xx = yy * N_row( i ) ;
         leq->add_to_vector( xx, i, 0 ) ;
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_grad_row_D_col_S( PDE_LocalEquation* leq,
                           PDE_LocalFEcell const* fe, 
                           double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_grad_row_D_col_S" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) == 
                     fe->field( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row )->nb_components() ==
                     fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( IMPLIES(
                     FE::geometry() == FE::axisymmetrical &&
                        fe->field( PDE_LocalFE::row )->nb_components() != 1,
                     fe->field_calculation_is_handled(
                        fe->field( PDE_LocalFE::row ), PDE_LocalFE::N ) ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;

   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   double c_w = coef*fe->weight_of_IP() ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   doubleArray2D const& dN_col = fe->dNs_at_IP( col ) ;

   for( size_t i=0 ; i<nb_nodes ; i++ )
   {
      for( size_t j=0 ; j<nb_nodes ; j++ )
      {
         for( size_t jc=0 ; jc<nb_dims ; jc++ )
         {
            double x = c_w * dN_row( i, jc ) ;
            double y = x* dN_col( j, jc ) ;
            for( size_t ic=0 ; ic<nb_dims ; ic++ )
            {
               leq->add_to_matrix( x*dN_col( j, ic ), i, j, ic, jc ) ;
               leq->add_to_matrix( y, i, j, ic, ic) ;
            }
         }
      }
   }

   if( FE::geometry() == FE::axisymmetrical )
   {
      double radius = fe->coordinates_of_IP()->coordinate(0) ;
      doubleVector const& N_row = fe->Ns_at_IP( row ) ;
      doubleVector const& N_col = fe->Ns_at_IP( col ) ;
      double yy = 2. * c_w / radius / radius ;
      
      for( size_t i=0 ; i<nb_nodes ; ++i )
      {
         for( size_t j=i ; j<nb_nodes ; ++j )
         {
            double xx = yy * N_row( i ) * N_col( j ) ;
            leq->add_to_matrix( xx, i, j, 0, 0 ) ;
            if( i != j ) leq->add_to_matrix( xx, j, i, 0, 0 ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_row_vvgrad_col( PDE_LocalEquation* leq,
                         PDE_LocalFEcell const* fe,
                         doubleVector const& aa,
                         double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_row_vvgrad_col" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::N ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::col ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == 
                     fe->field( PDE_LocalFE::col )->nb_components() ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col )->nb_components() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( aa.size() == fe->nb_space_dimensions() ) ;
   
   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_row_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_col_nodes = fe->nb_basis_functions( col ) ;
   double c_w = fe->weight_of_IP() * coef ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   doubleArray2D const& dN_col = fe->dNs_at_IP( col ) ;
   doubleVector const& N_row = fe->Ns_at_IP( row ) ;
   
   for( size_t j=0 ; j<nb_col_nodes ; ++j )
   {
      double aa_grad = 0. ;
      for( size_t d=0 ; d<nb_dims ; ++d )
      {
         aa_grad += aa(d) * dN_col( j, d ) ;
      }
      for( size_t i=0 ; i<nb_row_nodes ; ++i )
      {
         double xx = c_w * N_row( i ) * aa_grad ;
         for( size_t ic=0 ; ic<nbc ; ++ic )
         {
            leq->add_to_matrix( xx, i, j, ic, ic ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_grad_row_vv_otimes_col( PDE_LocalEquation* leq,
                                 PDE_LocalFEcell const* fe,
                                 doubleVector const& aa,
                                 double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_grad_row_vv_otimes_col" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::col ), PDE_LocalFE::N ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                     fe->nb_basis_functions( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == 
                     fe->field( PDE_LocalFE::col )->nb_components() ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col )->nb_components() == 
                     fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( aa.size() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double c_w = fe->weight_of_IP() * coef ;
  
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   doubleVector const& N_col = fe->Ns_at_IP( col ) ;

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      for( size_t j=0 ; j<nb_nodes ; ++j )
      {
         for( size_t ic=0 ; ic<nbc ; ++ic )
         {
            for( size_t jc=0 ; jc<nbc ; ++jc )
            {
               double xx = c_w * aa( ic ) * N_col( j ) * dN_row( i, jc ) ;
               leq->add_to_matrix( xx, i, j, ic, jc ) ;
            }
         }
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_vvgrad_row_col( PDE_LocalEquation* leq,
                         PDE_LocalFEcell const* fe,
                         doubleVector const& aa,
                         double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_vvgrad_row_col" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::col ), PDE_LocalFE::N ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                     fe->nb_basis_functions( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == 
                     fe->field( PDE_LocalFE::col )->nb_components() ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col )->nb_components() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( aa.size() == fe->nb_space_dimensions() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_row_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_col_nodes = fe->nb_basis_functions( col ) ;
   double c_w = fe->weight_of_IP() * coef ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   doubleVector const& N_col = fe->Ns_at_IP( col ) ;

   for( size_t i=0 ; i<nb_row_nodes ; ++i )
   {
      double aagrad = 0.0 ;
      for( size_t d=0 ; d<nb_dims ; ++d )
      {
         aagrad +=  aa(d) * dN_row( i, d ) ;
      }
      aagrad *= c_w ;
      for( size_t j=0 ; j<nb_col_nodes ; ++j )
      {
         double xx = aagrad * N_col( j ) ;

         for( size_t ic=0 ; ic<nbc ; ++ic )
         {
            leq->add_to_matrix( xx, i, j, ic, ic ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_vvgrad_row_vvgrad_col( PDE_LocalEquation* leq,
                                PDE_LocalFEcell const* fe,
                                doubleVector const& aa,
                                double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_vvgrad_row_vvgrad_col" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) == 
                     fe->field( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == 
                     fe->field( PDE_LocalFE::col )->nb_components() ) ;
   PEL_CHECK_PRE( aa.size() == fe->nb_space_dimensions() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double c_w = fe->weight_of_IP() * coef ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   doubleVector aa_grad = doubleVector( nb_nodes ) ;
   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      for( size_t d=0 ; d<nb_dims ; ++d )
      {
         aa_grad(i) +=  aa(d) * dN_row( i, d ) ;
      }
   }
   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      for( size_t j=0 ; j<nb_nodes ; ++j )
      {
         double xx = c_w * aa_grad(i) * aa_grad(j) ;
         for( size_t ic=0 ; ic<nbc ; ++ic )
         {
            leq->add_to_matrix( xx, i, j, ic, ic ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_vvgrad_row_lapl_col( PDE_LocalEquation* leq,
                              PDE_LocalFEcell const* fe,
                              doubleVector const& aa,
                              double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_vvgrad_row_lapl_col" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian  ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) == 
                     fe->field( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::col ), PDE_LocalFE::d2N ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == 
                     fe->field( PDE_LocalFE::col )->nb_components() ) ;
   PEL_CHECK_PRE( aa.size() == fe->nb_space_dimensions() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double c_w = fe->weight_of_IP() * coef ;
   
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   doubleArray3D const& d2N_col = fe->d2Ns_at_IP( col ) ;

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      double aagrad = 0.0 ;
      for( size_t d=0 ; d<nb_dims ; ++d )
      {
         aagrad +=  aa(d) * dN_row( i, d ) ;
      }
      aagrad *= c_w ;
      for( size_t j=0 ; j<nb_nodes ; ++j )
      {
         double xx = 0.0 ; //??????????????????????? xx ne depend pas de i
         for( size_t d=0 ; d<nb_dims ; ++d )
         {
            xx += d2N_col( j, d, d ) ;
         }
         xx *= aagrad ;

         for( size_t ic=0 ; ic<nbc ; ++ic )
         {
            leq->add_to_matrix( xx, i, j, ic, ic ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_vvgrad_row_div_D_col( PDE_LocalEquation* leq,
                               PDE_LocalFEcell const* fe,
                               doubleVector const& aa,
                               double coef,
                               doubleVector const& dcoef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_vvgrad_row_div_D_col" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) == 
                     fe->field( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col )->nb_components()==
                     fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::col ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::col ), PDE_LocalFE::d2N ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                  fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                  fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                  fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                  fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( aa.size() == fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( dcoef.size() == fe->nb_space_dimensions() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double w = fe->weight_of_IP() ;

   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   doubleArray2D const& dN_col = fe->dNs_at_IP( col ) ;
   doubleArray3D const& d2N_col = fe->d2Ns_at_IP( col ) ;

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      double aagrad = 0.0 ;
      for( size_t d=0 ; d<nb_dims ; ++d )
      {
         aagrad +=  aa(d) * dN_row( i, d ) ;
      }
      aagrad *= w ;
      for( size_t j=0 ; j<nb_nodes ; ++j )
      {
         double divNj = 0.0 ;
         double gradNjgradmu = 0.0 ;
         for( size_t d=0 ; d<nb_dims ; ++d ) //???????? independence avec i
         {
            divNj += d2N_col( j, d, d ) ;
            gradNjgradmu += dN_col( j, d ) * dcoef( d ) ;
         }

         for( size_t ic=0 ; ic<nbc ; ++ic )
         {
            for( size_t jc=0 ; jc<nbc ; ++jc )
            {
               double xx = coef * d2N_col( j, ic, jc )
    		                + dcoef(jc) * dN_col( j, ic ) ;
               if( ic == jc )
               {
                  xx += coef*divNj + gradNjgradmu ;
               }
               xx *= aagrad ;
               leq->add_to_matrix( xx, i, j, ic, jc ) ;
            }
         }
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_vvgrad_row( PDE_LocalEquation* leq,
                     PDE_LocalFEcell const* fe,
                     doubleVector const& aa,
                     double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_vvgrad_row" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row )->nb_components() == 1 ) ;   
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( aa.size() == fe->nb_space_dimensions() ) ;

   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double c_w = fe->weight_of_IP() * coef ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      double xx = 0.0 ;
      for( size_t d=0 ; d<nb_dims ; ++d )
      {
         xx +=  aa(d) * dN_row( i, d ) ;
      }
      xx *= c_w ;
      leq->add_to_vector( xx, i ) ;
   }
}

//----------------------------------------------------------------------
void
FE:: add_vvgrad_row( PDE_LocalEquation* leq,
                     PDE_LocalFEcell const* fe,
                     doubleVector const& aa,
                     doubleVector const& coef_1,
                     double coef_2 )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_vvgrad_row" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ; 
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( aa.size() == fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( coef_1.size() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double w = coef_2 * fe->weight_of_IP() ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      double yy = 0.0 ;
      for( size_t d=0 ; d<nb_dims ; ++d )
      {
         yy +=  aa(d) * dN_row( i, d ) ;
      }
      yy *= w ;
      for( size_t ic=0 ; ic<nbc ; ++ic )
      {
         double xx = yy * coef_1( ic ) ;
         leq->add_to_vector( xx, i, ic ) ;
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_row_div_col( PDE_LocalEquation* leq,
                      PDE_LocalFEcell const* fe,
                      double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_row_div_col" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 
                  fe->field( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row )->nb_components() == 1 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col )->nb_components() ==
                  fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::N ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::col ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->field_calculation_is_handled(
                              fe->field( PDE_LocalFE::col ),
                              PDE_LocalFE::N ) ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                  fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                  fe->nb_basis_functions( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                  fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                  fe->field( PDE_LocalFE::col )->nb_components() ) ;

   size_t nb_row_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_col_nodes = fe->nb_basis_functions( col ) ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   double c_w = fe->weight_of_IP() * coef ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   doubleArray2D const& dN_col = fe->dNs_at_IP( col ) ;
   doubleVector const& N_row = fe->Ns_at_IP( row ) ;

   for( size_t i=0 ; i<nb_row_nodes ; ++i )
   {
      for( size_t j=0 ; j<nb_col_nodes ; ++j )
      {
         for( size_t jc=0 ; jc<nb_dims ; ++jc )
         {
            double xx = c_w * N_row( i ) * dN_col( j, jc ) ;
            leq->add_to_matrix( xx, i, j, 0, jc ) ; 
         }
      }
   }

   if( FE::geometry() == FE::axisymmetrical )
   {
      // Axisymmetrical : div(f) = df/dr+df/dz+f/r
      double radius = fe->coordinates_of_IP()->coordinate(0) ;
      double yy = c_w / radius ;
      doubleVector const& N_col = fe->Ns_at_IP( col ) ;
      for( size_t i=0 ; i<nb_row_nodes ; ++i )
      {
         for( size_t j=0 ; j<nb_col_nodes ; ++j )
         {
            double xx = yy * N_row( i ) * N_col( j ) ;
            leq->add_to_matrix( xx, i, j, 0, 0 ) ;
         }
      }
   }
}
//----------------------------------------------------------------------
void
FE:: add_row_grad_col( PDE_LocalEquation* leq,
                       PDE_LocalFEcell const* fe,
                       double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_row_grad_col" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 
                  fe->field( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col )->nb_components() == 1 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row )->nb_components() ==
                  fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::N ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::col ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                  fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                  fe->nb_basis_functions( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                  fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                  fe->field( PDE_LocalFE::col )->nb_components() ) ;

   size_t nb_row_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_col_nodes = fe->nb_basis_functions( col ) ;
   size_t nb_dims = fe->nb_space_dimensions() ;

   double c_w = fe->weight_of_IP() * coef ;

   doubleArray2D const& dN_col = fe->dNs_at_IP( col ) ;
   doubleVector const& N_row = fe->Ns_at_IP( row ) ;

   for( size_t i=0 ; i<nb_row_nodes ; ++i )
   {
      for( size_t j=0 ; j<nb_col_nodes ; ++j )
      {
         for( size_t ic=0 ; ic<nb_dims ; ++ic )
         {
            double xx = c_w * N_row( i ) * dN_col( j, ic ) ;
            leq->add_to_matrix( xx, i, j, ic, 0 ) ; 
         }
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_div_row_div_col( PDE_LocalEquation* leq,
                          PDE_LocalFEcell const* fe,
                          double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_div_row_div_col" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row )->nb_components() == 
                     fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col )->nb_components() ==
                     fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->field_calculation_is_handled(
                              fe->field( PDE_LocalFE::row ),
                              PDE_LocalFE::N ) ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::col ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->field_calculation_is_handled(
                              fe->field( PDE_LocalFE::col ),
                              PDE_LocalFE::N ) ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                     fe->nb_basis_functions( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                    fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                     fe->field( PDE_LocalFE::col )->nb_components() ) ;

   size_t nb_row_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_col_nodes = fe->nb_basis_functions( col ) ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   double c_w = fe->weight_of_IP() * coef ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   doubleArray2D const& dN_col = fe->dNs_at_IP( col ) ;
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;

   if( FE::geometry() == FE::cartesian )
   {
      // Cartesian : div(f) = df/dx+df/dy+df/dz
      for( size_t i=0 ; i<nb_row_nodes ; ++i )
      {
         for( size_t j=0 ; j<nb_col_nodes ; ++j )
         {
            for( size_t ic=0 ; ic<nb_dims ; ++ic )
            {
               for( size_t jc=0 ; jc<nb_dims ; ++jc )
               {
                  double xx = dN_row( i, ic ) * dN_col( j, jc ) * c_w ;
                  leq->add_to_matrix( xx, i, j, ic, jc ) ; 
               }
            }
         }
      }
   }
   else if( FE::geometry() == FE::axisymmetrical )
   {
      // Axisymmetrical : div(f) = df/dr+df/dz+f/r
      double radius = fe->coordinates_of_IP()->coordinate(0) ;
      doubleVector const& N_col = fe->Ns_at_IP( col ) ;
      doubleVector const& N_row = fe->Ns_at_IP( row ) ;
      for( size_t i=0 ; i<nb_row_nodes ; ++i )
      {
         for( size_t j=0 ; j<nb_col_nodes ; ++j )
         {
            for( size_t ic=0 ; ic<nb_dims ; ++ic )
            {
               double divi = dN_row( i, ic ) ;
               if( ic == 0 )
               {
                  divi += N_row( i ) / radius ;
               }
               for( size_t jc=0 ; jc<nb_dims ; ++jc )
               {
                  double divj = dN_col( j, jc ) ;
                  if( jc == 0 )
                  {
                     divj += N_col( j ) / radius  ;
                  }
                  double xx = divi * divj * c_w ;
                  leq->add_to_matrix( xx, i, j, ic, jc ) ; 
               }
            }
         }
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_row_vv_col( PDE_LocalEquation* leq,
                     PDE_LocalFEcell const* fe,
                     doubleVector const& aa,
                     double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_row_vv_col" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 
                  fe->field( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row )->nb_components() == 1 ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::N ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::col ), PDE_LocalFE::N ) ) ;
   PEL_CHECK_PRE( aa.size() ==
                  fe->field( PDE_LocalFE::col )->nb_components() ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                  fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                  fe->nb_basis_functions( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                  fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() ==
                  fe->field( PDE_LocalFE::col )->nb_components() ) ;
   size_t nb_row_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_col_nodes = fe->nb_basis_functions( col ) ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   double c_w = fe->weight_of_IP() * coef ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   doubleVector const& N_col = fe->Ns_at_IP( col ) ;
   doubleVector const& N_row = fe->Ns_at_IP( row ) ;

   for( size_t i=0 ; i<nb_row_nodes ; ++i )
   {
      for( size_t j=0 ; j<nb_col_nodes ; ++j )
      {
         for( size_t jc=0 ; jc<nb_dims ; ++jc )
         {
            double xx = c_w * N_row( i ) * aa( jc ) * N_col( j ) ;
            leq->add_to_matrix( xx, i, j, 0, jc ) ; 
         }
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_lapl_row_vvgrad_col( PDE_LocalEquation* leq,
                              PDE_LocalFEcell const* fe,
                              doubleVector const& aa,
                              double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE::add_lapl_row_vvgrad_col " ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) == 
                  fe->field( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::d2N ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == 
                     fe->field( PDE_LocalFE::col )->nb_components() ) ;
   PEL_CHECK_PRE( aa.size() == fe->nb_space_dimensions() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double c_w = fe->weight_of_IP() * coef ;
 
   doubleArray2D const& dN_col = fe->dNs_at_IP( col ) ;
   doubleArray3D const& d2N_row = fe->d2Ns_at_IP( row ) ;

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      double xx = 0.0 ;  //??????????????????????? xx ne depend pas de i
      for( size_t d=0 ; d<nb_dims ; ++d )
      {
         xx +=  d2N_row( i, d, d ) ;
      }
      xx *= c_w  ; 
      for( size_t j=0 ; j<nb_nodes ; ++j )
      {
         double aagrad = 0.0 ;
         for( size_t d=0 ; d<nb_dims ; ++d )
         {
            aagrad +=  aa( d ) * dN_col( j, d ) ;
         }
         aagrad *= xx ;	 

         for( size_t ic=0 ; ic<nbc ; ++ic )
         {
            leq->add_to_matrix( aagrad, i, j, ic, ic ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_row_lapl_col( PDE_LocalEquation* leq,
                       PDE_LocalFEcell const* fe,
                       double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE::add_row_lapl_col " ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) == 
                     fe->field( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::N ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::col ), PDE_LocalFE::d2N ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                     fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                     fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == 
                     fe->field( PDE_LocalFE::col )->nb_components() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double c_w = fe->weight_of_IP() * coef ;
  
   doubleArray3D const& d2N_col = fe->d2Ns_at_IP( col ) ;
   doubleVector const& N_row = fe->Ns_at_IP( row ) ;

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      double xx = 0.0 ; 
      for( size_t d=0 ; d<nb_dims ; ++d )
      {
         xx += d2N_col( i, d, d ) ;
      }
      xx *= c_w ;
      for( size_t j=0 ; j<nb_nodes ; ++j )
      {
         double yy = xx * N_row( j )  ; 

         for( size_t ic=0 ; ic<nbc ; ++ic )
         {
            leq->add_to_matrix( yy, j, i, ic, ic ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_lapl_row_lapl_col( PDE_LocalEquation* leq,
                            PDE_LocalFEcell const* fe,
                            double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE::add_lapl_row_lapl_col " ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) == 
                  fe->field( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::d2N ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::col ), PDE_LocalFE::d2N ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                  fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                  fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                  fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == 
                  fe->field( PDE_LocalFE::col )->nb_components() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double c_w = fe->weight_of_IP() * coef ;
   
   doubleArray3D const& d2N_col = fe->d2Ns_at_IP( col ) ;
   doubleArray3D const& d2N_row = fe->d2Ns_at_IP( row ) ;

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      double xx = 0.0 ; 
      for( size_t d=0 ; d<nb_dims ; ++d )
      {
         xx += d2N_row( i, d, d ) ;
      }
      xx *= c_w ;
      for( size_t j=0 ; j<nb_nodes ; ++j )
      {
         double yy = 0. ; 
         for( size_t d=0 ; d<nb_dims ; ++d )
         {
            yy += d2N_col( j, d, d ) ;
         }
         yy *= xx ;
         for( size_t ic=0 ; ic<nbc ; ++ic )
         {
            leq->add_to_matrix( yy, j, i, ic, ic ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_lapl_row_col( PDE_LocalEquation* leq,
                       PDE_LocalFEcell const* fe,
                       double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE::add_lapl_row_col " ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) == 
                  fe->field( PDE_LocalFE::col ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::N ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::col ), PDE_LocalFE::d2N ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                  fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_columns() == 
                  fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                  fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( leq->nb_column_sub_indices() == 
                  fe->field( PDE_LocalFE::col )->nb_components() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double c_w = fe->weight_of_IP() * coef ;
  
   doubleArray3D const& d2N_col = fe->d2Ns_at_IP( col ) ;
   doubleVector const& N_row = fe->Ns_at_IP( row ) ;

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      double xx = 0.0 ; 
      for( size_t d=0 ; d<nb_dims ; ++d )
      {
         xx += d2N_col( i, d, d ) ;
      }
      xx *= c_w  ; 
      for( size_t j=0 ; j<nb_nodes ; ++j )
      {
         double yy = N_row( j ) * xx ;	 

         for( size_t ic=0 ; ic<nbc ; ++ic )
         {
            leq->add_to_matrix( yy, j, i, ic, ic ) ;
         }
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_lapl_row( PDE_LocalEquation* leq,
                   PDE_LocalFE const* fe,
                   double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_lapl_row" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row )->nb_components() == 1 ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::col ), PDE_LocalFE::d2N ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                  fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                  fe->field( PDE_LocalFE::row )->nb_components() ) ;

   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   double c_w = fe->weight_of_IP() * coef ;
 
   doubleArray3D const& d2N_col = fe->d2Ns_at_IP( col ) ;

   for( size_t i=0 ; i<nb_nodes ; i++ )
   {
      double xx = 0.0 ; 
      for( size_t d=0 ; d<nb_dims ; ++d )
      {
         xx += d2N_col( i, d, d ) ;
      }
      xx *= c_w  ; 
      leq->add_to_vector( xx, i ) ;
   }
}

//----------------------------------------------------------------------
void
FE:: add_lapl_row( PDE_LocalEquation* leq,
                   PDE_LocalFE const* fe,
                   doubleVector const& coef_1,
                   double coef_2 )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_lapl_row" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::col ), PDE_LocalFE::d2N ) ) ;
   PEL_CHECK_PRE( leq != 0 ) ;
   PEL_CHECK_PRE( leq->nb_rows() == 
                  fe->nb_basis_functions( PDE_LocalFE::row ) ) ;
   PEL_CHECK_PRE( leq->nb_row_sub_indices() == 
                  fe->field( PDE_LocalFE::row )->nb_components() ) ;
   PEL_CHECK_PRE( coef_1.size() == 
                  fe->field( PDE_LocalFE::row )->nb_components() ) ;

   size_t nbc = fe->field( row )->nb_components() ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   size_t nb_nodes = fe->nb_basis_functions( row ) ;
   double w = coef_2 * fe->weight_of_IP() ;
   
   doubleArray3D const& d2N_col = fe->d2Ns_at_IP( col ) ;

   for( size_t i=0 ; i<nb_nodes ; ++i )
   {
      double yy = 0.0 ; 
      for( size_t d=0 ; d<nb_dims ; ++d )
      {
         yy += d2N_col( i, d, d ) ;
      }
      yy *= w  ; 

      for( size_t ic=0 ; ic<nbc ; ++ic )
      {
         double xx = yy * coef_1( ic ) ;
         leq->add_to_vector( xx, i, ic ) ;
      }
   }
}

//----------------------------------------------------------------------
void
FE:: add_D_u_D_u_at_IP( double& target,
                        PDE_DiscreteField const* u,
                        size_t u_level,
                        PDE_LocalFE const* fe,
                        double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_D_u_D_u_at_IP" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( u != 0 ) ;
   PEL_CHECK_PRE( u->nb_components() == fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( u->storage_depth() > u_level ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled( u, PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( IMPLIES(
                     FE::geometry() == FE::axisymmetrical,
                     fe->field_calculation_is_handled(
                        u, PDE_LocalFE::N ) ) ) ;

   size_t const nb_comps = u->nb_components() ;
   size_t const nb_dims = fe->nb_space_dimensions() ;

   doubleArray2D grad( nb_comps, nb_dims ) ;
   for( size_t i=0 ; i<nb_dims ; ++i )
   {
      for( size_t j=0 ; j<nb_comps ; ++j )
      {
         grad(i,j) = fe->gradient_at_IP( u, u_level, i, j ) ;
      }
   }

   double xx = 0. ;
   for( size_t i=0 ; i<nb_dims ; ++i )
   {
      for( size_t j=i+1 ; j<nb_comps ; ++j )
      {
         double dij = grad(i,j)+grad(j,i) ;
         xx += 2.*dij*dij ;
      }
      double dii = 2.*grad(i,i) ;
      xx += dii*dii ;
   }
   
   if( FE::geometry() == FE::axisymmetrical )
   {
      double radius = fe->coordinates_of_IP()->coordinate(0) ;
      double dii = 2.*fe->value_at_IP( u, u_level, 0 )/radius ;
      xx += dii*dii ;
   }

   target += coef*xx ;
}

//----------------------------------------------------------------------
void
FE:: add_D_u_D_u_at_pt( double& target,
                        PDE_DiscreteField const* u,
                        size_t u_level,
                        PDE_LocalFE const* fe,
                        double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_D_u_D_u_at_pt" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->calculation_point() != 0 ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->calculation_point()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( u != 0 ) ;
   PEL_CHECK_PRE( u->nb_components() == fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( u->storage_depth() > u_level ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled( u, PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( IMPLIES(
                     FE::geometry() == FE::axisymmetrical,
                     fe->field_calculation_is_handled(
                        u, PDE_LocalFE::N ) ) ) ;

   size_t const nb_comps = u->nb_components() ;
   size_t const nb_dims = fe->nb_space_dimensions() ;

   doubleArray2D grad( nb_comps, nb_dims ) ;
   for( size_t i=0 ; i<nb_dims ; ++i )
   {
      for( size_t j=0 ; j<nb_comps ; ++j )
      {
         grad(i,j) = fe->gradient_at_pt( u, u_level, i, j ) ;
      }
   }

   double xx = 0. ;
   for( size_t i=0 ; i<nb_dims ; ++i )
   {
      for( size_t j=i+1 ; j<nb_comps ; ++j )
      {
         double dij = grad(i,j)+grad(j,i) ;
         xx += 2.*dij*dij ;
      }
      double dii = 2.*grad(i,i) ;
      xx += dii*dii ;
   }
   
   if( FE::geometry() == FE::axisymmetrical )
   {
      double radius = fe->calculation_point()->coordinate(0) ;
      double dii = 2.*fe->value_at_pt( u, u_level, 0 )/radius ;
      xx += dii*dii ;
   }

   target += coef*xx ;
}

//----------------------------------------------------------------------
void
FE:: add_div_u_at_IP( double& target,
                      PDE_DiscreteField const* u,
                      size_t u_level,
                      PDE_LocalFE const* fe,
                      double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_div_u_at_IP" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( u != 0 ) ;
   PEL_CHECK_PRE( u->nb_components() == fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( u->storage_depth() > u_level ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled( u, PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( IMPLIES(
                     FE::geometry() == FE::axisymmetrical,
                     fe->field_calculation_is_handled(
                        u, PDE_LocalFE::N ) ) ) ;

   size_t const nb_dims = fe->nb_space_dimensions() ;

   double xx = 0. ;
   for( size_t i=0 ; i<nb_dims ; ++i )
   {
      xx += fe->gradient_at_IP( u, u_level, i, i ) ;
   }
   
   if( FE::geometry() == FE::axisymmetrical )
   {
      double radius = fe->coordinates_of_IP()->coordinate(0) ;
      xx += fe->value_at_IP( u, u_level, 0 )/radius ;
   }

   target += coef*xx ;
}

//----------------------------------------------------------------------
void
FE:: add_div_u_at_pt( double& target,
                      PDE_DiscreteField const* u,
                      size_t u_level,
                      PDE_LocalFE const* fe,
                      double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE:: add_div_u_at_pt" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->calculation_point() != 0 ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->calculation_point()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( u != 0 ) ;
   PEL_CHECK_PRE( u->nb_components() == fe->nb_space_dimensions() ) ;
   PEL_CHECK_PRE( u->storage_depth() > u_level ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled( u, PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( IMPLIES(
                     FE::geometry() == FE::axisymmetrical,
                     fe->field_calculation_is_handled(
                        u, PDE_LocalFE::N ) ) ) ;

   size_t const nb_dims = fe->nb_space_dimensions() ;

   double xx = 0. ;
   for( size_t i=0 ; i<nb_dims ; ++i )
   {
      xx += fe->gradient_at_pt( u, u_level, i, i ) ;
   }
   
   if( FE::geometry() == FE::axisymmetrical )
   {
      double radius = fe->calculation_point()->coordinate(0) ;
      xx += fe->value_at_pt( u, u_level, 0 )/radius ;
   }

   target += coef*xx ;
}
