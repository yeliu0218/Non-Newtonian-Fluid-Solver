#include <FE_TensorFormAssembling.hh>
#include <UT_Viscoplastic.hh>
#include <AS_Viscoplastic.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Communicator.hh>
#include <PEL_Exec.hh>

#include <GE_Point.hh>
#include <GE_Mpolyhedron.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_ReferenceElement.hh>

#include <FE.hh>

#include <doubleVector.hh>

PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
PDE_LocalFE::field_id const col = PDE_LocalFE::col ;


/* Assemble the elementary matrix at the level of the element
   for the term div(lambda) where lambda is a tensor stored as a vector
-----------------------------------------------------------------------*/
void FE_TensorFormAssembling::add_grad_row_col( PDE_LocalEquation* leq,
		PDE_LocalFEcell const* fe, double coef )
{
   PEL_LABEL( "FE:: add_grad_row_col_2D" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;

   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;

   PEL_CHECK_PRE( fe->nb_space_dimensions() == 2 ||
   	fe->nb_space_dimensions() == 3 ) ;
   PEL_CHECK_PRE( IMPLIES( fe->field( PDE_LocalFE::row )->nb_components() == 2,
   	fe->field( PDE_LocalFE::col )->nb_components() == 3 ) ) ;
   PEL_CHECK_PRE( IMPLIES( fe->field( PDE_LocalFE::row )->nb_components() == 3,
   	fe->field( PDE_LocalFE::col )->nb_components() == 6 ) ) ;

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

   size_t nbc_row = fe->field(row)->nb_components() ;
   size_t nb_row_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_col_nodes = fe->nb_basis_functions( col ) ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   double c_w = fe->weight_of_IP() * coef ;

   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   doubleVector const& N_col = fe->Ns_at_IP( col ) ;

   for( size_t i=0 ; i<nb_row_nodes ; ++i )
   {
      for( size_t j=0 ; j<nb_col_nodes ; ++j )
      {
         for( size_t ir=0 ; ir<nbc_row ; ++ir )
         {
        	 for( size_t ic=0 ; ic<nbc_row ; ++ic )
        	 {
        		 double xx = dN_row( i, ic ) * N_col(j) ;
        		 xx *= c_w ;
        		 size_t k=UT_Viscoplastic::strain_rate_component_number(
								 nb_dims,ir,ic).first;
        		 leq->add_to_matrix( xx, i, j, ir, k ) ;
        	 }
         }
      }
   }
}

/* Assemble the elementary matrix at the level of the element
   for the term div(lambda) where lambda is a tensor stored as a vector
-----------------------------------------------------------------------*/
void FE_TensorFormAssembling::add_row_col( PDE_LocalEquation* leq,
		PDE_LocalFEcell const* fe, double coef )
{
   PEL_LABEL( "FE_TensorFormAssembling:: add_row_col" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;

   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;

   PEL_CHECK_PRE( fe->nb_space_dimensions() == 2 ||
   	fe->nb_space_dimensions() == 3 ) ;
   PEL_CHECK_PRE( IMPLIES( fe->field( PDE_LocalFE::row )->nb_components() == 2,
   	fe->field( PDE_LocalFE::col )->nb_components() == 3 ) ) ;
   PEL_CHECK_PRE( IMPLIES( fe->field( PDE_LocalFE::row )->nb_components() == 3,
   	fe->field( PDE_LocalFE::col )->nb_components() == 6 ) ) ;

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

   size_t nbc_row = fe->field(row)->nb_components() ;
   size_t nbc_col = fe->field(col)->nb_components() ;

   size_t nb_row_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_col_nodes = fe->nb_basis_functions( col ) ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   double c_w = fe->weight_of_IP() * coef ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   doubleVector const& N_row = fe->Ns_at_IP( row ) ;
   doubleVector const& N_col = fe->Ns_at_IP( col ) ;

//   PEL::out() << "nbc_row = " << nbc_row << " nbc_col = " << nbc_col << std::endl;
//   PEL::out() << "nb_row_nodes = " << nb_row_nodes << " nb_col_nodes = " << nb_col_nodes << std::endl;
//   PEL::out() << "nb_dims = " << nb_dims << std::endl;
//   nbc_row = 3 nbc_col = 3
//   nb_row_nodes = 4 nb_col_nodes = 4
//   nb_dims = 2

   for( size_t i=0 ; i<nb_row_nodes ; ++i )
   {
      for( size_t j=0 ; j<nb_col_nodes ; ++j )
      {
          for( size_t ir=0 ; ir<nbc_row ; ++ir )
          {
        	  for( size_t ic=0 ; ic<nbc_col; ++ic )
        	  {
        		  double xx = N_row( i ) * N_col(j) ;
        		  xx *= c_w ;
        		  size_t k=UT_Viscoplastic::strain_rate_component_number(
						  nb_dims,ir,ic).first;
        		  leq->add_to_matrix( xx, i, j, ir, k ) ;
        	  }
          }
      }
   }
}



/* Assemble the elementary matrix at the level of the element
   for the lumped mass term
-------------------------------------------------------------*/
void FE_TensorFormAssembling:: add_lumped_row_col( PDE_LocalEquation* leq,
   	  	         PDE_LocalFEcell const* fe,
  		         double coef )
{
   PEL_LABEL( "MY_NavierStokes:: add_lumped_row_col" ) ;
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
   double measure = fe->polyhedron()->measure() ;
   double node_contribution = coef * measure ;
   if( FE::geometry() == FE::axisymmetrical )
   {
     node_contribution *= fe->polyhedron()->center()->coordinate(0) ;
   }

   if (fe->reference_element(fe->field(row))->name()=="PDE_2D_P1isoP2_6nodes")
   {
     node_contribution /= 3.;
     for( size_t i=0 ; i<nb_nodes ; ++i )
     {
       double xx = node_contribution;
       if ((i==0)||(i==3)||(i==5)) xx *= 0.25;
       else xx *= 0.75;
       for( size_t ic=0 ; ic<nbc ; ++ic )
	 leq->add_to_matrix( xx, i, i, ic, ic ) ;
     }
   }
   else if (fe->reference_element(fe->field(row))->name()=="PDE_2D_P2_6nodes")
   {
     double xx = node_contribution / 6.;
     for( size_t i=0 ; i<nb_nodes ; ++i )
       for( size_t ic=0 ; ic<nbc ; ++ic )
	 leq->add_to_matrix( xx, i, i, ic, ic ) ;
   }
   else if (fe->reference_element(fe->field(row))->name()
   	=="PDE_2D_Q1isoNonConfA_4nodes")
   {
     double xx = node_contribution / 4.;
     for( size_t i=0 ; i<nb_nodes ; ++i )
       for( size_t ic=0 ; ic<nbc ; ++ic )
	 leq->add_to_matrix( xx, i, i, ic, ic ) ;
   }
   else if (fe->reference_element(fe->field(row))->name()=="PDE_3D_P2_10nodes")
   {
     double xx = node_contribution / 10.;
     for( size_t i=0 ; i<nb_nodes ; ++i )
       for( size_t ic=0 ; ic<nbc ; ++ic )
	 leq->add_to_matrix( xx, i, i, ic, ic ) ;
   }
   else
   {
     cout << "Mass lumping formula not implemented yet for element "
     	<< fe->reference_element(fe->field(row))->name() << endl;
     PEL_Error::exit();
   }

}

//by Leo

void FE_TensorFormAssembling::add_grad_row_col_A( PDE_LocalEquation* leq,
		PDE_LocalFEcell const* fe, double coef )
{
   PEL_LABEL( "FE:: add_grad_row_col_A" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;

   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;

   PEL_CHECK_PRE( fe->nb_space_dimensions() == 2 ||
   	fe->nb_space_dimensions() == 3 ) ;
   PEL_CHECK_PRE( IMPLIES( fe->field( PDE_LocalFE::row )->nb_components() == 2,
   	fe->field( PDE_LocalFE::col )->nb_components() == 4 ) ) ;
   PEL_CHECK_PRE( IMPLIES( fe->field( PDE_LocalFE::row )->nb_components() == 3,
   	fe->field( PDE_LocalFE::col )->nb_components() == 6 ) ) ;

   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::dN ) ) ;
   PEL_CHECK_PRE( fe->field_calculation_is_handled(
                     fe->field( PDE_LocalFE::row ), PDE_LocalFE::N ) ) ;
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

   size_t nbc_row = fe->field(row)->nb_components() ;
   size_t nb_row_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_col_nodes = fe->nb_basis_functions( col ) ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   double c_w = fe->weight_of_IP() * coef ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= fe->coordinates_of_IP()->coordinate(0) ;
   }
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   doubleVector const& N_col = fe->Ns_at_IP( col ) ;
   doubleVector const& N_row = fe->Ns_at_IP( row ) ;

   for( size_t i=0 ; i<nb_row_nodes ; ++i )
   {
      for( size_t j=0 ; j<nb_col_nodes ; ++j )
      {
         for( size_t ir=0 ; ir<nbc_row ; ++ir )
         {
         	 for( size_t ic=0 ; ic<nbc_row ; ++ic )
        	 {
        		 double xx = dN_row( i, ic ) * N_col(j) ;
        		 xx *= c_w ;
        		 size_t k=AS_Viscoplastic::strain_rate_component_number(
								 nb_dims,ir,ic).first;
        		 leq->add_to_matrix( xx, i, j, ir, k ) ;
        	 }
         }
         double yy = N_row(i)* N_col(j)*coef*fe->weight_of_IP();
         //double yy = N_row(i)* N_col(j)*c_w/fe->coordinates_of_IP()->coordinate(0); 
         leq->add_to_matrix( yy, i, j, 0, 3 ) ;
      }
   }
}

void FE_TensorFormAssembling::add_grad_row_col_Leo( PDE_LocalEquation* leq,
		PDE_LocalFEcell const* fe, double coef )
{
   PEL_LABEL( "FE:: add_grad_row_col_2D" ) ;
   PEL_CHECK_PRE( FE::geometry() == FE::cartesian ||
                  FE::geometry() == FE::axisymmetrical) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( fe->valid_IP() ) ;

   PEL_CHECK_PRE( fe->field( PDE_LocalFE::row ) != 0 ) ;
   PEL_CHECK_PRE( fe->field( PDE_LocalFE::col ) != 0 ) ;

   PEL_CHECK_PRE( fe->nb_space_dimensions() == 2 ||
   	fe->nb_space_dimensions() == 3 ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->nb_space_dimensions() == 2 ) ) ;
   PEL_CHECK_PRE( IMPLIES( FE::geometry() == FE::axisymmetrical,
                           fe->coordinates_of_IP()->coordinate(0) >= 0. ) ) ;
   PEL_CHECK_PRE( IMPLIES( fe->field( PDE_LocalFE::row )->nb_components() == 2,
   	fe->field( PDE_LocalFE::col )->nb_components() == 3 ) ) ;
   PEL_CHECK_PRE( IMPLIES( fe->field( PDE_LocalFE::row )->nb_components() == 3,
   	fe->field( PDE_LocalFE::col )->nb_components() == 6 ) ) ;

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

   size_t nbc_row = fe->field(row)->nb_components() ;
   size_t nb_row_nodes = fe->nb_basis_functions( row ) ;
   size_t nb_col_nodes = fe->nb_basis_functions( col ) ;
   size_t nb_dims = fe->nb_space_dimensions() ;
   double c_w = fe->weight_of_IP() * coef ;
   if( FE::geometry() == FE::axisymmetrical )
   {
      c_w *= (fe->coordinates_of_IP()->coordinate(0)) ;
   }
   doubleArray2D const& dN_row = fe->dNs_at_IP( row ) ;
   doubleVector const& N_col = fe->Ns_at_IP( col ) ;
   doubleVector const& N_row = fe->Ns_at_IP( row ) ;

   for( size_t i=0 ; i<nb_row_nodes ; ++i )
   {
      for( size_t j=0 ; j<nb_col_nodes ; ++j )
      {
         for( size_t ir=0 ; ir<nbc_row ; ++ir )
         {
        	 for( size_t ic=0 ; ic<nbc_row ; ++ic )
        	 {
        		 double xx = dN_row( i, ic ) * N_col(j) ;
        		 xx *= c_w ;
        		 size_t k=UT_Viscoplastic::strain_rate_component_number(
								 nb_dims,ir,ic).first;
        		 leq->add_to_matrix( xx, i, j, ir, k ) ;
        	 }
         }
         double yy = -N_row(i)*N_col(j)*fe->weight_of_IP()*coef;         
         leq->add_to_matrix( yy, i, j, 0, 0 ) ;
         leq->add_to_matrix( yy, i, j, 0, 1 ) ;
      }
   }
}
