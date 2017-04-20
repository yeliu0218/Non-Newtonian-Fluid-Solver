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

#include <PDE_BFvalues_TEST.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_VectorIterator.hh>
#include <doubleVector.hh>
#include <stringVector.hh>

#include <GE_Matrix.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QuadratureRule.hh>

#include <PDE_CellFE.hh>
#include <PDE_BFvalues.hh>
#include <PDE_DomainBuilder.hh>
#include <PDE_GridFE.hh>
#include <PDE_ReferenceElement.hh>

#include <ios>
#include <iomanip>
#include <iostream>
#include <limits>

using std::cout ;
using std::endl ;
using std::ios_base ;
using std::string ;

//---------------------------------------------------------------------------
PDE_BFvalues_TEST*
PDE_BFvalues_TEST:: registered_test = new PDE_BFvalues_TEST() ;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
PDE_BFvalues_TEST:: PDE_BFvalues_TEST( void )
//---------------------------------------------------------------------------
   : PEL_ObjectTest( "PDE_BFvalues", "PDE_BFvalues_TEST" )
{
}

//---------------------------------------------------------------------------
PDE_BFvalues_TEST:: ~PDE_BFvalues_TEST( void )
//---------------------------------------------------------------------------
{
}

// la regle de quadrature permet d'avoir un ensemble de points sur la maille
// de reference
//---------------------------------------------------------------------------
void
PDE_BFvalues_TEST:: process_one_test( PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_BFvalues_TEST:: process_one_test" ) ;

   bool eq ;
   double diff, val ;
   double dx = exp->double_data( "delta_x" ) ;
   double dbl_eps = exp->double_data( "dbl_epsilon" ) ;
   double dbl_min = exp->double_data( "dbl_minimum" ) ;

   stringVector const& elms = exp->stringVector_data( "elements" ) ;
   stringVector const& qrs  = exp->stringVector_data( "quadrature_rules" ) ;
   PEL_ASSERT( elms.size() == qrs.size() ) ;

   PDE_DomainBuilder* domain = PDE_DomainBuilder::create( 0, exp, "me" );
   PDE_GridFE const* grid = domain->finite_element_grid() ;
   PEL_VectorIterator* it = PEL_VectorIterator::create( domain,
                                                        grid->cells() ) ;

   PDE_BFvalues* bf     = PDE_BFvalues::create( domain ) ;
   PDE_BFvalues* bf_eps = PDE_BFvalues::create( domain ) ;

   size_t dim = grid->nb_space_dimensions() ;
   GE_Point* pt = GE_Point::create( domain, dim ) ;
   GE_Point* pt_eps = GE_Point::create( domain, dim ) ;
   GE_Point* pt_ref_eps = GE_Point::create( domain, dim ) ;
   GE_Matrix* tr_jac = GE_Matrix::create( domain, dim, dim ) ;
   GE_Matrix* tr_jac_eps = GE_Matrix::create( domain, dim, dim ) ;
   doubleArray3D* hessian = new doubleArray3D( dim, dim, dim ) ;

   int order = PDE_BFvalues::N ;
   order |= PDE_BFvalues::dN ;
   order |= PDE_BFvalues::d2N ;

   for( size_t e=0 ; e<elms.size() ; ++e )
   {
      PDE_ReferenceElement const* elm = 
                                   PDE_ReferenceElement::object( elms( e ) ) ;
      GE_QuadratureRule const* qr = GE_QuadratureRule::object( qrs( e ) ) ;

      bool ok_tr_jac  = true ;
      bool ok_mea     = true ;
      bool ok_dN      = true ;
      bool ok_d2N     = true ;

      for( it->start() ; it->is_valid() ; it->go_next() )
      {
	 GE_Mpolyhedron const* poly = 
                        static_cast<PDE_CellFE*>( it->item() )->polyhedron() ;

	 double itg_mea = 0.0 ;
	 for( size_t i=0 ; i<qr->nb_points() ; ++i )
	 {
	    GE_Point const* pt_ref = qr->point( i ) ;
	    poly->apply_mapping( pt_ref, pt ) ;

	    poly->build_tr_mapping_derivative( pt_ref, tr_jac ) ;
	    tr_jac->compute_determinant() ;

	    doubleArray3D const* hess = 0 ;
	    bool nonzero = false ;
	    poly->build_mapping_hessian( pt_ref, hessian, nonzero ) ;
	    if( nonzero ) hess = hessian ;

	    bf->re_initialize( elm, order, pt_ref, tr_jac, hess ) ;

	    for( size_t d=0 ; d<dim ; ++d )
	    {
	       doubleVector xx( pt_ref->coordinate_vector() ) ;
	       xx( d ) += dx ; // !!! maille de reference
	       pt_ref_eps->set_coordinates( xx ) ;
	       poly->apply_mapping( pt_ref_eps, pt_eps ) ;
	       for( size_t j=0 ; j<dim ; ++j )
	       {
		  diff = ( pt_eps->coordinate(j) - pt->coordinate(j) ) / dx ;
		  val = tr_jac->item( d, j ) ;
		  eq = PEL::double_equality( val, diff, dbl_eps, dbl_min ) ;
		  if( !eq ) display_error( "tr_jac", val, diff ) ;
		  ok_tr_jac = ok_tr_jac && eq ;
	       }
	    }

	    for( size_t d=0 ; d<dim ; ++d )
	    {
	       doubleVector xx( pt->coordinate_vector() ) ;
	       xx( d ) += dx ; // !!! maille reelle
	       pt_eps->set_coordinates( xx ) ;
	       poly->apply_inverse_mapping( pt_eps, pt_ref_eps ) ;

	       poly->build_tr_mapping_derivative( pt_ref_eps, tr_jac_eps ) ;
	       tr_jac_eps->compute_determinant() ;
	       bf_eps->re_initialize( elm, order, pt_ref_eps, tr_jac_eps, 0 ) ;

	       for( size_t n=0 ; n<elm->nb_nodes() ; ++n )
	       {
		  diff = ( bf_eps->N_at_pt( n ) - bf->N_at_pt( n ) ) / dx ;
		  val  = bf->dN_at_pt( n, d ) ;
		  eq = PEL::double_equality( val, diff, dbl_eps, dbl_min ) ;
		  if( !eq ) display_error( "dN", val, diff ) ;
		  ok_dN = ok_dN && eq ;
		  for( size_t d2=0 ; d2<dim ; ++d2 )
		  {
		     diff = ( bf_eps->dN_at_pt( n, d2 ) - 
			      bf->dN_at_pt( n, d2 ) ) / dx ;
		     val = bf->d2N_at_pt( n, d2, d ) ;
		     eq = PEL::double_equality( val, diff, dbl_eps, dbl_min  ) ;
		     if( !eq ) display_error( "d2N", val, diff ) ;
		     ok_d2N = ok_d2N && eq ;
		  }
	       }
	    }
	    itg_mea += qr->weight(i) * tr_jac->determinant() ;
	 }
	 double poly_mea = poly->measure() ;
	 eq = PEL::double_equality( itg_mea, poly_mea, 1.e-6, 1.e-50 ) ;
	 if( !eq ) display_error( "measure", itg_mea, poly_mea ) ;
	 ok_mea = ok_mea && eq ;
      }
      notify_one_test_result( "tr_jac", ok_tr_jac ) ;
      notify_one_test_result( "measure", ok_mea ) ;
      notify_one_test_result( "dN  diff", ok_dN ) ;
      notify_one_test_result( "d2N diff", ok_d2N ) ;

   }

   delete hessian ;
   domain->destroy() ;
}

//-----------------------------------------------------------------------
void 
PDE_BFvalues_TEST:: display_error( std::string const& mesg, 
                                   double val, double diff ) const
//-----------------------------------------------------------------------
{
   ios_base::fmtflags original_flags = out().flags() ;
   out().setf( ios_base::uppercase | ios_base::scientific ) ;
   out() << std::setprecision( 10 ) ;

   out() << mesg
         << std::setw( 20 ) << val 
         << std::setw( 20 ) << diff
         << std::setw( 20 ) << val-diff
         << endl ;

   out().flags( original_flags ) ;
}

