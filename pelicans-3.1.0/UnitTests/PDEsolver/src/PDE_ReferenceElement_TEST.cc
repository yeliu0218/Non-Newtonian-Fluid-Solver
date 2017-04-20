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

#include <PDE_ReferenceElement_TEST.hh>

#include <doubleArray2D.hh>
#include <stringVector.hh>

#include <GE_Point.hh>
#include <GE_ReferencePolyhedron.hh>

#include <LA_DenseMatrix.hh>
#include <LA_GaussLU_DS.hh>
#include <LA_SeqVector.hh>

#include <PEL.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Randomizer.hh>
#include <PEL_Variable.hh>

#include <PDE_ReferenceElement.hh>

#include <iostream>

//---------------------------------------------------------------------------
PDE_ReferenceElement_TEST*
PDE_ReferenceElement_TEST:: REGISTRATOR = new PDE_ReferenceElement_TEST() ;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
PDE_ReferenceElement_TEST:: PDE_ReferenceElement_TEST( void )
//---------------------------------------------------------------------------
      : PEL_ObjectTest( "PDE_ReferenceElement", "PDE_ReferenceElement_TEST" )
      , RANDOMIZER( PEL_Randomizer::create( this, 1 ) )
        
{
}

//---------------------------------------------------------------------------
PDE_ReferenceElement_TEST:: ~PDE_ReferenceElement_TEST( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
PDE_ReferenceElement_TEST:: process_one_test( PEL_ModuleExplorer const* texp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElement_TEST:: process_one_test" ) ;
   PEL_CHECK( texp!=0 ) ;

   std::string const& tname = texp->string_data( "name" ) ;
   PDE_ReferenceElement const* ref = PDE_ReferenceElement::object(tname) ;

   out() << "| ... " << tname << " : " << std::endl ;

   // Geometrical test :
   {
      PEL_ModuleExplorer const* e = texp->create_subexplorer( 0, "geometry" ) ;
      process_geometry( tname, ref, e ) ;
      e->destroy() ; e = 0 ;
   }

   // BF value at nodes :
   if( texp->has_module( "bf_value_at_nodes" ) )
   {
      PEL_ModuleExplorer const* e =
                texp->create_subexplorer( 0, "bf_value_at_nodes" ) ;
      process_bf_value_at_nodes( ref, e ) ;
      e->destroy() ; e = 0 ;
   }
   
   // BF value at random points :
   if( texp->has_module( "bf_value_at_random_points" ) )
   {
      PEL_ModuleExplorer const* e =
                texp->create_subexplorer( 0, "bf_value_at_random_points" ) ;
      process_bf_value_at_random_points( ref, e ) ;
      e->destroy() ; e = 0 ;
   }
   
   // dBF value at random points :
   if( texp->has_module( "dbf_value_at_random_points" ) )
   {
      PEL_ModuleExplorer const* e =
                texp->create_subexplorer( 0, "dbf_value_at_random_points" ) ;
      process_dbf_value_at_random_points( ref, e ) ;
      e->destroy() ; e = 0 ;
   }   
 
   // d2BF value at random points :
   if( texp->has_module( "d2bf_value_at_random_points" ) )
   {
      PEL_ModuleExplorer const* e =
                texp->create_subexplorer( 0, "d2bf_value_at_random_points" ) ;
      process_d2bf_value_at_random_points( ref, e ) ;
      e->destroy() ; e = 0 ;
   }

   // Value at random points :
   if( texp->has_module( "value_at_random_points" ) )
   {
      PEL_ModuleExplorer* e =
                texp->create_subexplorer( 0, "value_at_random_points" ) ;
      e->start_module_iterator() ;
      for( ; e->is_valid_module() ; e->go_next_module() )
      {
         PEL_ModuleExplorer* ee = e->create_subexplorer( 0 ) ;
         process_value_at_random_points( ref, ee ) ;
         ee->destroy() ; ee = 0 ;
      }
      e->destroy() ; e = 0 ;
   }
}

//---------------------------------------------------------------------------
void
PDE_ReferenceElement_TEST:: process_geometry( std::string const& name,
                                               PDE_ReferenceElement const* ref,
                                               PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElement_TEST:: process_geometry" ) ;
   PEL_CHECK( !name.empty() ) ;
   PEL_CHECK( ref!=0 ) ;
   PEL_CHECK( exp!=0 ) ;

   // Name :
   bool const n = ref->name()==name ;
   notify_one_test_result( "name", n ) ;

   // Dimension :
   bool const dim = ref->dimension()==(size_t) exp->int_data( "dimension" ) ;
   notify_one_test_result( "dimension", dim  ) ;

   // Reference polyhedron :
   bool const poly =
      ref->reference_polyhedron()->name()==
                                 exp->string_data( "reference_polyhedron" ) ;
   notify_one_test_result( "reference_polyhedron", poly ) ;
   
   // Number of nodes :
   bool const nb_nodes = ref->nb_nodes()==(size_t) exp->int_data( "nb_nodes" ) ;
   notify_one_test_result( "nb_nodes", nb_nodes ) ;

   // Data checker:
   if( exp->has_entry( "node_locations" ) )
   {
      if( ref->dimension()==0 )
      {
         PEL_Error::object()->raise_data_error(
            exp, "node_locations", "test is not expected in 0D" ) ;
      }
   }
   
   // Node locations :
   if( ref->dimension()!=0 )
   {
      doubleArray2D const& node_locs =
                  exp->doubleArray2D_data( "node_locations" ) ;
      if( node_locs.index_bound(0)!=ref->nb_nodes() ||
          node_locs.index_bound(1)!=ref->dimension() )
      {
         
         PEL_Error::object()->raise_data_error(
            exp, "node_locations", "bad size for the table" ) ;
      }
      bool ok_nodes = true ;
      GE_Point* pt = GE_Point::create( 0, ref->dimension() ) ;
      doubleVector coords(ref->dimension()) ;
      for( size_t i=0 ; i<ref->nb_nodes() ; ++i )
      {
         node_locs.extract_section( 0, i, coords ) ;
         pt->set_coordinates( coords ) ;
         ok_nodes &= ref->node_location( i )->is_equal( pt ) ;
         if( !ref->node_location( i )->is_equal( pt ) )
         {
            out() << "Expected node number " << i << " : " ;
            pt->print( out(), 0 ) ;
            out() << std::endl ;
            out() << "Node number " << i << " : " ;
            ref->node_location( i )->print( out(), 0 ) ;
            out() << std::endl ;
         }
      }
      pt->destroy() ; pt = 0 ;
      notify_one_test_result( "node_locations", ok_nodes ) ;
   }
}

//---------------------------------------------------------------------------
void
PDE_ReferenceElement_TEST:: process_bf_value_at_nodes(
                                         PDE_ReferenceElement const* ref,
                                         PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElement_TEST:: process_bf_value_at_nodes" ) ;
   PEL_CHECK( ref!=0 ) ;
   PEL_CHECK( exp!=0 ) ;

   double d_eps = exp->double_data( "dbl_epsilon" ) ;
   double d_min = exp->double_data( "dbl_minimum" ) ;

   bool okii = true ;
   bool okij = true ;
   for( size_t j=0 ; j<ref->nb_nodes() ; j++ )
   {
      GE_Point const* node = ref->node_location(j) ;
      for( size_t i=0 ; i<ref->nb_nodes() ; i++ )
      {
         double const Ni = ref->N_local(i,node) ;
         if( i==j )
         {
            bool const eq = PEL::double_equality( Ni, 1., d_eps, d_min ) ;
            okii &= eq ;
            if( !eq )
            {

               out() << "x_" << j << " = " ;
               node->print( out(), 0 ) ;
               out() << std::endl ;
               out() << "N_" << i << "(x_" << j << ") = "
                     << Ni << std::endl ;
            }
         }
         else
         {
            bool const eq = PEL::double_equality( Ni, 0., d_eps, d_min ) ;
            okij &= eq ;
            if( !eq )
            {
               out() << "x_" << j << " = " ;
               node->print( out(), 0 ) ;
               out() << std::endl ;
               out() << "N_" << i << "(x_" << j << ") = "
                     << Ni << std::endl ;
            }
         }
      }    
   }
   notify_one_test_result( "bf_at_nodes::Ni(xi)=1", okii ) ;
   notify_one_test_result( "bf_at_nodes::Ni(xj)=0", okij ) ;
}

//---------------------------------------------------------------------------
void
PDE_ReferenceElement_TEST:: process_bf_value_at_random_points(
                                         PDE_ReferenceElement const* ref,
                                         PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElement_TEST:: process_bf_value_at_random_points" ) ;
   PEL_CHECK( ref!=0 ) ;
   PEL_CHECK( exp!=0 ) ;

   PEL_ContextSimple* ct = PEL_ContextSimple::create( 0 ) ;
   PEL_DoubleVector* X = PEL_DoubleVector::create( ct, ref->dimension() ) ;
   ct->extend( PEL_Variable::object( "DV_X" ), X ) ;  
   
   size_t const nb_pts = (size_t) exp->int_data( "nb_random_points" ) ;
   PEL_DataWithContext* bfs = exp->abstract_data( 0, "basis_functions", ct ) ;
   if( bfs->data_type()!=PEL_Data::DoubleVector )
   {
      PEL_Error::object()->raise_bad_data_type(
            exp, "basis_functions", PEL_Data::DoubleVector  ) ;
   }
   if( !bfs->value_can_be_evaluated() )
   {
      PEL_Error::object()->raise_not_evaluable(
            exp, "basis_functions", bfs->undefined_variables() ) ;
   }
   double const d_eps = exp->double_data( "dbl_epsilon" ) ;
   double const d_min = exp->double_data( "dbl_minimum" ) ;
         
   GE_Point* pt = GE_Point::create( 0, ref->dimension() ) ;
   doubleVector coords( ref->dimension() ) ;
   bool ok_val = true ;
   initialize_randomizer() ;
   for( size_t ipt=0 ; ipt<nb_pts ; ++ipt )
   {
      // Search for a random point in the reference element :
      {
         set_random_point( ref, coords ) ;
         pt->set_coordinates( coords ) ;
         X->set( coords ) ;
      }

      // Value of the basis fonctions at this point :
      {
         doubleVector const& exp_val = bfs->to_double_vector() ;
         if( exp_val.size()!=ref->nb_nodes() )
         {
            PEL_Error::object()->raise_data_error(
               exp, "basis_functions", "bad size for the table" ) ;
         }
         else
         {
            for( size_t ibf=0 ; ibf<ref->nb_nodes() ; ++ibf )
            {
               double const val = ref->N_local( ibf, pt ) ;
               bool const eq = PEL::double_equality( val, exp_val(ibf),
                                                     d_eps, d_min ) ;
               ok_val &= eq ;
               if( !eq )
               {
                  out() << "Point : " ;
                  pt->print( out(), 0 ) ;
                  out() << std::endl ;
                  out() << "   N_exp(" << ibf << ",pt) = "
                        << exp_val(ibf) << std::endl ;
                  out() << "   N_get(" << ibf << ",pt) = "
                        << val << std::endl ;
                  if( PEL::abs( exp_val(ibf) )>d_min )
                  {
                     out() << "   epsilon : "
                           << PEL::abs( PEL::min( 1.-val/exp_val(ibf),
                                                  1.+val/exp_val(ibf) ) )
                           << std::endl ;
                  }
               }
            }
         }
      }
   }
   notify_one_test_result( "bf_at_random_points", ok_val ) ;
   bfs->destroy() ; bfs = 0 ;
   ct->destroy() ; ct = 0 ;
   pt->destroy() ; pt = 0 ;
}

//---------------------------------------------------------------------------
void
PDE_ReferenceElement_TEST:: process_dbf_value_at_random_points(
                                         PDE_ReferenceElement const* ref,
                                         PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElement_TEST:: process_dbf_value_at_random_points" ) ;
   PEL_CHECK( ref!=0 ) ;
   PEL_CHECK( exp!=0 ) ;
   
   size_t const nb_pts = (size_t) exp->int_data( "nb_random_points" ) ;
   double const d_eps = exp->double_data( "dbl_epsilon" ) ;
   double const d_min = exp->double_data( "dbl_minimum" ) ;
   double const epsilon = 1.E-8 ;
         
   GE_Point* pt = GE_Point::create( 0, ref->dimension() ) ;
   GE_Point* pt_eps = GE_Point::create( 0, ref->dimension() ) ;
   doubleVector coords( ref->dimension() ) ;
   bool ok_dval = true ;
   initialize_randomizer() ;
   for( size_t ipt=0 ; ipt<nb_pts ; ++ipt )
   {
      // Search for a random point in the reference element :
      {
         set_random_point( ref, coords ) ;
         pt->set_coordinates( coords ) ;
      }

      // Derivative value of the basis fonctions at this point :
      {
         for( size_t ibf=0 ; ibf<ref->nb_nodes() ; ++ibf )
         {
            double const val = ref->N_local( ibf, pt ) ;
            for( size_t d=0 ; d<ref->dimension() ; ++d )
            {
               pt_eps->set( pt ) ;
               double eps = epsilon ;
               pt_eps->set_coordinate( d, pt->coordinate(d)+eps ) ;
               if( !ref->reference_polyhedron()->contains( pt_eps ) )
               {
                  eps = -epsilon ;
                  pt_eps->set_coordinate( d, pt->coordinate(d)+eps ) ;
               }
               double const val_eps = ref->N_local( ibf, pt_eps ) ;
               double const dval = (val_eps-val)/eps ;
               double const dval_theo = ref->dN_local( ibf, d, pt ) ;
               bool const eq = PEL::double_equality( dval, dval_theo,
                                                     d_eps, d_min ) ;
               ok_dval &= eq ;
               if( !eq )
               {
                  out() << "Point : " ;
                  pt->print( out(), 0 ) ;
                  out() << std::endl ;
                  out() << "   dN_num(" << ibf << "," << d << ",pt) = "
                        << dval << std::endl ;
                  out() << "   dN_get(" << ibf << "," << d << ",pt) = "
                        << dval_theo << std::endl ;
                  if( PEL::abs( dval_theo )>d_min )
                  {
                     out() << "   epsilon : "
                           << PEL::abs( PEL::min( 1.-dval/dval_theo,
                                                  1.+dval/dval_theo ) )
                           << std::endl ;
                  }
               }
            }
         }
      }
   }
   notify_one_test_result( "dbf_at_random_points", ok_dval ) ;
   pt->destroy() ; pt = 0 ;
   pt_eps->destroy() ; pt_eps = 0 ;
}

//---------------------------------------------------------------------------
void
PDE_ReferenceElement_TEST:: process_d2bf_value_at_random_points(
                                         PDE_ReferenceElement const* ref,
                                         PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElement_TEST:: process_d2bf_value_at_random_points" ) ;
   PEL_CHECK( ref!=0 ) ;
   PEL_CHECK( exp!=0 ) ;
   
   size_t const nb_pts = (size_t) exp->int_data( "nb_random_points" ) ;
   double const d_eps = exp->double_data( "dbl_epsilon" ) ;
   double const d_min = exp->double_data( "dbl_minimum" ) ;
   double const epsilon = 1.E-8 ;
         
   GE_Point* pt = GE_Point::create( 0, ref->dimension() ) ;
   GE_Point* pt_eps = GE_Point::create( 0, ref->dimension() ) ;
   doubleVector coords( ref->dimension() ) ;
   bool ok_d2val = true ;
   initialize_randomizer() ;
   for( size_t ipt=0 ; ipt<nb_pts ; ++ipt )
   {
      // Search for a random point in the reference element :
      {
         set_random_point( ref, coords ) ;
         pt->set_coordinates( coords ) ;
      }

      // Derivative value of the basis fonctions at this point :
      {
         for( size_t ibf=0 ; ibf<ref->nb_nodes() ; ++ibf )
         {
            for( size_t c=0 ; c<ref->dimension() ; ++c )
            {
               double const val = ref->dN_local( ibf, c, pt ) ;
               for( size_t d=0 ; d<ref->dimension() ; ++d )
               {
                  pt_eps->set( pt ) ;
                  double eps = epsilon ;
                  pt_eps->set_coordinate( d, pt->coordinate(d)+eps ) ;
                  if( !ref->reference_polyhedron()->contains( pt_eps ) )
                  {
                     eps = -epsilon ;
                     pt_eps->set_coordinate( d, pt->coordinate(d)+eps ) ;
                  }
                  
                  double const val_eps = ref->dN_local( ibf, c, pt_eps ) ;
                  double const dval = (val_eps-val)/eps ;
                  double const dval_theo = ref->d2N_local( ibf, c, d, pt ) ;
                  bool const eq = PEL::double_equality( dval, dval_theo,
                                                        d_eps, d_min ) ;
                  ok_d2val &= eq ;
                  if( !eq )
                  {  
                     out() << "Point : " ;
                     pt->print( out(), 0 ) ;
                     out() << std::endl ;
                     out() << "   d2N_num(" << ibf << "," << c << "," << d
                           << ",pt) = " << dval << std::endl ;
                     out() << "   d2N_get(" << ibf << "," << c << "," << d
                           << ",pt) = "
                           << dval_theo << std::endl ;
                     if( PEL::abs( dval_theo )>d_min )
                     {
                        out() << "   epsilon : "
                              << PEL::abs( PEL::min( 1.-dval/dval_theo,
                                                     1.+dval/dval_theo ) )
                              << std::endl ;
                     }
                  }
               }
            }
         }
      }
   }
   notify_one_test_result( "d2bf_at_random_points", ok_d2val ) ;
   pt->destroy() ; pt = 0 ;
   pt_eps->destroy() ; pt_eps = 0 ;
}

//---------------------------------------------------------------------------
void
PDE_ReferenceElement_TEST:: process_value_at_random_points(
                                         PDE_ReferenceElement const* ref,
                                         PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElement_TEST:: process_value_at_random_points" ) ;
   PEL_CHECK( ref!=0 ) ;
   PEL_CHECK( exp!=0 ) ;

   size_t const dim = ref->dimension() ;
   size_t const nb_nodes = ref->nb_nodes() ;

   if( dim==0 )
   {
      PEL_Error::object()->raise_plain(
         "No \"value_at_random_points\" test with 0D element" ) ;
   }
   
   size_t const nb_pts = (size_t) exp->int_data( "nb_random_points" ) ;
   double const d_eps = exp->double_data( "dbl_epsilon" ) ;
   double const d_min = exp->double_data( "dbl_minimum" ) ;
   
   PEL_ContextSimple * ct = PEL_ContextSimple::create( 0 ) ;
   PEL_DoubleVector* X = PEL_DoubleVector::create( ct, dim  ) ;
   ct->extend( PEL_Variable::object( "DV_X" ), X ) ;
   PEL_DataWithContext* formula = exp->abstract_data( 0, "function", ct ) ;
   if( formula->data_type()!=PEL_Data::Double )
   {
      PEL_Error::object()->raise_bad_data_type(
         exp, "function", PEL_Data::Double  ) ;
   }
   if( !formula->value_can_be_evaluated(0) )
   {
      PEL_Error::object()->raise_not_evaluable(
         exp, "function", formula->undefined_variables(0) ) ;
   }
   PEL_DataWithContext* dformula = exp->abstract_data( 0, "gradient", ct ) ;
   if( dformula->data_type()!=PEL_Data::DoubleVector )
   {
      PEL_Error::object()->raise_bad_data_type(
         exp, "gradient", PEL_Data::DoubleVector  ) ;
   }
   if( !dformula->value_can_be_evaluated(0) )
   {
      PEL_Error::object()->raise_not_evaluable(
         exp, "gradient", dformula->undefined_variables(0) ) ;
   }
   PEL_DataWithContext* d2formula = exp->abstract_data( 0, "second_derivative", ct ) ;
   if( d2formula->data_type()!=PEL_Data::DoubleArray2D )
   {
      PEL_Error::object()->raise_bad_data_type(
         exp, "second_derivative", PEL_Data::DoubleArray2D ) ;
   }
   if( !d2formula->value_can_be_evaluated(0) )
   {
      PEL_Error::object()->raise_not_evaluable(
         exp, "second_derivative", d2formula->undefined_variables(0) ) ;
   }
   
   
   doubleVector coords( dim ) ;
   GE_Point* pt = GE_Point::create( 0, dim ) ;

   LA_Solver* solver = LA_GaussLU_DS::create( 0, 1.E-30 ) ;
   
   bool ok_val = true ;
   bool ok_dval = true ;
   bool ok_d2val = true ;

   // Nodes value :
   doubleVector Ni( nb_nodes ) ;
   {
      LA_DenseMatrix* mat = LA_DenseMatrix::create( 0, nb_nodes, nb_nodes ) ;
      LA_SeqVector* f = LA_SeqVector::create( 0, nb_nodes ) ;
      LA_SeqVector* x = LA_SeqVector::create( 0, nb_nodes ) ;
      for( size_t i=0 ; i<nb_nodes ; ++i )
      {
         GE_Point const* node = ref->node_location(i) ;
         for( size_t j=0 ; j<nb_nodes ; ++j )
         {
            mat->set_item( i, j, ref->N_local( j, node ) ) ;
         }
         X->set( node->coordinate_vector() ) ;
         f->set_item( i, formula->to_double() ) ;
      }
      mat->synchronize() ;
      f->synchronize() ;
      solver->set_matrix( mat ) ;
      solver->solve( f, x ) ;
      solver->unset_matrix() ;
      for( size_t i=0 ; i<nb_nodes ; ++i )
      {
         Ni( i ) = x->item( i ) ;
      }
      mat->destroy() ; mat = 0 ;
      f->destroy() ; f = 0 ;
      x->destroy() ; x = 0 ;
   }

   initialize_randomizer() ;
   
   for( size_t ipt=0 ; ipt<nb_pts ; ++ipt )
   {
      // Search for a random point in the reference element :
      {
         set_random_point( ref, coords ) ;
         pt->set_coordinates( coords ) ;
         X->set( coords ) ;
      }
      
      // Values at this point :
      {
         double val = 0. ;
         doubleVector dval( dim ) ;
         dval.set( 0. ) ;
         doubleArray2D d2val( dim, dim ) ;
         d2val.set( 0. ) ;
         for( size_t i=0 ; i<nb_nodes ; ++i )
         {
            val += Ni(i)*ref->N_local( i, pt ) ;
            for( size_t d=0 ; d<dim ; ++d )
            {
               dval(d) += Ni(i)*ref->dN_local( i, d, pt ) ;
               for( size_t d2=0 ; d2<dim ; ++d2 )
               {
                  d2val(d,d2) += Ni(i)*ref->d2N_local( i, d, d2, pt ) ;
               }
            }
         }

         // Check value :
         double const val_theo = formula->to_double() ;
         bool const eq = PEL::double_equality( val, val_theo,
                                               d_eps, d_min ) ;
         ok_val &= eq ;
         if( !eq )
         {  
            out() << "Point : " ;
            pt->print( out(), 0 ) ;
            out() << std::endl ;
            out() << "   val_computed = " << val << std::endl ;
            out() << "   val_formula  = " << val_theo << std::endl ;
            if( PEL::abs( val_theo )>d_min )
            {
               out() << "   epsilon : "
                     << PEL::abs( PEL::min( 1.-val/val_theo,
                                            1.+val/val_theo ) )
                     << std::endl ;
            }
         }

         // Check derivative :
         doubleVector const& dval_theo = dformula->to_double_vector() ;
         if( dval_theo.size()!=ref->dimension() )
         {
            PEL_Error::object()->raise_data_error(
               exp, "gradient",  "bad size for the table" ) ;
         }
         for( size_t c=0 ; c<dim ; ++c )
         {
            bool const eq2 = PEL::double_equality( dval(c), dval_theo(c),
                                                   d_eps, d_min ) ;
            ok_dval &= eq2 ;
            if( !eq2 )
            {  
               out() << "Point : " ;
               pt->print( out(), 0 ) ;
               out() << std::endl ;
               out() << "   grad_computed(" << c << ") = " << dval(c) << std::endl ;
               out() << "   grad_formula(" << c << ")  = " << dval_theo(c) << std::endl ;
               if( PEL::abs( dval_theo(c) )>d_min )
               {
                  out() << "   epsilon : "
                        << PEL::abs( PEL::min( 1.-dval(c)/dval_theo(c),
                                               1.+dval(c)/dval_theo(c) ) )
                        << std::endl ;
               }
            }
         }

         // Check second derivative :
         doubleArray2D const& d2val_theo = d2formula->to_double_array2D() ;
         if( d2val_theo.index_bound(0)!=dim ||
             d2val_theo.index_bound(1)!=dim )
         {
            PEL_Error::object()->raise_data_error(
               exp, "second_derivative",  "bad size for the table" ) ;
         }
         for( size_t c=0 ; c<dim ; ++c )
         {
            for( size_t d=0 ; d<dim ; ++d )
            {
               bool const eq2 = PEL::double_equality( d2val(c,d), d2val_theo(c,d),
                                                      d_eps, d_min ) ;
               ok_d2val &= eq2 ;
               if( !eq2 )
               {  
                  out() << "Point : " ;
                  pt->print( out(), 0 ) ;
                  out() << std::endl ;
                  out() << "   d2_computed(" << c << "," << d << ") = "
                        << d2val(c,d) << std::endl ;
                  out() << "   d2_formula(" << c << "," << d << ")  = "
                        << d2val_theo(c,d) << std::endl ;
                  if( PEL::abs( d2val_theo(c,d) )>d_min )
                  {
                     out() << "   epsilon : "
                           << PEL::abs( PEL::min( 1.-d2val(c,d)/d2val_theo(c,d),
                                                  1.+d2val(c,d)/d2val_theo(c,d) ) )
                           << std::endl ;
                  }
               }
            }
         }
      }
   }
   notify_one_test_result( exp->name()+"_value", ok_val ) ;
   notify_one_test_result( exp->name()+"_derivative", ok_dval ) ;
   notify_one_test_result( exp->name()+"_second_derivative", ok_d2val ) ;
   formula->destroy() ; formula = 0 ;
   dformula->destroy() ; dformula = 0 ;
   d2formula->destroy() ; d2formula = 0 ;
   ct->destroy() ; ct = 0 ;
   pt->destroy() ; pt = 0 ;
   solver->destroy() ; solver = 0 ;
}

//---------------------------------------------------------------------------
void
PDE_ReferenceElement_TEST:: initialize_randomizer( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElement_TEST:: initialize_randomizer" ) ;
   RANDOMIZER->start() ;
}

//---------------------------------------------------------------------------
void
PDE_ReferenceElement_TEST:: set_random_point( PDE_ReferenceElement const* ref,
                                              doubleVector& coords )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ReferenceElement_TEST:: set_random_point" ) ;
   PEL_CHECK( ref!=0 ) ;
   PEL_CHECK( coords.size()==ref->dimension() ) ;

   // Search for a random point in the reference element :
   coords.set( 0. ) ;
   double sum = 0. ;
   for( size_t v=0 ; v<ref->reference_polyhedron()->nb_vertices() ; ++v )
   {
      double const alpha = RANDOMIZER->item() ;
      sum += alpha ;
      RANDOMIZER->go_next() ;
      GE_Point const* vert = ref->reference_polyhedron()->vertex(v) ;
      for( size_t d=0 ; d<ref->dimension() ; ++d )
      {
         coords(d) += alpha*vert->coordinate(d) ;
      }
   }
   for( size_t d=0 ; d<ref->dimension() ; ++d )
   {
      coords(d) /= sum ;
   }
}
