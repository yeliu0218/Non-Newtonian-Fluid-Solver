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

#include <FE_FieldReader.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Communicator.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_Double.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>
#include <PEL_Vector.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <GE_Color.hh>
#include <GE_Point.hh>

#include <boolVector.hh>

#include <iostream>
#include <fstream>
#include <sstream>

struct FE_FieldReader_ERROR
{
   static void n0( std::string const& n ) ;
   static void n1( std::string const& n ) ;
   static void n2( std::string const& n ) ;
   static void n3( std::string const& f_name, size_t f_level,
                   std::string const& n ) ;
   static void n4( std::string const& f_name,
                   std::string const& n ) ;
   static void n5( std::string const& f_name,
                   std::string const& n ) ;
   static void n6( std::string const& f_name,
                   size_t i_node,
                   GE_Point const* pt1,
                   GE_Point const* pt2,
                   std::string const& n ) ;
} ;


FE_FieldReader const* FE_FieldReader::PROTOTYPE = new FE_FieldReader() ;

//----------------------------------------------------------------------
FE_FieldReader:: FE_FieldReader( void )
//----------------------------------------------------------------------
   : FE_OneStepIteration( "FE_FieldReader" )
   , DOM( 0 )
   , DBLE_EPS( PEL::bad_double() )
   , DBLE_MIN( PEL::bad_double() )
   , FIELDS( 0 )
   , EXPRS( 0 )
   , FIELD_LEVELS( 0 )
   , FIELD_RNAMES( 0 )
   , FIELD_RLEVELS( 0 )
   , CTX( 0 )
   , VAL( 0 )
{
}

//----------------------------------------------------------------------
FE_OneStepIteration*
FE_FieldReader:: create_replica( PEL_Object* a_owner,
                                 PDE_DomainAndFields const* dom,
                                 FE_SetOfParameters const* prms,
                                 PEL_ModuleExplorer* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldReader:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_FieldReader* result = new FE_FieldReader( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_FieldReader:: FE_FieldReader( PEL_Object* a_owner,
                                 PDE_DomainAndFields const* dom,
                                 FE_SetOfParameters const* prms,
                                 PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , DOM( dom )
   , DBLE_EPS( exp->has_entry( "dbl_epsilon" )
                  ? exp->double_data( "dbl_epsilon" ) : 1.E-6 )
   , DBLE_MIN( exp->has_entry( "dbl_minimum" )
                  ? exp->double_data( "dbl_minimum" ) : 1.E-8 )
   , OFILENAME( exp->string_data( "file_name" ) )
   , FIELDS( PEL_Vector::create( this, 0 ) )
   , EXPRS( PEL_Vector::create( this, 0 ) )
   , FIELD_LEVELS( 0 )
   , FIELD_RNAMES( 0 )
   , FIELD_RLEVELS( 0 )
   , CTX( 0 )
   , VAL( 0 )
{
   PEL_LABEL( "FE_FieldReader:: FE_FieldReader" ) ;

   // Node coordinates tolerance:
   {
      if( exp->has_entry( "dbl_epsilon" ) )
      {
         exp->test_data( "dbl_epsilon", "dbl_epsilon>0." ) ;
         exp->set_default( "dbl_epsilon", "1.e-6" ) ;
      }
      if( exp->has_entry( "dbl_minimum" ) )
      {
         exp->test_data( "dbl_minimum", "dbl_minimum>0." ) ;
         exp->set_default( "dbl_minimum", "1.e-8" ) ;
      }
   }
   
   PEL_Communicator const* com = PEL_Exec::communicator() ;
   if( com->nb_ranks()>1 )
   {
      std::ostringstream rank ;
      rank << "." << com->rank() ;
      OFILENAME += rank.str() ;
   }

   std::ifstream file( OFILENAME.c_str(), std::ios::in  ) ;
   if( !file ) FE_FieldReader_ERROR::n0( OFILENAME ) ;
   file.close() ;

   // Context:
   {
      PEL_ContextSimple* ctx_dummy =
                          PEL_ContextSimple::create( this ) ;
      VAL = PEL_Double::create( ctx_dummy, 0. ) ;
      ctx_dummy->extend( PEL_Variable::object( "DS_val" ), VAL ) ;
      CTX = ctx_dummy ;
   }

   // Fields:
   {
      PDE_SetOfDiscreteFields const* sdf = dom->set_of_discrete_fields() ;
      PEL_ModuleExplorer* sexp = exp->create_subexplorer( 0, "discrete_fields" ) ;
      for( sexp->start_module_iterator() ;
           sexp->is_valid_module() ;
           sexp->go_next_module() )
      {
         PEL_ModuleExplorer* ee = sexp->create_subexplorer( 0 ) ;
         std::string const& f_name = ee->string_data( "name" ) ;
         if( !sdf->has( f_name ) )
         {
            PEL_Error::object()->raise_data_error(
               ee, "name","   unknown discrete field \""+f_name+"\"" ) ;
         }
         PDE_DiscreteField* df = sdf->item( f_name ) ;
         FIELDS->append( df ) ;
         if( ee->has_entry( "level" ) )
         {
            int level = ee->int_data( "level" ) ;
            if( level<0 || level>=(int) df->storage_depth() )
            {
               PEL_Error::object()->raise_data_error(
                  ee, "level","   bad level for discrete field \""+f_name+"\"" ) ;
            }
            FIELD_LEVELS.append( (size_t) level ) ;
         }
         else
         {
            FIELD_LEVELS.append( PEL::bad_index() ) ;
         }
         int rlevel = 0 ;
         if( ee->has_entry( "stored_field_level" ) )
         {
            rlevel = ee->int_data( "stored_field_level" ) ;
            ee->set_default( "stored_field_level", "0" ) ;
            if( rlevel<0 )
            {
               PEL_Error::object()->raise_data_error(
                  ee, "stored_field_level","   bad level for discrete field \""+f_name+"\"" ) ;
            }
         }
         FIELD_RLEVELS.append( (size_t) rlevel ) ;
         std::string n = f_name ;
         if( ee->has_entry( "stored_field_name" ) )
         {
            n = ee->string_data( "stored_field_name" ) ;
         }
         FIELD_RNAMES.append( n ) ;
         if( ee->has_entry( "expression" ) )
         {
            PEL_Data* d = ee->abstract_data( EXPRS, "expression", CTX ) ;
            if( !d->value_can_be_evaluated( 0 ) )
            {
               PEL_Error::object()->raise_not_evaluable(
                  ee, "expression",
                  d->undefined_variables( 0 ) ) ;
            }
            if( d->data_type() != PEL_Data::Double )
            {
               PEL_Error::object()->raise_bad_data_value(
                  ee, "expression", "    a double value is expected" ) ;
            }
            EXPRS->append( d ) ;
         }
         else
         {
            EXPRS->resize( EXPRS->index_limit()+1 ) ;
         }
         ee->destroy() ; ee = 0 ;
      }
      sexp->destroy() ; sexp = 0 ;
   }   
}

//----------------------------------------------------------------------
FE_FieldReader:: ~FE_FieldReader( void ) 
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
FE_FieldReader:: do_before_time_stepping( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldReader:: do_before_time_stepping" ) ;
   PEL_CHECK_PRE( do_before_time_stepping_PRE( t_it ) ) ;

   start_total_timer( "FE_FieldReader:: do_before_time_stepping" ) ;
   
   PEL_Module* root = PEL_Module::create( 0, "ROOT", OFILENAME ) ;
   
   if( !root->has_module( "Root" ) ) FE_FieldReader_ERROR::n1( OFILENAME ) ;
   PEL_ModuleExplorer* exp =
               PEL_ModuleExplorer::create( root, root->module( "Root" )  ) ;
   
   if( !exp->has_entry( "nb_ranks" ) ) FE_FieldReader_ERROR::n1( OFILENAME ) ;
   size_t const nb_ranks = exp->int_data( "nb_ranks" ) ;
   if( !exp->has_entry( "rank" ) ) FE_FieldReader_ERROR::n1( OFILENAME ) ;
   size_t const rank = exp->int_data( "rank" ) ;
   PEL_Communicator const* com = PEL_Exec::communicator() ;
   if( nb_ranks != com->nb_ranks() || rank != com->rank() )
   {
      FE_FieldReader_ERROR::n2( OFILENAME ) ;
   }

   for( size_t i=0 ; i<FIELDS->index_limit() ; ++i )
   {
      restore_field( i, exp ) ;
   }
   
   root->destroy() ; root = 0 ;
   
   stop_total_timer() ;
}

//----------------------------------------------------------------------
void
FE_FieldReader:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldReader:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;
}

//----------------------------------------------------------------------
void
FE_FieldReader:: restore_field( size_t i_field, PEL_ModuleExplorer* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldReader:: restore_field" ) ;
   PEL_CHECK( i_field < FIELDS->index_limit() ) ;
   PEL_CHECK( exp != 0 ) ;
   
   PDE_DiscreteField* f =
              static_cast<PDE_DiscreteField*>( FIELDS->at( i_field ) ) ;
   size_t const nb_nodes = f->nb_nodes() ;
   size_t const nb_comps = f->nb_components() ;
   size_t const nb_sp_dims = DOM->nb_space_dimensions() ;
   size_t const level =
      ( FIELD_LEVELS( i_field ) == PEL::bad_index()
                                       ? 0 : FIELD_LEVELS( i_field ) ) ;
   
   std::string const& f_name = FIELD_RNAMES( i_field ) ;
   size_t const f_level = FIELD_RLEVELS( i_field ) ;
   
   PEL::out() << "         restoring field \"" << f->name() << "\"" << std::endl ;

   PEL_ModuleExplorer* ee = 0 ;

   for( exp->start_module_iterator() ;
        ee == 0 && exp->is_valid_module() ;
        exp->go_next_module() )
   {
      ee = exp->create_subexplorer( 0 ) ;
      if( !ee->has_entry( "name" ) )
                             FE_FieldReader_ERROR::n1( OFILENAME ) ;
      if( !ee->has_entry( "level" ) )
                             FE_FieldReader_ERROR::n1( OFILENAME ) ;
      if( ee->string_data( "name" ) != f_name ||
          ee->int_data( "level" ) != (int) f_level )
      {
         ee->destroy() ; ee = 0 ;
      }
   }

   if( ee == 0 ) FE_FieldReader_ERROR::n3( f_name, f_level, OFILENAME ) ;

   if( !ee->has_entry( "nb_components" ) )
                             FE_FieldReader_ERROR::n1( OFILENAME ) ;
   if( !ee->has_entry( "nb_nodes" ) )
                             FE_FieldReader_ERROR::n1( OFILENAME ) ;
   if( ee->int_data( "nb_components" ) != (int) nb_comps )
                     FE_FieldReader_ERROR::n4( f_name, OFILENAME ) ;
   if( ee->int_data( "nb_nodes" ) != (int) nb_nodes )
                     FE_FieldReader_ERROR::n4( f_name, OFILENAME ) ;
   if( !ee->has_entry( "node_values" ) )
                             FE_FieldReader_ERROR::n1( OFILENAME ) ;

   // Set field values:
   {
      PEL_Data* d = static_cast<PEL_Data*>( EXPRS->at( i_field ) ) ;
      doubleVector const& values =
                            ee->doubleVector_data( "node_values" ) ;
      if( values.size() != nb_comps*nb_nodes )
                             FE_FieldReader_ERROR::n1( OFILENAME ) ;
      size_t idx = 0 ;
      for( size_t i=0 ; i<nb_nodes ; ++i )
      {
         for( size_t ic=0 ; ic<nb_comps ; ++ic )
         {
            double val = values( idx++ ) ;
            if( d != 0 )
            {
               VAL->set( val ) ;
               val = d->to_double() ;
            }
            f->set_DOF_value( level, i, val, ic ) ;
         }
      }
      if( FIELD_LEVELS( i_field ) == PEL::bad_index() )
      {
         for( size_t l=1 ; l<f->storage_depth() ; ++l )
         {
            f->copy_DOFs_value( 0, l ) ;
         }
      }
   }

   // Check node locations:
   {
      if( !ee->has_entry( "node_locations" ) )
                             FE_FieldReader_ERROR::n1( OFILENAME ) ;
      doubleVector const& loc =
                         ee->doubleVector_data( "node_locations" ) ;
      if( loc.size() != nb_nodes*nb_sp_dims )
                             FE_FieldReader_ERROR::n1( OFILENAME ) ;
      if( !ee->has_entry( "has_node_locations" ) )
                             FE_FieldReader_ERROR::n1( OFILENAME ) ;
      boolVector const& has_loc =
                       ee->boolVector_data( "has_node_locations" ) ;
      if( has_loc.size() != nb_nodes )
                             FE_FieldReader_ERROR::n1( OFILENAME ) ;
      PDE_LocalFEcell* cFE = DOM->create_LocalFEcell( 0 ) ;
      cFE->include_color( GE_Color::halo_color() ) ;
      cFE->require_field_calculation( f, PDE_LocalFE::node ) ;
      for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
      {
         size_t const nb_loc_nodes = cFE->nb_local_nodes( f ) ;
         for( size_t i=0 ; i<nb_loc_nodes ; ++i )
         {
            size_t const i_node = cFE->global_node( f, i ) ;
            if( has_loc( i_node ) )
            {
               size_t const idx = i_node*nb_sp_dims ;
               GE_Point const* node = cFE->local_node_location( f, i ) ;
               for( size_t ic=0 ; ic<nb_sp_dims ; ++ic )
               {
                  bool ok = PEL::double_equality( node->coordinate(ic),
                                                  loc( idx+ic ),
                                                  DBLE_EPS, DBLE_MIN ) ;
                  if( !ok )
                  {
                     GE_Point* n_ini = GE_Point::create( this, nb_sp_dims ) ;
                     for( size_t iic=0 ; iic<nb_sp_dims ; ++iic )
                     {
                        n_ini->set_coordinate( iic, loc( idx+iic ) ) ;
                     }
                     FE_FieldReader_ERROR::n6(
                        f->name(), i_node, node, n_ini, OFILENAME ) ;
                     
                  }
               }
            }
         }
      }
      cFE->destroy() ; cFE = 0 ;
   }
   
   ee->destroy() ; ee = 0 ;
      
}

//internal--------------------------------------------------------------
void 
FE_FieldReader_ERROR:: n0( std::string const& n )
//internal--------------------------------------------------------------
{
   std::string mess ;
   mess += "*** FE_FieldReader error:\n" ;
   mess += "    unable to open file \"" ;
   mess += n ;
   mess += "\"" ;
   PEL_Error::object()->raise_plain( mess ) ;
}

//internal--------------------------------------------------------------
void 
FE_FieldReader_ERROR:: n1( std::string const& n )
//internal--------------------------------------------------------------
{
   std::string mess ;
   mess += "*** FE_FieldReader error:\n" ;
   mess += "    file \"" ;
   mess += n ;
   mess += "\" has an invalid structure" ;
   PEL_Error::object()->raise_plain( mess ) ;
}

//internal--------------------------------------------------------------
void 
FE_FieldReader_ERROR:: n2( std::string const& n )
//internal--------------------------------------------------------------
{
   std::string mess ;
   mess += "*** FE_FieldReader error:\n" ;
   mess += "    file \"" ;
   mess += n ;
   mess += "\" has incompatible parallel data\n" ;
   mess += "    (number of processes or index of current process)" ;
   PEL_Error::object()->raise_plain( mess ) ;
}

//internal--------------------------------------------------------------
void 
FE_FieldReader_ERROR:: n3( std::string const& f_name, size_t f_level,
                           std::string const& n )
//internal--------------------------------------------------------------
{
   std::ostringstream mess ;
   mess << "*** FE_FieldReader error:\n" ;
   mess << "    file \"" << n << "\"\n" ;
   mess << "    unable to find field saving for:\n" ;
   mess << "       field name: \"" << f_name << "\"\n" ;
   mess << "       field level: " << f_level ;
   PEL_Error::object()->raise_plain( mess.str() ) ;
}

//internal--------------------------------------------------------------
void 
FE_FieldReader_ERROR:: n4( std::string const& f_name,
                           std::string const& n )
//internal--------------------------------------------------------------
{
   std::ostringstream mess ;
   mess << "*** FE_FieldReader error:\n" ;
   mess << "    file \"" << n << "\"\n" ;
   mess << "    incompatible number of components:\n" ;
   mess << "       field name: \"" << f_name << "\"" ;
   PEL_Error::object()->raise_plain( mess.str() ) ;
}

//internal--------------------------------------------------------------
void 
FE_FieldReader_ERROR:: n5( std::string const& f_name,
                           std::string const& n )
//internal--------------------------------------------------------------
{
   std::ostringstream mess ;
   mess << "*** FE_FieldReader error:\n" ;
   mess << "    file \"" << n << "\"\n" ;
   mess << "    incompatible number of nodes:\n" ;
   mess << "       field name: \"" << f_name << "\"" ;
   PEL_Error::object()->raise_plain( mess.str() ) ;
}

//internal--------------------------------------------------------------
void 
FE_FieldReader_ERROR:: n6( std::string const& f_name,
                           size_t i_node,
                           GE_Point const* pt1,
                           GE_Point const* pt2,
                           std::string const& n )
//internal--------------------------------------------------------------
{
   std::ostringstream mess ;
   mess << "*** FE_FieldReader error:\n" ;
   mess << "    file \"" << n << "\"\n" ;
   mess << "    incompatible node location:\n" ;
   mess << "       field name: \"" << f_name << "\"" ;
   mess << "       node: " << i_node << std::endl ;
   mess << "       current location: " ;
   pt1->print( mess, 0 ) ;
   mess << std::endl ;
   mess << "       saving location: " ;
   pt2->print( mess, 0 ) ;
   PEL_Error::object()->raise_plain( mess.str() ) ;
}
