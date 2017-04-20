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

#include <FE_FieldSaver.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_BoolVector.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Int.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_String.hh>
#include <PEL_Vector.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <GE_Color.hh>
#include <GE_Point.hh>

#include <FE_TimeIterator.hh>

#include <iostream>
#include <fstream>
#include <sstream>

struct FE_FieldSaver_ERROR
{
   static void n0( std::string const& n ) ;
} ;

FE_FieldSaver const* FE_FieldSaver::PROTOTYPE = new FE_FieldSaver() ;

//----------------------------------------------------------------------
FE_FieldSaver:: FE_FieldSaver( void )
//----------------------------------------------------------------------
   : FE_OneStepIteration( "FE_FieldSaver" )
   , DOM( 0 )
   , DBLE_EPS( PEL::bad_double() )
   , DBLE_MIN( PEL::bad_double() )
   , FIELDS( 0 )
   , FIELD_LEVELS( 0 )
   , OTYPE( FE_FieldSaver::none )
   , OFILEFORMAT()
   , OFILEBNAME()
   , OFILENAME1()
   , OFILENAME2()
   , OFILENAME()
   , SAVING_NUMBER( PEL::bad_index() )
   , NEXT_SAVING_TIME( PEL::bad_double() )
   , SAVING_TIMES( 0 )
{
}

//----------------------------------------------------------------------
FE_OneStepIteration*
FE_FieldSaver:: create_replica( PEL_Object* a_owner,
                                PDE_DomainAndFields const* dom,
                                FE_SetOfParameters const* prms,
                                PEL_ModuleExplorer* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldSaver:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;
   
   FE_FieldSaver* result = new FE_FieldSaver( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_FieldSaver:: FE_FieldSaver( PEL_Object* a_owner,
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
   , FIELDS( PEL_Vector::create( this, 0 ) )
   , FIELD_LEVELS( 0 )
   , OTYPE( FE_FieldSaver::none )
   , OFILEFORMAT()
   , OFILEBNAME()
   , OFILENAME1()
   , OFILENAME2()
   , OFILENAME()
   , SAVING_NUMBER( PEL::bad_index() )
   , NEXT_SAVING_TIME( PEL::bad_double() )
   , SAVING_TIMES( 0 )
{
   PEL_LABEL( "FE_FieldSaver:: FE_FieldSaver" ) ;

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
     
   // Post-processing:
   {
      PEL_ModuleExplorer* sexp = exp->create_subexplorer( 0, "post_processing" ) ;
      
      std::string const& w_type = sexp->string_data( "type" ) ;
      if( w_type == "cycles_in_separate_files" )
      {
         OTYPE = FE_FieldSaver::per_one_cycle ;
         OFILEBNAME = sexp->string_data( "file_basename" ) ;
      }
      else if( w_type == "last_two_cycles" )
      {
         OTYPE = FE_FieldSaver::last_two_cycles ;
         OFILENAME1 = sexp->string_data( "file_name_0" ) ;
         OFILENAME2 = sexp->string_data( "file_name_1" ) ;
         PEL_Communicator const* com = PEL_Exec::communicator() ;
         if( com->nb_ranks()>1 )
         {
            std::ostringstream rank ;
            rank << "." << com->rank() ;
            OFILENAME1 += rank.str() ;
            OFILENAME2 += rank.str() ;
         }
         OFILENAME = OFILENAME1 ;
      }
      else
      {
         PEL_Error::object()->raise_bad_data_value(
            sexp, "type",
            "   - \"cycles_in_separate_files\"\n"
            "   - \"last_two_cycles\"" ) ;
      }

      OFILEFORMAT = "hybrid" ;
      if( sexp->has_entry( "output_format" ) )
      {
         OFILEFORMAT = sexp->string_data( "output_format" ) ;
         sexp->set_default( "output_format", "hybrid" ) ;
         sexp->test_data_in( "output_format", "hybrid,text" ) ;
      }
      
      if( sexp->has_entry( "saving_times" ) )
      {
         SAVING_TIMES = sexp->doubleVector_data( "saving_times" ) ;
      }
      
      sexp->destroy() ; sexp = 0 ;
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
         int level = 0 ;
         if( ee->has_entry( "level" ) )
         {
            level = ee->int_data( "level" ) ;
            ee->set_default( "level", "0" ) ;
            if( level<0 || level>=(int) df->storage_depth() )
            {
               PEL_Error::object()->raise_data_error(
                  ee, "level","   bad level for discrete field \""+f_name+"\"" ) ;
            }
         }
         FIELD_LEVELS.append( (size_t) level ) ;
         ee->destroy() ; ee = 0 ;
      }
      sexp->destroy() ; sexp = 0 ;
   }
}

//----------------------------------------------------------------------
FE_FieldSaver:: ~FE_FieldSaver( void ) 
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
FE_FieldSaver:: do_before_time_stepping( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldSaver:: do_before_time_stepping" ) ;
   PEL_CHECK_PRE( do_before_time_stepping_PRE( t_it ) ) ;
   
   // Initialize post_processing:
   if( SAVING_TIMES.size() != 0 )
   {
      if( !t_it->table_of_times_is_valid( SAVING_TIMES ) )
      {
         t_it->raise_invalid_table_of_times( "FE_FieldSaver", "saving_times", 
                                             SAVING_TIMES ) ;
      }
      SAVING_NUMBER = 0 ;
      if( SAVING_TIMES.has( t_it->time() ) )
      {
         save_field_value( ++SAVING_NUMBER, t_it ) ;
      }
      NEXT_SAVING_TIME = t_it->next_time_in_table( SAVING_TIMES ) ;
   }
}

//----------------------------------------------------------------------
void
FE_FieldSaver:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldSaver:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;
}

//----------------------------------------------------------------------
void
FE_FieldSaver:: do_after_time_adaptation( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldSaver:: do_after_time_adaptation" ) ;
   PEL_CHECK_PRE( do_after_time_adaptation_PRE( t_it ) ) ;
      
   if( SAVING_TIMES.size() != 0 )
   {
      if( FE_TimeIterator::greater_or_equal( t_it->time(),
                                             NEXT_SAVING_TIME ) )
      {
         start_total_timer( "FE_FieldSaver:: do_after_time_adaptation" ) ;
         
         save_field_value( ++SAVING_NUMBER, t_it ) ;
         NEXT_SAVING_TIME = t_it->next_time_in_table( SAVING_TIMES ) ;
         
         stop_total_timer() ;      
      }
   }
}

//----------------------------------------------------------------------
void
FE_FieldSaver:: save_other_than_time_and_fields(
                     FE_TimeIterator const* t_it,  PDE_ResultSaver* rs )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldSaver:: save_other_than_time_and_fields" ) ;
   PEL_CHECK_PRE( save_other_than_time_and_fields_PRE( t_it, rs ) ) ;

   if( SAVING_TIMES.size() == 0 )
   {
      start_total_timer( "FE_FieldSaver:: save_other_than_time_and_fields" ) ;
      
      save_field_value( rs->cycle_number(), t_it ) ;
      
      stop_total_timer() ;
   }      
}

//----------------------------------------------------------------------
void
FE_FieldSaver:: initialize_file( size_t i_cycle, double time )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldSaver:: initialize_file" ) ;

   PEL_Communicator const* com = PEL_Exec::communicator() ;
      
   if( OTYPE == per_one_cycle )
   {
      PEL_ASSERT( i_cycle<99999 ) ;
      std::ostringstream tmp ;
      tmp << i_cycle ;
      std::string nb_string = tmp.str() ;
      OFILENAME = OFILEBNAME+"_00000" ;
      OFILENAME.replace( OFILENAME.length()-nb_string.length(),
                         nb_string.length(), nb_string ) ;
      OFILENAME += ".pel" ;
      if( com->nb_ranks()>1 )
      {
         std::ostringstream rank ;
         rank << "." << com->rank() ;
         OFILENAME += rank.str() ;
      }
   }
   else if( OTYPE == last_two_cycles )
   {
      if( OFILENAME == OFILENAME1 )
      {
         OFILENAME = OFILENAME2 ;
      }
      else if( OFILENAME == OFILENAME2 )
      {
         OFILENAME = OFILENAME1 ;
      }
   }

   std::ofstream file( OFILENAME.c_str(), std::ios::out | std::ios::trunc ) ;
   if( !file ) FE_FieldSaver_ERROR:: n0( OFILENAME ) ;
   file << "//" << std::endl ;
   file << "//  FE_FieldSaver generated file" << std::endl ;
   file << "//" << std::endl ;
   file << "//      time: " << time << std::endl ;
   if( com->nb_ranks()>1 )
   {
      file << "//      rank: " << com->rank() << std::endl ;
   }
   file << "//" << std::endl ;
   file.close() ;
   if( OFILEFORMAT=="hybrid" )
   {
      std::string const bin_file_name =  OFILENAME+ ".bin" ;
      std::ofstream file_bin( bin_file_name.c_str(),
                              std::ios::out |
                              std::ios::binary | 
                              std::ios::trunc ) ;
      if( !file_bin ) FE_FieldSaver_ERROR:: n0( bin_file_name ) ;
      file_bin.close() ;
   }
}

//----------------------------------------------------------------------
void
FE_FieldSaver:: add_field( size_t i_field, PEL_Module* mod )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldSaver:: add_field" ) ;
   PEL_CHECK( i_field < FIELDS->index_limit() ) ;
   PEL_CHECK( mod != 0 ) ;

   PDE_DiscreteField const* f =
         static_cast<PDE_DiscreteField const*>( FIELDS->at( i_field ) ) ;
   size_t const level = FIELD_LEVELS( i_field ) ;
   size_t const nb_nodes = f->nb_nodes() ;
   size_t const nb_comps = f->nb_components() ;
   size_t const nb_sp_dims = DOM->nb_space_dimensions() ;

   std::ostringstream n ;
   n << "discrete_field#" << i_field ;
   PEL_Module* m = PEL_Module::create( mod, n.str() ) ;
   mod->add_module( m ) ;

   m->add_entry( "name", PEL_String::create( m, f->name() ) ) ;
   m->add_entry( "nb_components", PEL_Int::create( m, nb_comps ) ) ;
   m->add_entry( "nb_nodes", PEL_Int::create( m, nb_nodes ) ) ;
   m->add_entry( "level", PEL_Int::create( m, level ) ) ;
   
   doubleVector values( nb_comps*nb_nodes, PEL::bad_double() ) ;
   {
      size_t idx = 0 ;
      for( size_t i=0 ; i<nb_nodes ; ++i )
      {
         for( size_t ic=0 ; ic<nb_comps ; ++ic )
         {
            values( idx++ ) = f->DOF_value( level, i, ic ) ;
         }
      }
   }
   m->add_entry( "node_values", PEL_DoubleVector::create( m, values ) ) ;

   doubleVector loc( nb_sp_dims*nb_nodes, PEL::bad_double() ) ;
   boolVector found( nb_nodes, false ) ;
   boolVector has_loc( nb_nodes, true ) ;
   {
      PDE_LocalFEcell* cFE = DOM->create_LocalFEcell( 0 ) ;
      cFE->include_color( GE_Color::halo_color() ) ;
      cFE->require_field_calculation( f, PDE_LocalFE::node ) ;
      for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
      {
         size_t const nb_loc_nodes = cFE->nb_local_nodes( f ) ;
         for( size_t i=0 ; i<nb_loc_nodes ; ++i )
         {
            size_t const i_node = cFE->global_node( f, i ) ;
            size_t const idx = i_node*nb_sp_dims ;
            GE_Point const* node = cFE->local_node_location( f, i ) ;
            if( found(i_node) )
            {
               bool ok = true ;
               for( size_t ic=0 ; ok && ic<nb_sp_dims ; ++ic )
               {
                  ok = PEL::double_equality( node->coordinate(ic),
                                             loc( idx+ic ),
                                             DBLE_EPS, DBLE_MIN ) ;
               }
               if( !ok ) // Periodic node...
               {
                  has_loc(i_node) = false ;
                  for( size_t ic=0 ; ic<nb_sp_dims ; ++ic )
                  {
                     loc( idx+ic ) = 0. ;
                  }
               }
            }
            else
            {
               for( size_t ic=0 ; ic<nb_sp_dims ; ++ic )
               {
                  loc( idx+ic ) = node->coordinate(ic) ;
               }
            }      
            found(i_node) = true ;
         }
      }
      for( size_t i_node=0 ; i_node<nb_nodes ; ++i_node )
      {
         if( !found(i_node) )
         {
            has_loc(i_node) = false ;
            size_t const idx = i_node*nb_sp_dims ;
            for( size_t ic=0 ; ic<nb_sp_dims ; ++ic )
            {
               loc( idx+ic ) = 0. ;
            }
         }
      }
      cFE->destroy() ; cFE = 0 ;
   }
   m->add_entry( "node_locations", PEL_DoubleVector::create( m, loc ) ) ;
   m->add_entry( "has_node_locations", PEL_BoolVector::create( m, has_loc ) ) ;
}

//----------------------------------------------------------------------
void
FE_FieldSaver:: save_field_value(
                           size_t i_cycle, FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_FieldSaver:: save_field_value" ) ;
   
   initialize_file( i_cycle, t_it->time() ) ;

   PEL_Module* mod = PEL_Module::create( 0, "Root" ) ;
   
   mod->add_entry( "time", PEL_Double::create( mod, t_it->time() ) ) ;

   PEL_Communicator const* com = PEL_Exec::communicator() ;
   mod->add_entry( "nb_ranks", PEL_Int::create( mod, com->nb_ranks() ) ) ;
   mod->add_entry( "rank", PEL_Int::create( mod, com->rank() ) ) ;
   
   for( size_t i=0 ; i<FIELDS->index_limit() ; ++i )
   {
      add_field( i, mod ) ;
   }
   
   mod->write( OFILENAME, OFILEFORMAT ) ;

   mod->destroy() ; mod = 0 ;
}

//internal--------------------------------------------------------------
void 
FE_FieldSaver_ERROR:: n0( std::string const& n )
//internal--------------------------------------------------------------
{
   std::string mess ;
   mess += "*** FE_FieldSaver arror:\n" ;
   mess += "    unable to open file \"" ;
   mess += n ;
   mess += "\" for writing" ;
   PEL_Error::object()->raise_plain( mess ) ;
}

