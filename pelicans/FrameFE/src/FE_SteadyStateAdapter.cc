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

#include <FE_SteadyStateAdapter.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <FE_TimeIterator.hh>

#include <fstream>
#include <iostream>
#include <iomanip>

FE_SteadyStateAdapter const*
FE_SteadyStateAdapter::PROTOTYPE = new FE_SteadyStateAdapter(
                                               "FE_SteadyStateAdapter" ) ;

//-------------------------------------------------------------------------
FE_SteadyStateAdapter:: FE_SteadyStateAdapter(
                                         std::string const& concrete_name )
//-------------------------------------------------------------------------
   : FE_TimeIteratorAdapter( concrete_name )
   , FIELDS( 0 )
   , LEVEL1( PEL::bad_index() )
   , LEVEL0( PEL::bad_index() )
   , FIELDS_TABLE( 0 )
   , FIELDS_ERROR( 0 )
   , OFILE_NAME( )
{
}

//-------------------------------------------------------------------------
FE_SteadyStateAdapter:: FE_SteadyStateAdapter(
                                     FE_TimeIterator* t_it,
                                     PDE_DomainAndFields const* dom,
                                     FE_SetOfParameters const* prms,
                                     PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : FE_TimeIteratorAdapter( t_it )
   , FIELDS( dom->set_of_discrete_fields() )
   , LEVEL1( (size_t) exp->int_data( "initial_level" ) )
   , LEVEL0( (size_t) exp->int_data( "current_level" ) )
   , FIELDS_TABLE( exp->stringVector_data( "discrete_fields" ) )
   , FIELDS_ERROR( exp->doubleVector_data( "minimal_error" ) )
   , OFILE_NAME( )
{
   for( size_t i=0 ; i<FIELDS_TABLE.size() ; ++i )
   {
      if( !FIELDS->has( FIELDS_TABLE(i) ) )
      {
         PEL_Error::object()->raise_bad_data_value(
            exp, "discrete_fields",
            "   unknown discrete field of name \""+FIELDS_TABLE(i)+"\"" ) ;
      }
      PDE_DiscreteField* f = FIELDS->item( FIELDS_TABLE(i) ) ;
      if( f->storage_depth() <= LEVEL1 ||
          f->storage_depth() <= LEVEL0 )
      {
         PEL_Error::object()->raise_plain(
            "*** FE_SteadyStateAdapter error:\n"
            "    bad storage depth for discrete field of name \""+FIELDS_TABLE(i)+"\"" ) ;
      }
   }
   if( FIELDS_ERROR.size() != FIELDS_TABLE.size() )
   {
      PEL_Error::object()->raise_bad_data_value(
         exp, "minimal_error",
         "   same size than \"discrete_fields\" is expected" ) ;
   }
   if( exp->has_module( "post_processing" ) )
   {
      PEL_ModuleExplorer const* eee =
                            exp->create_subexplorer( 0, "post_processing" ) ;
      OFILE_NAME = eee->string_data( "file_name" ) ;
      bool banner = true ;
      if( eee->has_entry( "banner" ) )
      {
         banner = eee->bool_data( "banner" ) ;
      }
      if( PEL_Exec::communicator()->rank() == 0 )
      {
         initialize_cv_file( banner ) ;
      }
            eee->destroy() ; eee = 0 ;
   }
}

//-------------------------------------------------------------------------
FE_SteadyStateAdapter:: ~FE_SteadyStateAdapter( void )
//-------------------------------------------------------------------------
{
   if( this == PROTOTYPE )
   {
      PROTOTYPE = 0 ;
   }
}

//-------------------------------------------------------------------------
void
FE_SteadyStateAdapter:: define_parameters_for_next_iteration( 
                           bool& finished, bool& restart, double& next_dt )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SteadyStateAdapter:: define_parameters_for_next_iteration" ) ;
   PEL_CHECK( define_parameters_for_next_iteration_PRE( finished, restart, next_dt ) ) ;

   PEL_Communicator const* com = PEL_Exec::communicator() ;
   doubleVector values( FIELDS_TABLE.size()+1 ) ;
   values(0) =  time_iterator()->time() ;
   finished = true ;
   for( size_t i=0 ; i<FIELDS_TABLE.size() ; ++i )
   {
      PDE_DiscreteField const* f = FIELDS->item( FIELDS_TABLE(i) ) ;
      double v_max = -PEL::bad_double() ;
      double err = -PEL::bad_double() ;
      for( size_t ic=0 ; ic<f->nb_components() ; ++ic )
      {
         for( size_t i_node=0 ; i_node<f->nb_nodes() ; ++i_node )
         {
            double const v_cur =
               f->DOF_value( LEVEL0, i_node, ic ) ;
            double const v_ini =
               f->DOF_value( LEVEL1, i_node, ic ) ;
            if( PEL::abs(v_ini) > v_max ) v_max = PEL::abs(v_ini) ;
            double dv = PEL::abs( v_cur-v_ini ) ;
            if( dv > err ) err = dv ;
         }
      }
      v_max = com->max( v_max ) ;
      err = com->max( err ) ;
      err /= PEL::max( 1., v_max ) ;
      if( err > FIELDS_ERROR(i) )
      {
         finished = false ;
      }
      values( i+1 ) = err ;
   }
   if( finished )
   {
      PEL_Error::object()->display_info(
         "*** FE_SteadyStateAdapter : steady state reached" ) ;
   }
   if( com->rank() == 0 && !OFILE_NAME.empty() )
   {
      save_in_cv_file( values ) ;
   }   
   
   PEL_CHECK( define_parameters_for_next_iteration_POST( finished, restart, next_dt ) ) ;
}

//-------------------------------------------------------------------------
FE_TimeIteratorAdapter*
FE_SteadyStateAdapter:: create_replica(
                                   FE_TimeIterator* t_it,
                                   PDE_DomainAndFields const* dom,
                                   FE_SetOfParameters const* prms,
                                   PEL_ModuleExplorer const* exp  ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SteadyStateAdapter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( t_it, dom, prms, exp ) ) ;
   
   FE_SteadyStateAdapter* result =
                      new FE_SteadyStateAdapter( t_it, dom, prms, exp ) ;
   
   PEL_CHECK( create_replica_POST( result, t_it, dom, prms, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
void
FE_SteadyStateAdapter:: initialize_cv_file( bool banner ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_SteadyStateAdapter:: initialize_cv_file" ) ;
   
   std::ofstream file( OFILE_NAME.c_str(),
                       std::ios::out | std::ios::trunc ) ;
   if( !file )
   {
      std::string mess ;
      mess += "*** FE_SteadyStateAdapter error:\n" ;
      mess += "    Saving failure : unable to open file \"" ;
      mess += OFILE_NAME ;
      mess += "\" for writing" ;
      PEL_Error::object()->raise_plain( mess ) ;
   }
   if( banner )
   {
      file << "#" << std::endl ;
      file << "# FE_SteadyStateAdapter generated file" << std::endl ;
      file << "#" << std::endl ;
      file << "#        time" ;
      for( size_t i=0 ; i<FIELDS_TABLE.size() ; ++i )
      {
         std::string s = FIELDS_TABLE(i) ;
         if( s.size()>10 )
         {
            s.erase( 10, s.size() ) ;
         }
         file << " # " << std::setw( 10 ) << s ;
      }
      file << std::endl ;
      file << "#" << std::endl ;
   }
   file.close() ;
}

//-------------------------------------------------------------------------
void
FE_SteadyStateAdapter::save_in_cv_file( doubleVector const& values ) const
//-------------------------------------------------------------------------
{
   std::ofstream os( OFILE_NAME.c_str(), std::ios::out | std::ios::app ) ;
   if( !os )
   {
      std::string mess ;
      mess += "*** FE_SteadyStateAdapter error:\n" ;
      mess += "    Saving failure : unable to open file \"" ;
      mess += OFILE_NAME ;
      mess += "\" for writing" ;      
   }
   os << std::setiosflags( std::ios::scientific )
      << std::setprecision( 4 ) ;
   for( size_t i=0 ; i<values.size() ; ++i )
   {
      os << "   " << values(i) ;
   }
   os << std::endl ;
   os.close() ;
}
