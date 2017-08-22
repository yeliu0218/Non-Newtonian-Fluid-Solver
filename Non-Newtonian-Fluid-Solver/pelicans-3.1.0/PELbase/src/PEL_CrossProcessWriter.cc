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

#include <PEL_CrossProcessWriter.hh>

#include <size_t_array2D.hh>
#include <string>
#include <stringVector.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <doubleArray2D.hh>
#include <intArray2D.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleComparatorFloat.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Int.hh>
#include <PEL_IntVector.hh>
#include <PEL_IntArray2D.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_DoubleArray2D.hh>
#include <PEL_List.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_String.hh>
#include <PEL_StringVector.hh>
#include <PEL_assertions.hh>
#include <boolVector.hh>

#include <iostream>
#include <fstream>
#include <sstream>
using std::string ;


PEL_CrossProcessWriter const*
PEL_CrossProcessWriter::PROTOTYPE = new PEL_CrossProcessWriter() ;

struct PEL_CrossProcessWriter_ERROR
{
   static void n1( std::string const& fname,
                   size_t ic,
                   double val1, double val2 ) ;
   static void n3( std::string const& name,
                   double val1, double val2 ) ;
} ;

//----------------------------------------------------------------------
bool
PEL_CrossProcessWriter:: is_parallel_writer( void ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_CrossProcessWriter:: PEL_CrossProcessWriter( void )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( "PEL_CrossProcessWriter" )
   , SUB_WRITERS( 0 )
   , COM( 0 )
   , COORDS_CMP( 0 )
   , DBL_EPSI( PEL::bad_double() )
   , DBL_MINI( PEL::bad_double() )
   , loc2glob_VERTS( 0 )
   , loc2glob_FACES( 0 )
   , loc2glob_CELLS( 0 )
   , NB_VERTS( PEL::bad_index() )
   , NB_FACES( PEL::bad_index() )
   , NB_CELLS( PEL::bad_index() )
   , HALO_COLOR_INDEX( PEL::bad_index() )
{
   PEL_LABEL( "PEL_CrossProcessWriter:: PEL_CrossProcessWriter" ) ;
}

//----------------------------------------------------------------------
PEL_CrossProcessWriter*
PEL_CrossProcessWriter:: create_replica( PEL_Object* a_owner,
                                     PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CrossProcessWriter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PEL_CrossProcessWriter* result = new PEL_CrossProcessWriter( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}


//----------------------------------------------------------------------
PEL_CrossProcessWriter:: PEL_CrossProcessWriter( PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingWriter( a_owner )
   , SUB_WRITERS( PEL_List::create( this ) )
   , COM( PEL_Exec::communicator() )
   , COORDS_CMP( PEL_DoubleComparatorFloat::create( this, 1.e-12 ) )
   , DBL_EPSI( exp->has_entry( "dbl_epsilon" ) ?
                 exp->double_data( "dbl_epsilon" ) : 1.E-4 )
   , DBL_MINI( exp->has_entry( "dbl_minimum" ) ?
                 exp->double_data( "dbl_minimum" ) : 1.E-8 )
   , loc2glob_VERTS( 0 )
   , loc2glob_FACES( 0 )
   , loc2glob_CELLS( 0 )
   , NB_VERTS( PEL::bad_index() )
   , NB_FACES( PEL::bad_index() )
   , NB_CELLS( PEL::bad_index() )
   , HALO_COLOR_INDEX( -12 )
{
   if( DBL_EPSI <= 0. )
   {
      PEL_Error::object()->raise_bad_data_value(
         exp, "dbl_epsilon", "A positive value is expected" ) ;
   }
   if( DBL_MINI <= 0. )
   {
      PEL_Error::object()->raise_bad_data_value(
         exp, "dbl_minimum", "A positive value is expected" ) ;
   }
   
   if( COM->rank()==0 )
   {
      if( exp->has_entry( "sub_writers" ) )
      {
         stringVector const& subs = exp->stringVector_data( "sub_writers" ) ;
         for( size_t i=0 ; i<subs.size() ; ++i )
         {
            for( size_t j=i+1 ; j<subs.size() ; ++j )
            {
               if( subs(i) == subs(j) )
               {
                  PEL_Error::object()->raise_bad_data_value(
                     exp, "sub_writers",
                     "The writer \""+subs(i)+"\" is defined twice." ) ;
               }
            }
            PEL_DataOnMeshingWriter* sub =
               PEL_DataOnMeshingWriter::make( this, subs(i), exp ) ;
            if( sub->is_parallel_writer() )
            {
               PEL_Error::object()->raise_plain(
                  "Parallel writer \"" + subs(i) +"\"\n"
                  "can't be used as sub-writer of PEL_CrossProcessWriter instance" ) ;
            }
            SUB_WRITERS->append( sub ) ;
         }
      }
      else
      {
         // for the error message
         exp->stringVector_data( "sub_writers" ) ;
      }
   }
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_CrossProcessWriter:: ~PEL_CrossProcessWriter( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CrossProcessWriter:: ~PEL_CrossProcessWriter" ) ;
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//----------------------------------------------------------------------
void
PEL_CrossProcessWriter:: write_cycle( PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CrossProcessWriter:: write_cycle" ) ;
   PEL_CHECK_PRE( write_cycle_PRE( exp ) ) ;   

   PEL_Module* mod = collect_data( 0, exp ) ;

   if( COM->rank()==0 )
   {
      PEL_ModuleExplorer* explorer = PEL_ModuleExplorer::create( mod, mod ) ;
      for( size_t i=0 ; i<SUB_WRITERS->index_limit() ; i++ )
      {
         PEL_DataOnMeshingWriter * writer =
            static_cast<PEL_DataOnMeshingWriter*>( SUB_WRITERS->at( i ) ) ;
         writer->write_cycle( explorer ) ;
      }
      mod->destroy() ; mod=0 ;
   }
   else PEL_ASSERT( mod == 0 ) ;
}

//----------------------------------------------------------------------
PEL_Module*
PEL_CrossProcessWriter:: collect_data( PEL_Object* a_owner,
                                       PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CrossProcessWriter:: collect_data" ) ;
   
   size_t rank = COM->rank() ;

   PEL_Module* result = 0 ;
   if( rank == 0 )
   {
      result = PEL_Module::create( a_owner, exp->name() ) ;
      int cycle = exp->int_data( "cycle_number" ) ;
      result->add_entry( "cycle_number", PEL_Int::create( result, cycle ) ) ;
   }
                      
   // First, we save common scalar variables

   PEL_ModuleExplorer* sexp ;
   if( exp->has_module( "variables" ) )
   {
      sexp = exp->create_subexplorer( 0, "variables" ) ;
      size_t nb ;
      PEL_Module* variables = collect_variables( result, sexp, nb ) ;
      if( rank == 0 )
      {
         result->add_module( variables ) ;
         result->add_entry( "nb_variables",
                            PEL_Int::create( result, (int)nb ) ) ;
      }
      else PEL_ASSERT( variables == 0 ) ;
      sexp->destroy() ;
   }
   
   if( exp->has_module( "meshing" ) )
   {
      sexp = exp->create_subexplorer( 0, "meshing" ) ;
      PEL_Module* meshing = collect_meshing( result, sexp ) ;
      if( rank == 0 )
      {
         result->add_module( meshing ) ;
      }
      else PEL_ASSERT( meshing == 0 ) ;
      sexp->destroy() ;
   }
   
   if( exp->has_module( "fields" ) )
   {
      sexp = exp->create_subexplorer( 0, "fields" ) ;
      size_t nbf ;
      PEL_Module* fields = collect_field( result, sexp, nbf ) ;
      if( rank == 0 )
      {
         result->add_module( fields ) ;
         result->add_entry( "nb_fields",
                            PEL_Int::create( result, (int)nbf ) ) ;
      }
      else PEL_ASSERT( fields == 0 ) ;
      sexp->destroy() ;
   }
   if( exp->has_module( "integration_domain" ) )
   {
      sexp = exp->create_subexplorer( 0, "integration_domain" ) ;
      size_t nb ;
      PEL_Module* idom = collect_variables( result, sexp, nb ) ;
      if( rank == 0 )
      {
         result->add_module( idom ) ;
      }
      else PEL_ASSERT( idom == 0 ) ;
      sexp->destroy() ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Module*
PEL_CrossProcessWriter:: collect_variables( PEL_Object* a_owner,
                                            PEL_ModuleExplorer* exp,
                                            size_t& nb )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CrossProcessWriter:: collect_variables" ) ;
   
   PEL_Module* result = 0 ;
   size_t last = COM->nb_ranks()-1 ;
   size_t rank = COM->rank() ;
   nb = 0 ;

   if( rank == 0 )
   {
      result = PEL_Module::create( a_owner, exp->name() ) ;
   }
   
   exp->start_entry_iterator() ;
   for( ; exp->is_valid_entry() ; exp->go_next_entry() )
   {
      PEL_Data const* data = exp->data(0) ;
      string const& key = exp->keyword() ;
      bool same0 = false ;
      if( data->data_type()==PEL_Data::Double )
      {
         same0 = true ;
         double val = data->to_double() ;
         if( rank == 0  )
         {
            for( size_t i=1 ; i<=last ; i++ )
            {
               COM->send( i, val ) ;
               bool same ;
               COM->receive( i, same ) ;
               same0 = same0 && same ;
            }
         }
         else
         {
            double val0 ;
            COM->receive( 0, val0 ) ;
            bool same = PEL::double_equality( val0, val, DBL_EPSI, DBL_MINI ) ;
            if( !same )
            {
               PEL_CrossProcessWriter_ERROR:: n3( key, val0, val ) ;
            }
            COM->send( 0, same ) ;
         }
      }
      else if( data->data_type()==PEL_Data::DoubleVector )
      {
         same0 = true ;
         doubleVector const& val = data->to_double_vector() ;
         if( rank == 0  )
         {
            for( size_t i=1 ; i<=last ; i++ )
            {
               COM->send( i, val ) ;
               bool same ;
               COM->receive( i, same ) ;
               same0 = same0 && same ;
            }
         }
         else
         {
            doubleVector val0(0) ;
            COM->receive( 0, val0 ) ;
            bool same = val0.size()==val.size() ;
            for( size_t i=0 ; i<val.size() && same ; i++ )
            {
               same &= PEL::double_equality( val0(i), val(i),
                                             DBL_EPSI, DBL_MINI ) ;
               if( !same )
               {
                  PEL_CrossProcessWriter_ERROR:: n3( key, val0(i), val(i) ) ;
               }
            }
            COM->send( 0, same ) ;
         }
      }
      else if( data->data_type()==PEL_Data::DoubleArray2D )
      {
         same0 = true ;
         doubleArray2D const& val = data->to_double_array2D() ;
         if( rank == 0  )
         {
            for( size_t i=1 ; i<=last ; i++ )
            {
               COM->send( i, val ) ;
               bool same ;
               COM->receive( i, same ) ;
               same0 = same0 && same ;
            }
         }
         else
         {
            doubleArray2D val0(0,0) ;
            COM->receive( 0, val0 ) ;
            bool same = val0.index_bound(0)==val.index_bound(0) &&
                                     val0.index_bound(1)==val.index_bound(1) ;
            for( size_t i=0 ; i<val.index_bound(0) && same ; i++ )
            {
               for( size_t j=0 ; j<val.index_bound(1) && same ; j++ )
               {
                  same &= PEL::double_equality( val0(i,j), val(i,j),
                                                DBL_EPSI, DBL_MINI ) ;
                  if( !same )
                  {
                     PEL_CrossProcessWriter_ERROR:: n3( key, val0(i,j), val(i,j) ) ;
                  }
               }
            }
            COM->send( 0, same ) ;
         }
      }
      else if( data->data_type()==PEL_Data::Int )
      {
         same0 = true ;
         int val = data->to_int() ;
         if( rank == 0  )
         {
            for( size_t i=1 ; i<=last ; i++ )
            {
               COM->send( i, val ) ;
               bool same ;
               COM->receive( i, same ) ;
               same0 = same0 && same ;
            }
         }
         else
         {
            int val0 ;
            COM->receive( 0, val0 ) ;
            bool same = val0==val ;
            COM->send( 0, same ) ;
         }
      }
      else if( data->data_type()==PEL_Data::IntVector )
      {
         same0 = true ;
         intVector const& val = data->to_int_vector() ;
         if( rank == 0  )
         {
            for( size_t i=1 ; i<=last ; i++ )
            {
               COM->send( i, val ) ;
               bool same ;
               COM->receive( i, same ) ;
               same0 = same0 && same ;
            }
         }
         else
         {
            intVector val0(0) ;
            COM->receive( 0, val0 ) ;
            bool same = val0.size()==val.size() ;
            for( size_t i=0 ; i<val.size() && same ; i++ )
               same &= val0(i)==val(i) ;
            COM->send( 0, same ) ;
         }
      }
      else
      {
         std::ostringstream m ;
         m << "*** PEL_CrossProcessWriter : saving of \""
           << key << "\" values." << std::endl ;
         m << "    type \"" << data->data_type() << "\" not supported." << std::endl ;
         PEL_Error::object()->raise_plain( m.str() ) ;
      }

      if( same0 )
      {
         nb++ ;
         if( rank == 0 )
            result->add_entry( key, data->create_clone( result ) ) ;
      }
      else
      {
         std::ostringstream m ;
         m << "*** PEL_CrossProcessWriter : saving of \""
           << key << "\" values." << std::endl ;
         m << "    Not the same values on processes." << std::endl ;
         PEL_Error::object()->raise_plain( m.str() ) ;
      }
      data->destroy() ;
   }
   return result ;
}

//----------------------------------------------------------------------
PEL_Module*
PEL_CrossProcessWriter:: collect_meshing( PEL_Object* a_owner,
                                          PEL_ModuleExplorer* exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CrossProcessWriter:: collect_meshing" ) ;
   
   PEL_Module* result = 0 ;
   
   size_t rank = COM->rank() ;

   stringVector const& color_table = exp->stringVector_data( "color_table" ) ;
   intArray2D const& color_table_connectivity =
                           exp->intArray2D_data( "color_table_connectivity" ) ;
   std::string const& halo_color_name = exp->string_data( "halo_color_name" ) ;
   if( color_table.has( halo_color_name ) )
   {
      HALO_COLOR_INDEX = color_table.index_of( halo_color_name ) ;
   }
      
   //**Recovering cross-process vertices coordinates
   doubleArray2D verts = exp->doubleArray2D_data( "vertices" ) ;
   COM->merge( COORDS_CMP, verts, loc2glob_VERTS ) ;
   size_t glo_nb_verts = ( rank == 0 ) ? verts.index_bound(1) : 0 ;
   
   //**Recovering cross-process vertices color
   intVector glob_vert_color( glo_nb_verts ) ;
   compute_global_vector_of_indices( glob_vert_color, 
                                     exp->intVector_data( "vertex_color" ), 
                                     loc2glob_VERTS, HALO_COLOR_INDEX ) ;

   //**Recovering cross-process cell vertices
   intArray2D const& c2v = exp->intArray2D_data( "cell2vertex" ) ;
   size_t loc_nb_verts_per_cell = c2v.index_bound( 0 ) ;
   size_t loc_nb_cells = c2v.index_bound( 1 ) ;
   size_t glo_nb_verts_per_cell = COM->max( loc_nb_verts_per_cell ) ;   
   doubleArray2D d_cell2vert( glo_nb_verts_per_cell, loc_nb_cells ) ;
   for( size_t f = 0 ; f < loc_nb_cells ; f++ )
   {
      for( size_t i = 0 ; i < loc_nb_verts_per_cell ; i++ )
      {
         d_cell2vert( i, f ) = (double)( loc2glob_VERTS( c2v( i, f ) ) ) ; 
      }
   }
   COM->merge( COORDS_CMP, d_cell2vert, loc2glob_CELLS ) ;
   size_t glo_nb_cells = ( rank == 0 ) ? d_cell2vert.index_bound( 1 ) : 0 ;

   //**Recovering cross-process number of vertices per cell
   intVector const& cell_nb_verts = exp->intVector_data( "cell_nb_vertices" ) ;
   intVector glob_cell_nb_verts( glo_nb_cells ) ;
   compute_global_vector_of_indices( glob_cell_nb_verts, cell_nb_verts,
                                     loc2glob_CELLS, -100 ) ;

   //**Recovering cross-process cell colors
   intVector glob_cell_color( glo_nb_cells ) ;
   compute_global_vector_of_indices( glob_cell_color, 
                                     exp->intVector_data( "cell_color" ), 
                                     loc2glob_CELLS, HALO_COLOR_INDEX ) ;
   
   //**Recovering cross-process face vertices
   intArray2D const& f2v = exp->intArray2D_data( "face2vertex" ) ;
   size_t nb_verts_per_face = f2v.index_bound( 0 ) ;
   size_t loc_nb_faces = f2v.index_bound( 1 ) ;
   size_t glo_nb_verts_per_face = COM->max( nb_verts_per_face ) ;   
   doubleArray2D d_face2vert( glo_nb_verts_per_face, loc_nb_faces ) ;
   for( size_t f = 0 ; f < loc_nb_faces ; f++ )
   {
      for( size_t i = 0 ; i < nb_verts_per_face ; i++ )
      {
         d_face2vert( i, f ) = (double)( loc2glob_VERTS( f2v( i, f ) ) ) ;
      }
   }
   COM->merge( COORDS_CMP, d_face2vert, loc2glob_FACES ) ;    
   size_t glo_nb_faces = (rank == 0 ) ? d_face2vert.index_bound( 1 ) : 0 ;
   
   //**Recovering cross-process number of vertices per face
   intVector const& loc_face_nb_verts = exp->intVector_data( "face_nb_vertices" ) ;
   intVector glob_face_nb_verts( glo_nb_faces ) ;
   compute_global_vector_of_indices( glob_face_nb_verts, loc_face_nb_verts,
                                     loc2glob_FACES, -100 ) ;

   //**Recovering cross-process face colors
   intVector glob_face_color( glo_nb_faces ) ;
   compute_global_vector_of_indices( glob_face_color,
                                     exp->intVector_data( "face_color" ), 
                                     loc2glob_FACES, HALO_COLOR_INDEX ) ;

   //**Recovering cell<-->face correspondence
   //cell-->face
   intArray2D const& cell2face = exp->intArray2D_data( "cell2face" ) ;
   intVector glob_cell_nb_faces( glo_nb_cells ) ;
   size_t glo_nb_face_max_per_cell = COM->max( cell2face.index_bound( 0 ) ) ;
   if( rank > 0 ) glo_nb_face_max_per_cell = 0 ;
   
   intArray2D glob_cell2face( glo_nb_face_max_per_cell, glo_nb_cells ) ;
   compute_global_cell2face( glob_cell_nb_faces, glob_cell2face, 
                             exp->intVector_data( "cell_nb_faces" ),
                             cell2face, loc2glob_CELLS ) ;
   //face-->cell
   intArray2D const& face2cell = exp->intArray2D_data( "face2cell" ) ;
   size_t glo_nb_cell_max_per_face = COM->max( face2cell.index_bound( 0 ) ) ;
   if( rank > 0 ) glo_nb_cell_max_per_face = 0 ;
   intArray2D glob_face2cell( glo_nb_cell_max_per_face, glo_nb_faces ) ;
   compute_global_face2cell( glob_face2cell, face2cell, 
                             loc2glob_CELLS, loc2glob_FACES ) ;

   //Construction of meshing module on process 0
   if( rank == 0 )
   {
      NB_VERTS = glo_nb_verts ;         
      NB_FACES = glo_nb_faces ;
      NB_CELLS = glo_nb_cells ;

      intArray2D glob_cell2vert( glo_nb_verts_per_cell, glo_nb_cells ) ;
      
      for( size_t j=0 ; j < glo_nb_cells ; j++ )
         for( size_t s=0 ; s<glo_nb_verts_per_cell ; s++ )
            glob_cell2vert(s,j) = (int)d_cell2vert(s,j) ;

      intArray2D glob_face2vert( glo_nb_verts_per_face, glo_nb_faces ) ;

      for( size_t j=0 ; j < glo_nb_faces ; j++ )
         for( size_t s=0 ; s<glo_nb_verts_per_face ; s++ )
            glob_face2vert(s,j) = (int)d_face2vert(s,j) ;

      result = PEL_DataOnMeshingWriter::create_meshing_module( 
              a_owner, 
              exp->string_data( "name" ), exp->int_data( "nb_sp_dims" ),
              verts,
              glob_cell_nb_verts, glob_cell2vert,
              glob_cell_nb_faces, glob_cell2face,
              glob_face_nb_verts, glob_face2vert, glob_face2cell,
              color_table, halo_color_name,
              color_table_connectivity,
              glob_vert_color, glob_cell_color, glob_face_color ) ;

   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_CrossProcessWriter:: compute_global_vector_of_indices( 
                                        intVector& glob_vect_indices,
                                        intVector const& loc_vect_indices,
                                        size_t_vector const& loc2glob,
                                        int index_to_ignore ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CrossProcessWriter:: compute_global_vector_of_indices" ) ;

   size_t rank = COM->rank() ;
   size_t nb_ranks = COM->nb_ranks() ;
   
   if( rank==0 )
   {
      size_t_vector a_loc2glob( loc2glob ) ;
      intVector a_loc_vect_indices( loc_vect_indices ) ;
      
      for( size_t iproc=0 ; iproc<nb_ranks ; ++iproc )
      {
         if( iproc>0 )
         {
            COM->receive( iproc, a_loc2glob ) ;
            COM->receive( iproc, a_loc_vect_indices ) ;
         }
         
         for( size_t j=0 ; j<a_loc2glob.size() ; ++j )
         {
            if( a_loc_vect_indices( j ) != index_to_ignore )
               glob_vect_indices( a_loc2glob( j ) ) = a_loc_vect_indices( j ) ;
         }
      }
   }
   else
   {
      COM->send( 0, loc2glob ) ;
      COM->send( 0, loc_vect_indices ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_CrossProcessWriter:: compute_global_cell2face( 
                                   intVector& glob_cell_nb_faces,
                                   intArray2D& glob_cell2face,
                                   intVector const& loc_cell_nb_faces,
                                   intArray2D const& loc_cell2face,
                                   size_t_vector const& loc2glob_cells ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CrossProcessWriter:: compute_global_cell2face" ) ;

   size_t rank = COM->rank() ;
   size_t nb_ranks = COM->nb_ranks() ;
   
   if( rank==0 )
   {
      size_t_vector a_loc2glob_cells( loc2glob_cells ) ;
      intVector a_loc_cell_nb_faces( loc_cell_nb_faces ) ;
      intArray2D a_loc_cell2face( loc_cell2face ) ;
      
      for( size_t iproc=0 ; iproc<nb_ranks ; ++iproc )
      {
         if( iproc>0 )
         {
            COM->receive( iproc, a_loc2glob_cells ) ;
            COM->receive( iproc, a_loc_cell_nb_faces ) ;
            COM->receive( iproc, a_loc_cell2face ) ;
         }
         
         for( size_t j=0 ; j<a_loc2glob_cells.size() ; ++j )
         {
            size_t nn = a_loc_cell_nb_faces( j ) ;
            glob_cell_nb_faces( a_loc2glob_cells( j ) ) = nn ;
            for( size_t i=0 ; i<nn ; ++i )
            {
               glob_cell2face( i, a_loc2glob_cells( j ) ) = 
                                                    a_loc_cell2face( i, j ) ;
            }
         }
      }
   }
   else
   {
      COM->send( 0, loc2glob_cells ) ;
      COM->send( 0, loc_cell_nb_faces ) ;
      COM->send( 0, loc_cell2face ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_CrossProcessWriter:: compute_global_face2cell( 
                                   intArray2D& glob_face2cell,
                                   intArray2D const& loc_face2cell,
                                   size_t_vector const& loc2glob_cells,
                                   size_t_vector const& loc2glob_faces ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CrossProcessWriter:: compute_global_face2cell" ) ;

   size_t const rank = COM->rank() ;
   size_t const nb_ranks = COM->nb_ranks() ;

   size_t const nn = glob_face2cell.index_bound( 0 ) ;
   
   if( rank==0 )
   {
      size_t_vector a_loc2glob_cells( loc2glob_cells ) ;
      size_t_vector a_loc2glob_faces( loc2glob_faces ) ;
      intArray2D a_loc_face2cell( loc_face2cell ) ;
      
      for( size_t iproc=0 ; iproc<nb_ranks ; ++iproc )
      {
         if( iproc>0 )
         {
            COM->receive( iproc, a_loc2glob_cells ) ;
            COM->receive( iproc, a_loc2glob_faces ) ;
            COM->receive( iproc, a_loc_face2cell ) ;
         }
         
         for( size_t j=0 ; j<a_loc2glob_faces.size() ; ++j )
         {
            for( size_t i=0 ; i<nn ; ++i )
            {
               int i_loc = a_loc_face2cell( i, j ) ;

               int i_glob = index_for_trash() ;
               if( i_loc != index_for_trash() ) 
               {
                  i_glob = a_loc2glob_cells( i_loc ) ;
               }
               glob_face2cell( i, a_loc2glob_faces( j ) ) =  i_glob ;
            }
         }
      }
   }
   else
   {
      COM->send( 0, loc2glob_cells ) ;
      COM->send( 0, loc2glob_faces ) ;
      COM->send( 0, loc_face2cell ) ;
   }
}

//----------------------------------------------------------------------
PEL_Module*
PEL_CrossProcessWriter:: collect_field( PEL_Object* a_owner,
                                        PEL_ModuleExplorer* exp,
                                        size_t& nbf )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CrossProcessWriter:: collect_field" ) ;

   size_t const size = COM->nb_ranks() ;
   size_t const last = size-1 ;
   size_t const rank = COM->rank() ;
   
   PEL_Module* result = 0 ;
   if( rank == 0 )
   {
      result = PEL_Module::create( a_owner, "fields" ) ;
   }
   nbf = 0 ;
   
   exp->start_module_iterator() ;
   for( ; exp->is_valid_module() ; exp->go_next_module() )
   {
      PEL_ModuleExplorer* sexp = exp->create_subexplorer( 0 ) ;
      std::string const&  fname = sexp->string_data( "name" ) ;
      std::string const&  location = sexp->string_data( "location" ) ;
      doubleArray2D const& value = sexp->doubleArray2D_data( "value" ) ;
      size_t const nbcomp = value.index_bound(0) ;
      
      doubleArray2D val(0,0) ;
      
      if( rank==0 )
      {
         if( location=="at_vertices" ) 
            val.re_initialize( nbcomp, NB_VERTS ) ;
         else if( location == "at_cell_centers" )
            val.re_initialize( nbcomp, NB_CELLS ) ;
         else if( location == "at_face_centers" )
            val.re_initialize( nbcomp, NB_FACES ) ;
         val.set( undefined_value() ) ;
      }
      else
      {
         COM->receive( rank-1, val ) ;
      }

      if( location=="at_vertices" )
      {
         PEL_ASSERT( value.index_bound(1) == loc2glob_VERTS.size() ) ;
         size_t n = value.index_bound(1) ;
         for( size_t i=0 ; i<n ; i++ )
         {
            size_t idx = loc2glob_VERTS( i ) ;
            for( size_t j=0 ; j<nbcomp ; j++ )
            {
               if( value(j,i) != undefined_value() )
               {
                  if( val( j, idx ) != undefined_value() &&
                      !PEL::double_equality( val( j, idx ), value(j,i),
                                             DBL_EPSI, DBL_MINI ) )
                  {
                     PEL_CrossProcessWriter_ERROR::n1( fname, j,
                                                       val( j, idx ),
                                                       value(j,i) ) ;
                  }
                  val( j, idx ) = value(j,i) ;
               }
            }
         }
      }
      else if( location == "at_cell_centers" )
      {
         PEL_ASSERT( value.index_bound(1) == loc2glob_CELLS.size() ) ;
         size_t n = value.index_bound(1) ;
         for( size_t i=0 ; i<n ; i++ )
         {
            size_t idx = loc2glob_CELLS( i ) ;
            for( size_t j=0 ; j<nbcomp ; j++ )
            {
               if( value(j,i) != undefined_value() )
               {
                  if( val( j, idx ) != undefined_value() &&
                      !PEL::double_equality( val( j, idx ), value(j,i),
                                             DBL_EPSI, DBL_MINI ) )
                  {
                     PEL_CrossProcessWriter_ERROR::n1( fname, j,
                                                       val( j, idx ),
                                                       value(j,i) ) ;
                  }
                  val( j, loc2glob_CELLS( i ) ) = value(j,i) ;
               }
            }
         }
      }
      else
      {
         PEL_ASSERT( value.index_bound( 1 ) == loc2glob_FACES.size() ) ;
         size_t n = value.index_bound(1) ;
         for( size_t i=0 ; i<n ; i++ )
         {
            size_t idx = loc2glob_FACES( i ) ;
            for( size_t j=0 ; j<nbcomp ; j++ )
            {
               if( value(j,i) != undefined_value() )
               {
                  if( val( j, idx ) != undefined_value() &&
                      !PEL::double_equality( val( j, idx ), value(j,i),
                                             DBL_EPSI, DBL_MINI ) )
                  {
                     PEL_CrossProcessWriter_ERROR::n1( fname, j,
                                                       val( j, idx ),
                                                       value(j,i) ) ;
                  }
                  val( j, loc2glob_FACES( i ) ) = value(j, i) ;
               }
            }
         }
      }
      
      if( size>1 )
      {
         COM->send( (rank+1) % size, val ) ;
         if( rank==0 )
         {
            COM->receive( last, val ) ;
         }
      }

      if( rank==0 )
      {
         PEL_Module* fm = PEL_DataOnMeshingWriter::create_field_module( 
                                               result,
                                               fname,
                                               sexp->string_data( "meshing" ),
                                               location,
                                               val ) ;
         result->add_module( fm ) ;
      }
      
      sexp->destroy() ;
      nbf++ ;
   }
   
   return result ;
   
}

//internal--------------------------------------------------------------
void 
PEL_CrossProcessWriter_ERROR:: n1( std::string const& fname,
                                   size_t ic,
                                   double val1, double val2 )
//internal--------------------------------------------------------------
{
   std::ostringstream m ;
   m << "*** PEL_CrossProcessWriter : saving of \""
     << fname << "\"" << std::endl ;
   m << "      the value of the " << ic << "-th component" << std::endl ;
   m << "      is different when computed from two processes." << std::endl ;
   m << "         first  value : " << val1 << std::endl ;
   m << "         second value : " << val2 << std::endl ;
   m << "   If it is not an error, do adjust de value of" << std::endl ;
   m << "   the keywords \"dbl_minimum\" and \"dbl_epsilon\"." << std::endl ;
   PEL_Error::object()->raise_plain( m.str() ) ;
}

//internal--------------------------------------------------------------
void 
PEL_CrossProcessWriter_ERROR:: n3( std::string const& name,
                                   double val1, double val2 )
//internal--------------------------------------------------------------
{
   std::ostringstream m ;
   m << "*** PEL_CrossProcessWriter : saving of \""
     << name << "\"" << std::endl ;
   m << "      the value is different when computed from two processes." << std::endl ;
   m << "         first  value : " << val1 << std::endl ;
   m << "         second value : " << val2 << std::endl ;
   m << "   If it is not an error, do adjust de value of" << std::endl ;
   m << "   the keywords \"dbl_minimum\" and \"dbl_epsilon\"." << std::endl ;
   PEL_Error::object()->raise_plain( m.str() ) ;
}
