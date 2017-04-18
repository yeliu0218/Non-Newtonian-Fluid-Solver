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

#include <PEL_TICreader.hh>

#include <PEL.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleArray2D.hh>
#include <PEL_Error.hh>
#include <PEL_IntArray2D.hh>
#include <PEL_Int.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_String.hh>
#include <PEL_TICio.hh>
#include <PEL_assertions.hh>

#include <intArray2D.hh>
#include <intVector.hh>
#include <doubleArray2D.hh>
#include <doubleVector.hh>
#include <stringVector.hh>

#include <iostream>
#include <sstream>

PEL_TICreader const* PEL_TICreader::PROTOTYPE = new PEL_TICreader() ;

struct PEL_TICreader_ERROR
{
   static void n0( std::string const& mesg ) ;
} ;

//----------------------------------------------------------------------
PEL_TICreader:: PEL_TICreader( void )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingReader( "PEL_TICreader" )
   , MESHING_EXP( 0 )
   , FIELDS_EXP( 0 )
   , I_DOM_EXP( 0 )
   , VAR_EXP( 0 )
{
}

//----------------------------------------------------------------------
PEL_TICreader*
PEL_TICreader:: create_replica( PEL_Object* a_owner,
                                PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICreader:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PEL_TICreader* result = new PEL_TICreader( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_TICreader:: PEL_TICreader( PEL_Object* a_owner,
                               PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : PEL_DataOnMeshingReader( a_owner )
   , MESHING_EXP( 0 )
   , FIELDS_EXP( 0 )
   , I_DOM_EXP( 0 )
   , VAR_EXP( 0 )
{
   PEL_LABEL( "PEL_TICreader:: PEL_TICreader" ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t i_cycle = exp->int_data( "cycle" ) ;
   
   std::string TIC_file = exp->string_data( "files_basename" )+".gene" ;
   PEL_Module* mod1 =
             PEL_TICio::create_from_gene_file( 0, "MAIN", TIC_file, 1 ) ;
   PEL_Module* mod =
       PEL_TICio::create_from_gene_file( 0, "MAIN", TIC_file, i_cycle ) ;

   PEL_ModuleExplorer* m_exp = restore_cycle( mod, i_cycle ) ;
   PEL_ModuleExplorer* m_exp1 = restore_cycle( mod1, 1 ) ;
   
   // Meshing :
   {
      PEL_ModuleExplorer* meshing_exp = 0 ;
      if( m_exp->has_entry( "XNOD" ) )
      {
         meshing_exp = m_exp ;
      }
      else if( m_exp1->has_entry( "XNOD" ) )
      {
         meshing_exp = m_exp1 ;
      }
      else
      {
         PEL_TICreader_ERROR::n0( "no meshing found" ) ;
      }
      MESHING_EXP = create_meshing( meshing_exp ) ;
   }

   // Fields :
   {
      size_t const nb_vertices =
         MESHING_EXP->doubleArray2D_data( "vertices" ).index_bound(1) ;
      size_t const nb_meshes =
         MESHING_EXP->intArray2D_data( "cell2vertex" ).index_bound(1) ;
      FIELDS_EXP = create_fields( m_exp, nb_vertices, nb_meshes )  ;
   }

   // Integration domain :
   {
      PEL_ModuleExplorer* idom_exp = 0 ;
      if( m_exp->has_entry( "XFS" ) )
      {
         idom_exp = m_exp ;
      }
      else if( m_exp1->has_entry( "XFS" ) )
      {
         idom_exp = m_exp1 ;
      }
      if( idom_exp!=0 )
      {
         I_DOM_EXP = create_idom( idom_exp ) ;
      }
   }

   // Variables :
   // ...

   mod->destroy() ; mod = 0 ;
   mod1->destroy() ; mod1 = 0 ;
      
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_TICreader:: ~PEL_TICreader( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICreader:: ~PEL_TICreader" ) ;
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}

//----------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_TICreader:: meshing( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICreader:: meshing" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_ModuleExplorer* result = MESHING_EXP ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( meshing_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_TICreader:: fields( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICreader:: fields" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_ModuleExplorer* result = FIELDS_EXP ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( fields_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_TICreader:: integration_domain( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICreader:: integration_domain" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_ModuleExplorer* result = I_DOM_EXP ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( integration_domain_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_TICreader:: variables( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICreader:: variables" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   PEL_ModuleExplorer* result = VAR_EXP ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( variables_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_TICreader:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( PEL_DataOnMeshingReader::invariant() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_TICreader:: restore_cycle( PEL_Module* m, size_t i_cycle )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICreader:: restore_cycle" ) ;
   PEL_CHECK( m!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_ModuleExplorer* result = 0 ;
   std::ostringstream os ;
   os << "cycle_" << (int) i_cycle ;

   if( m->has_module( os.str() ) )
   {
      result = PEL_ModuleExplorer::create( m, m->module( os.str() ) ) ;
   }
   else
   {
      std::ostringstream msg ;
      msg << "Cycle " << i_cycle << " not found" ;
      PEL_TICreader_ERROR::n0( msg.str() ) ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result!=0 && result->owner()==m ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_TICreader:: create_meshing( PEL_ModuleExplorer* m_exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICreader:: create_meshing" ) ;
   PEL_CHECK( m_exp!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Module* m = PEL_Module::create( this , "meshing" ) ;
   PEL_ModuleExplorer* result = PEL_ModuleExplorer::create( this, m ) ;

   // Number of space dimensions :
   size_t nb_sp_dims = PEL::bad_index() ;
   if( m_exp->has_entry( "YNOD" ) )
   {
      nb_sp_dims = 3 ;
   }
   else if( m_exp->has_entry( "ZNOD" ) )
   {
      nb_sp_dims = 2 ;
   }
   else if( m_exp->has_entry( "XNOD" ) )
   {
      nb_sp_dims = 1 ;
   }
   else
   {
      PEL_TICreader_ERROR::n0( "\"XNOD\" not found" ) ;
   }
   m->add_entry( "nb_sp_dims", PEL_Int::create( m, (int) nb_sp_dims ) ) ;

   // Vertices :
   doubleVector const& x_nodes = m_exp->doubleVector_data( "XNOD" ) ;
   size_t const nb_verts = x_nodes.size() ;
   doubleArray2D vertices( nb_sp_dims, nb_verts ) ;
   vertices.set_section( 0, 0, x_nodes ) ;
   if( nb_sp_dims>=2 )
   {
      if( !m_exp->has_entry( "ZNOD" ) )
      {
         PEL_TICreader_ERROR::n0( "\"ZNOD\" not found" ) ;
      }
      doubleVector const& z_nodes = m_exp->doubleVector_data( "ZNOD" ) ;
      if( z_nodes.size()!=x_nodes.size() )
      {
         PEL_TICreader_ERROR::n0( "\"ZNOD\" has bad dimension" ) ;
      }
      vertices.set_section( 0, 1, z_nodes ) ;
   }
   if( nb_sp_dims>=3 )
   {
      if( !m_exp->has_entry( "YNOD" ) )
      {
         PEL_TICreader_ERROR::n0( "\"YNOD\" not found" ) ;
      }
      doubleVector const& y_nodes = m_exp->doubleVector_data( "YNOD" ) ;
      if( y_nodes.size()!=x_nodes.size() )
      {
         PEL_TICreader_ERROR::n0( "\"YNOD\" has bad dimension" ) ;
      }
      vertices.set_section( 0, 2, y_nodes ) ;
   }
   m->add_entry( "vertices", PEL_DoubleArray2D::create( m, vertices ) ) ;

   // Connectivity :
   if( !m_exp->has_entry( "NVER" ) )
   {
      PEL_TICreader_ERROR::n0( "\"NVER\" not found" ) ;
   }
   size_t const mesh_size = m_exp->int_data( "NVER" ) ;
   if( !m_exp->has_entry( "NODE" ) )
   {
      PEL_TICreader_ERROR::n0( "\"NODE\" not found" ) ;
   }
   intVector const& nodes = m_exp->intVector_data( "NODE" ) ;
   size_t nb_meshes = nodes.size()/mesh_size ;
   intArray2D connectivity( mesh_size, nb_meshes ) ;
   size_t i_vert=0 ;
   for( size_t i=0 ; i<nb_meshes ; i++ ) 
   {
      for (size_t j=0 ; j<mesh_size ; j++ ) 
      {
	 connectivity( j, i ) = nodes( i_vert++ )-1 ;
      }
   }
   m->add_entry( "cell2vertex", PEL_IntArray2D::create( m, connectivity ) ) ;

   // Mesh type : mesh reference element
   std::string mesh_type = "" ;
   if( nb_sp_dims==1 )
   {
      if( mesh_size==2 )
      {
         mesh_type = "GE_Segment" ;
      }
      else
      {
	 PEL_TICreader_ERROR::n0( "\"NVER\" bad value" ) ;
      }
   }
   else if( nb_sp_dims==2 )
   {
      if( mesh_size==3 )
      {
         mesh_type = "GE_ReferenceTriangle" ;
      }
      else if( mesh_size==4 )
      {
         mesh_type = "GE_ReferenceSquare" ;
      }
      else
      {
         PEL_TICreader_ERROR::n0( "\"NVER\" bad value" ) ;
      }
   }
   else if( nb_sp_dims==3 )
   {
      if( mesh_size==4 )
      {
         mesh_type = "GE_ReferenceTetrahedron" ;
      }
      else if( mesh_size==8 )
      {
         mesh_type = "GE_ReferenceCube" ;
      }
      else
      {
         PEL_TICreader_ERROR::n0( "\"NVER\" bad value" ) ;
      }
   }
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result!=0 && result->owner()==this ) ;
   PEL_CHECK_POST( result->name()=="meshing" ) ;
   PEL_CHECK_POST( result->has_entry( "nb_sp_dims" ) ) ;
   PEL_CHECK_POST( result->has_entry( "vertices" ) ) ;
   PEL_CHECK_POST( result->has_entry( "cell2vertex" ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_TICreader:: create_idom( PEL_ModuleExplorer* m_exp )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICreader:: create_idom" ) ;
   PEL_CHECK( m_exp!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Module* m = PEL_Module::create( this, "integration_domain" ) ;
   PEL_ModuleExplorer* result = PEL_ModuleExplorer::create( this, m ) ;

   // Inner boundary :
   if( !m_exp->has_entry( "XFS" ) )
   {
      PEL_TICreader_ERROR::n0( "\"XFS\" not found" ) ;
   }
   if( !m_exp->has_entry( "ZFS" ) )
   {
      PEL_TICreader_ERROR::n0( "\"ZFS\" not found" ) ;
   }
   doubleVector const& xFS = m_exp->doubleVector_data( "XFS" ) ;
   doubleVector const& zFS = m_exp->doubleVector_data( "ZFS" ) ;
   if( xFS.size()!=zFS.size() )
   {
      PEL_TICreader_ERROR::n0( "\"XFS\" and \"ZFS\" invalid dimensions" ) ;
   }
   doubleArray2D inner_boundary( (size_t) 2, xFS.size() ) ;
   inner_boundary.set_section( 0, 0, xFS ) ;
   inner_boundary.set_section( 0, 1, zFS ) ;
   m->add_entry( "inner_boundary",
                 PEL_DoubleArray2D::create( m, inner_boundary ) ) ;

   // Polygon :
   if( !m_exp->has_entry( "XDOM" ) )
   {
      PEL_TICreader_ERROR::n0( "\"XDOM\" not found" ) ;
   }
   if( !m_exp->has_entry( "ZDOM" ) )
   {
      PEL_TICreader_ERROR::n0( "\"ZDOM\" not found" ) ;
   }
   doubleVector const& xDOM = m_exp->doubleVector_data( "XDOM" ) ;
   doubleVector const& zDOM = m_exp->doubleVector_data( "ZDOM" ) ;
   if( xDOM.size()!=zDOM.size() )
   {
      PEL_TICreader_ERROR::n0( "\"XDOM\" and \"ZDOM\" invalid dimensions" ) ;
   }
   doubleArray2D polygon( (size_t) 2, xDOM.size() ) ;
   polygon.set_section( 0, 0, xDOM ) ;
   polygon.set_section( 0, 1, zDOM ) ;
   m->add_entry( "polygon", PEL_DoubleArray2D::create( m, polygon ) ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result!=0 && result->owner()==this ) ;
   PEL_CHECK_POST( result->name()=="integration_domain" ) ;
   PEL_CHECK_POST( result->has_entry( "inner_boundary" ) ) ;
   PEL_CHECK_POST( result->has_entry( "polygon" ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_ModuleExplorer*
PEL_TICreader:: create_fields( PEL_ModuleExplorer* m_exp,
                               size_t nb_vertices,
                               size_t nb_meshes )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICreader:: create_fields" ) ;
   PEL_CHECK( m_exp!=0  ) ;
   PEL_CHECK( nb_vertices>0 && nb_meshes>0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   PEL_Module* m = PEL_Module::create( this , "fields" ) ;
   PEL_ModuleExplorer* result = PEL_ModuleExplorer::create( this, m ) ;
   
   for( m_exp->start_entry_iterator() ;
        m_exp->is_valid_entry() ;
        m_exp->go_next_entry() )
   {
      PEL_Data const* data = m_exp->data( 0 ) ;
      if( data->data_type()==PEL_Data::DoubleVector )
      {
         doubleVector const& v = data->to_double_vector() ;
         std::string const& name = m_exp->keyword() ;

         // Location :
         std::string location = "not_a_field" ;
         size_t dim = PEL::bad_index() ;
         if( (v.size()%nb_vertices)==0 )
         {
            location = "at_vertices" ;
            dim = nb_vertices ;
         }
         else if( (v.size()%nb_meshes)==0 )
         {
            location = "at_cell_centers" ;
            dim = nb_meshes ;
         }

         if( location!="not_a_field" )
         {
            PEL_Module* f = PEL_Module::create( m, name ) ;
            m->add_module( f ) ;

            // Name :
            f->add_entry( "name", PEL_String::create( f, name ) ) ;
            f->add_entry( "type", PEL_String::create( f, "field" ) ) ;

            // Location and value :
            f->add_entry( "location", PEL_String::create( f, location ) ) ;

            // Value
            size_t const nb_comps = v.size()/dim ;
            doubleArray2D value( nb_comps, dim ) ;
            size_t i=0 ;
            for( size_t i_dim=0 ; i_dim<dim ; ++i_dim )
            {
               for( size_t i_comp=0 ; i_comp<nb_comps ; ++i_comp )
               {
                  value( i_comp, i_dim ) = v(i++) ;
               }
            }
            f->add_entry( "value", PEL_DoubleArray2D::create( f, value ) ) ;
         }
      }
      data->destroy() ; data = 0 ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result!=0 && result->owner()==this ) ;
   PEL_CHECK_POST( result->name()=="fields" ) ;
   return( result ) ;
}

//internal--------------------------------------------------------------
void 
PEL_TICreader_ERROR:: n0( std::string const& mesg )
//internal--------------------------------------------------------------
{
   std::ostringstream err ;
   err << "*** " << "PEL_TICreader:: invalid result file" << std::endl ;
   err << "*** " << mesg << std::endl ;
   PEL_Error::object()->raise_plain( err.str() ) ;
}
