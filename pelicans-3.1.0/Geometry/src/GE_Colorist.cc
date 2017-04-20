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

#include <GE_Colorist.hh>

#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_KeywordDataIterator.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ModuleIterator.hh>
#include <PEL_Variable.hh>
#include <PEL_Vector.hh>
#include <GE_Color.hh>
#include <GE_Meshing.hh>

//----------------------------------------------------------------------------
GE_Colorist*
GE_Colorist:: create( PEL_Object* a_owner, PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Colorist:: create" ) ;
   GE_Colorist* result = new GE_Colorist( a_owner, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}



//----------------------------------------------------------------------------
GE_Colorist:: GE_Colorist( PEL_Object* a_owner, PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
   : PEL_Object( a_owner ),
     DIM( 0 ),
     VERTICES( 0 ),
     V_COLORS( PEL_Vector::create( this, 0 ) ),
     V_EXPRES( PEL_Vector::create( this, 0 ) ),
     C_COLORS( PEL_Vector::create( this, 0 ) ),
     C_EXPRES( PEL_Vector::create( this, 0 ) ),
     S_COLORS( PEL_Vector::create( this, 0 ) ),
     S_EXPRES( PEL_Vector::create( this, 0 ) ),
     CONTEXT( 0 ),
     COORDS( 0 )
{
   CONTEXT = PEL_ContextSimple::create( this ) ;
   COORDS = PEL_DoubleVector::create( CONTEXT, doubleVector( 0 ) ) ;
   CONTEXT->extend( PEL_Variable::object("DV_X"), COORDS ) ;

   // check_module( exp ) ;

   read_colors( exp, "vertices", V_COLORS, V_EXPRES ) ;
   read_colors( exp, "cells", C_COLORS, C_EXPRES ) ;
   read_colors( exp, "faces", S_COLORS, S_EXPRES ) ;

}

//---------------------------------------------------------------------------
GE_Colorist:: ~GE_Colorist( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
GE_Colorist:: initialize( GE_Meshing* meshing )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Colorist:: initialize" ) ;

   if( VERTICES == 0 )
   {
      VERTICES = PEL_Vector::create( this, meshing->nb_vertices() ) ;
   }
   else
   {
      VERTICES->re_initialize( meshing->nb_vertices() ) ;
   }
   DIM = meshing->nb_space_dimensions() ;

   meshing->start_vertex_iterator() ;
   size_t i=0 ;
   for( ; meshing->valid_vertex() ; meshing->go_next_vertex() )
   {
      PEL_DoubleVector* pt = PEL_DoubleVector::create( VERTICES, 
                                       meshing->vertex_coordinates() ) ;
      PEL_CHECK( pt->to_double_vector().size() == DIM ) ;
      VERTICES->set_at( i, pt ) ;
      ++i ;
   }
   PEL_ASSERT( i == meshing->nb_vertices() ) ;
}

//---------------------------------------------------------------------------
bool
GE_Colorist:: is_initialized( void ) const
//---------------------------------------------------------------------------
{
   return( VERTICES != 0 ) ;
}

//---------------------------------------------------------------------------
GE_Color const*
GE_Colorist:: vertex_color( doubleVector const& coordinates ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Colorist:: vertex_color" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;

   GE_Color const* result = entity_color( coordinates, V_COLORS, V_EXPRES ) ;

   return( result ) ;
}

//---------------------------------------------------------------------------
GE_Color const*
GE_Colorist:: cell_color( size_t_vector const& vertex_indices ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Colorist:: cell_color" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;

   doubleVector coords( DIM ) ;
   compute_center( vertex_indices, coords ) ;

   GE_Color const* result = entity_color( coords, C_COLORS, C_EXPRES ) ;

   return( result ) ;
}

//---------------------------------------------------------------------------
GE_Color const*
GE_Colorist:: face_color( size_t_vector const& vertex_indices ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Colorist:: face_color" ) ;
   PEL_CHECK_PRE( is_initialized() ) ;

   doubleVector coords( DIM ) ;
   compute_center( vertex_indices, coords ) ;

   GE_Color const* result = entity_color( coords, S_COLORS, S_EXPRES )  ;

   return( result ) ;
}

//---------------------------------------------------------------------------
void
GE_Colorist:: read_colors( PEL_ModuleExplorer const* exp,
                           std::string const& module_name,
                           PEL_Vector* colors, PEL_Vector* expressions ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Colorist:: read_colors" ) ;
   if( exp->has_module( module_name ) )
   {
      PEL_ModuleExplorer* sexp = exp->create_subexplorer( 0, 
                                                          module_name ) ;
      sexp->start_entry_iterator() ;
      for( ; sexp->is_valid_entry() ; sexp->go_next_entry() )
      {
	 std::string const& str = sexp->keyword() ;
	 GE_Color::extend( str ) ;
         GE_Color const* col = GE_Color::object( str ) ;
         PEL_Data const* formula = sexp->data( colors, CONTEXT ) ;
         if( !formula->value_can_be_evaluated( CONTEXT ) )
         {
            PEL_Error::object()->raise_not_evaluable(
                 sexp, str, formula->undefined_variables( CONTEXT ) ) ;
         }
         colors->append( const_cast<GE_Color*>( col ) ) ;
         expressions->append( const_cast<PEL_Data*>( formula ) ) ;
      }
      sexp->destroy() ;
   }
}

//---------------------------------------------------------------------------
void
GE_Colorist:: compute_center( size_t_vector const& vertex_indices,
                              doubleVector& result ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Colorist:: compute_center" ) ;
   PEL_CHECK( result.size() == DIM ) ;
   PEL_CHECK( FORALL( ( size_t i=0 ; i<DIM ; ++i ), result(i)==0.0 ) ) ;

   size_t nbv = vertex_indices.size() ;
   for( size_t i=0 ; i<nbv ; ++i )
   {
      PEL_Object* oo = VERTICES->at( vertex_indices(i) ) ;
      PEL_DoubleVector* vv = static_cast<PEL_DoubleVector*>( oo ) ;
      for( size_t d=0 ; d<DIM ; ++d )
      {
         result(d) += vv->to_double_vector()(d) ;
      }
   }
   for( size_t d=0 ; d<DIM ; ++d )
   {
      result(d) /= (double)nbv ;
   }
}

//---------------------------------------------------------------------------
GE_Color const*
GE_Colorist:: entity_color( doubleVector const& coordinates,
                            PEL_Vector const* colors,
                            PEL_Vector const* expressions ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "GE_Colorist:: entity_color" ) ;
   
   COORDS->set( coordinates ) ;

   GE_Color const* result = 0 ;
   size_t i=0 ;
   for( ; i<expressions->index_limit() ; ++i )
   {
      PEL_Data const* ff = static_cast<PEL_Data*>( expressions->at(i) ) ;
      if( ff->to_bool( CONTEXT ) )
      {
         result = static_cast<GE_Color*>( colors->at( i ) ) ;
         break ;
      }
   }
   if( result != 0 )
   {
      for( ++i ; i<expressions->index_limit() ; ++i )
      {
         PEL_Data const* ff = static_cast<PEL_Data*>( expressions->at(i) ) ;
         if( ff->to_bool( CONTEXT ) )
         {
            //????? mettre plus joli
	    PEL_Error::object()->raise_plain( "GE_Colorist : conflict" ) ;
	 }
      }
   }
   return( result ) ;
}

// //---------------------------------------------------------------------------
// void
// GE_Colorist:: check_module( PEL_ModuleExplorer const* mod ) const
// //---------------------------------------------------------------------------
// {
//    PEL_ModuleIterator* itm = mod->create_module_iterator( 0 ) ;
//    for( itm->start() ; itm->is_valid() ; itm->go_next() )
//    {
//       std::string const& name = itm->item()->name() ; 
//       if( name!="vertices" && name!="cells" && name!="sides" )
//       {
// 	 std::string mesg = "Within module GE_Colorist :\n" ;
//          mesg +=            "invalid module name : " + name ;
// 	 PEL_Error::object()->raise_plain( mesg ) ;
//       }
//    }
//    itm->destroy() ;
//    PEL_KeywordDataIterator* ite = mod->create_entry_iterator( 0 ) ;
//    for( ite->start() ; ite->is_valid() ; ite->go_next() ) 
//    {
//       std::string const& name = ite->item()->keyword() ; 
//       std::string mesg = "Within module GE_Colorist :\n" ;
//       mesg +=            "invalid entry of keyword : " + name ;
//       PEL_Error::object()->raise_plain( mesg ) ;
//    }
//    ite->destroy() ;
// }


