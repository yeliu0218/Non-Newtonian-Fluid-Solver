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

#include <PDE_ValuesOnMeshing.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_SetOfPoints.hh>

#include <PDE_ReferenceElement.hh>
#include <PDE_ResultReader.hh>

#include <PEL_Communicator.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>
#include <PEL_Vector.hh>

#include <stringVector.hh>

#include <sstream>
using std::endl ;

struct PDE_ValuesOnMeshing_ERROR
{
   static void n0( std::string const& name,
                   std::string const& m_v,
                   std::string const& col ) ;
   static void n1( std::string const& name,
                   std::string const& col1,
                   std::string const& col2 ) ;
   static void n2_1( std::string const& name,
                     GE_Mpolyhedron const* mesh_poly,
                     GE_Color const* mesh_color ) ;
   static void n2_2( std::string const& name,
                     GE_Point const* pt,
                     GE_Color const* color ) ;
   static void n3( std::string const& name,
                   std::string const& col1, size_t nbc1,
                   std::string const& col2, size_t nbc2 ) ;
   static void n4( std::string const& name,
                   stringVector const& undef_vars  ) ;
   static void n5( std::string const& name, GE_Point const* pt ) ;
   static void n6( std::string const& name ) ;
   static void n7( std::string const& name ) ;
   static void n8( std::string const& name,
                   size_t current, size_t expected ) ;
   static void n9( std::string const& name,
                   size_t current, size_t expected ) ;
} ;

//----------------------------------------------------------------------
PDE_ValuesOnMeshing*
PDE_ValuesOnMeshing:: create(
                   PEL_Object* a_owner,
                   std::string const& a_field_name,
                   size_t nb_field_components,
                   PEL_ModuleExplorer const* exp, 
                   GE_SetOfPoints const* vs,
                   PDE_ResultReader const* reader )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ValuesOnMeshing:: create" ) ;
   PEL_CHECK_PRE( !a_field_name.empty() ) ;
   PEL_CHECK_PRE( exp!=0 ) ;
   PEL_CHECK_PRE( vs!=0 ) ;

   PDE_ValuesOnMeshing* result = new PDE_ValuesOnMeshing( a_owner,
                                                          a_field_name,
                                                          nb_field_components,
                                                          exp, vs, reader ) ;

   PEL_CHECK_POST( result!=0 && result->owner()==a_owner ) ;
   PEL_CHECK_POST( result->field_name() == a_field_name ) ;
   PEL_CHECK_POST( result-> nb_components() == nb_field_components ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_ValuesOnMeshing:: PDE_ValuesOnMeshing(
                   PEL_Object* a_owner,
                   std::string const& a_field_name,
                   size_t nb_field_components,
                   PEL_ModuleExplorer const* exp, 
                   GE_SetOfPoints const* vs,
                   PDE_ResultReader const* reader )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , NAME( a_field_name )
   , NB_COMPS( nb_field_components )
   , VM( vs )
   , READER( reader )
   , COLORS( PEL_Vector::create( this, 0 ) )
   , FORMULAS( PEL_Vector::create( this, 0 ) )
   , USED_COLORS( 0 )
   , DEFAULT_FORMULA( 0 )
   , CTX( 0 )
   , COORDS( 0 )
   , POLY( 0 )
   , CURRENT_FORMULA( 0 )
   , ELM( 0 )
   , VVALS( 0, 0 )
{
   PEL_LABEL( "PDE_ValuesOnMeshing:: PDE_ValuesOnMeshing" ) ;

   // Context for formula :
   {
      PEL_ContextSimple* c = PEL_ContextSimple::create( this ) ;
      COORDS = PEL_DoubleVector::create( c, doubleVector(0) ) ;
      c->extend( PEL_Variable::object( "DV_X" ), COORDS ) ;
      CTX = c ;
   }

   // Type of field definition :
   std::string tt = exp->string_data( "type" ) ;
   if( tt == "vertex_defined" ) TYPE = vertex_defined ;
   else if( tt == "mesh_defined" ) TYPE = mesh_defined ;
   else if( tt == "uniformly_defined" ) TYPE = uniformly_defined ;
   else if( tt == "defined_by_PDE_ResultReader" ) TYPE = reader_defined ;
   else
   {
      std::string const val = "   \"vertex_defined\"\n"
                              "   \"mesh_defined\"\n"
                              "   \"uniformly_defined\"\n"
                              "   \"defined_by_PDE_ResultReader\"" ;
      PEL_Error::object()->raise_bad_data_value( exp, "type", val ) ;
   }
 
   if( TYPE == vertex_defined || TYPE == mesh_defined )
   {
      PEL_ModuleExplorer* sexp = exp->create_subexplorer( 0, "value" ) ;
      sexp->start_entry_iterator() ;
      for( ; sexp->is_valid_entry() ; sexp->go_next_entry() )
      {
         std::string const& str = sexp->keyword() ;
         PEL_DataWithContext const* form = sexp->data( CTX, CTX ) ;
         if( form->data_type() != PEL_Data::DoubleVector )
         {
            PEL_Error::object()->raise_bad_data_type(
               sexp, str, PEL_Data::DoubleVector ) ;
         }
         if( str == "default" )
         {
            DEFAULT_FORMULA = form ;
         }
         else
         {
            GE_Color const* col = GE_Color::object( str ) ;
            COLORS->append( const_cast<GE_Color*>( col ) ) ;
            FORMULAS->append( const_cast<PEL_DataWithContext*>( form ) ) ;
         }
      }
      sexp->destroy() ;
      USED_COLORS.re_initialize( COLORS->index_limit() ) ;
   }
   else if( TYPE == reader_defined )
   {
      if( READER == 0 )
      {
         PDE_ValuesOnMeshing_ERROR:: n6( NAME ) ;
      }
      if( !READER->has_field( NAME ) )
      {
         PDE_ValuesOnMeshing_ERROR:: n7( NAME ) ;
      }
      
      // Value outside the stored grid
      if( exp->has_entry( "out_of_stored_grid_value" ) )
      {
         DEFAULT_FORMULA = exp->abstract_data(
            CTX, "out_of_stored_grid_value", CTX ) ;
         if( DEFAULT_FORMULA->data_type() != PEL_Data::DoubleVector )
         {
            PEL_Error::object()->raise_bad_data_type(
               exp, "out_of_stored_grid_value", PEL_Data::DoubleVector ) ;
         }
      }
   }
   else if( TYPE == uniformly_defined )
   {
      DEFAULT_FORMULA = exp->abstract_data( CTX, "value", CTX ) ;
      if( DEFAULT_FORMULA->data_type() != PEL_Data::DoubleVector )
      {
         PEL_Error::object()->raise_bad_data_type( exp, "value", 
                                                   PEL_Data::DoubleVector ) ;
      }
   }
}

//----------------------------------------------------------------------
PDE_ValuesOnMeshing:: ~PDE_ValuesOnMeshing( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ValuesOnMeshing:: ~PDE_ValuesOnMeshing" ) ;

   for( size_t i=0 ; i<USED_COLORS.size() ; ++i )
   {
      bool used = USED_COLORS(i) ;
      used = PEL_Exec::communicator()->boolean_or( used ) ;
      if( !used )
      {
         std::string m_v = ( TYPE==vertex_defined ) ? "vertex" : "mesh" ;
         GE_Color const* col = static_cast<GE_Color*>( COLORS->at(i) ) ;
         PDE_ValuesOnMeshing_ERROR:: n0( NAME, m_v, col->name() ) ;
      }
   }
}

//----------------------------------------------------------------------
std::string const&
PDE_ValuesOnMeshing:: field_name( void ) const
//----------------------------------------------------------------------
{
   return( NAME ) ;
}

//----------------------------------------------------------------------
size_t
PDE_ValuesOnMeshing:: nb_components( void ) const
//----------------------------------------------------------------------
{
   return( NB_COMPS ) ;
}

//----------------------------------------------------------------------
PEL_DataWithContext const*
PDE_ValuesOnMeshing:: formula( GE_Color const* color )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ValuesOnMeshing:: formula" ) ;
   PEL_CHECK( color!=0 ) ;
   
   PEL_DataWithContext const* result = 0 ;
   size_t i=0 ;
   for( ; i<COLORS->index_limit() ; ++i )
   {
      GE_Color const* col = static_cast<GE_Color*>( COLORS->at(i) ) ;
      if( col->is_matching( color ) )
      {
         USED_COLORS(i) = true ;
         result = static_cast<PEL_DataWithContext*>( FORMULAS->at(i) ) ;
         break ;
      }
   }
   for( ++i ; i<COLORS->index_limit() ; ++i )
   {
      GE_Color const* col = static_cast<GE_Color*>( COLORS->at(i) ) ;
      if( col->is_matching( color ) )
      {
         PDE_ValuesOnMeshing_ERROR:: n1( NAME, col->name(), color->name() ) ;
      }
   }
   if( result == 0 ) result = DEFAULT_FORMULA ;

   PEL_CHECK_POST( IMPLIES( TYPE!=reader_defined, result!=0 ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_ValuesOnMeshing:: set_mesh( GE_Mpolyhedron const* mesh_poly,
                                GE_Color const* mesh_color )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ValuesOnMeshing:: set_mesh" ) ;
   PEL_CHECK_PRE( mesh_poly!=0 ) ;
   PEL_CHECK_PRE( mesh_color!=0 ) ;
   
   POLY = mesh_poly ;
   CURRENT_FORMULA = 0 ;
   ELM = 0 ;

   if( TYPE==mesh_defined ||
       TYPE==uniformly_defined )
   {
      CURRENT_FORMULA = formula( mesh_color ) ;
      if( CURRENT_FORMULA == 0 ) 
         PDE_ValuesOnMeshing_ERROR:: n2_1( NAME, mesh_poly, mesh_color ) ;
   }
   else if( TYPE==reader_defined )
   {
      CURRENT_FORMULA = formula( mesh_color ) ;
   }
   else if( TYPE==vertex_defined )
   {
      ELM = PDE_ReferenceElement::object_with_nodes_at_vertices( POLY ) ;
      PEL_CHECK( ELM->nb_nodes() == POLY->nb_vertices() ) ;
      GE_Point* node_pt = //?????? a chaque fois ???????
                         GE_Point::create( 0, POLY->nb_space_dimensions()  ) ;
      GE_Color const* col_0 = 0 ;
      for( size_t in=0 ; in<ELM->nb_nodes() ; ++in )
      {
         GE_Point const* pt_ref = ELM->node_location( in ) ;
         POLY->apply_mapping( pt_ref, node_pt ) ;
         PEL_ASSERT( VM->has( node_pt ) ) ;
         GE_Color const* col = VM->color( VM->index( node_pt ) ) ;
         if( in == 0 )
         {
            col_0 = col ;
         }
         COORDS->set( node_pt->coordinate_vector() ) ;
         PEL_DataWithContext const* fo = formula( col ) ;
         if( fo == 0 )
         {
            PDE_ValuesOnMeshing_ERROR:: n2_2( NAME, node_pt, col ) ;
         }
         doubleVector const& val = fo->to_double_vector() ;
         if( in == 0 )
         {
            VVALS.re_initialize( ELM->nb_nodes(), val.size() ) ;
         }
         else
         {
            PEL_CHECK( VVALS.index_bound(0) == ELM->nb_nodes() ) ;
            if( VVALS.index_bound(1) != val.size() )
            {
               PDE_ValuesOnMeshing_ERROR:: n3( NAME,
                                               col_0->name(),
                                               VVALS.index_bound(1),
                                               col->name(),
                                               val.size() ) ;
            }
         }
         VVALS.set_section( 0, in, val ) ;
      }
      node_pt->destroy() ; node_pt = 0 ;
   }
}

//----------------------------------------------------------------------
void
PDE_ValuesOnMeshing:: compute_value( GE_Point const* pt,
                                     doubleVector& result ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ValuesOnMeshing:: value" ) ;
   PEL_CHECK( IMPLIES( TYPE != reader_defined,
                       (ELM==0 && CURRENT_FORMULA!=0) || 
                       (ELM!=0 && CURRENT_FORMULA==0) ) ) ;
   PEL_CHECK( POLY->contains( pt ) ) ;

   if( TYPE == reader_defined )
   {
      if( READER->is_in_grid( pt ) )
      {
         result = READER->field_value( NAME, pt ) ;
         if( result.size() != NB_COMPS )
         {
            PDE_ValuesOnMeshing_ERROR:: n8( NAME, result.size(), NB_COMPS ) ;
         }
      }
      else if( CURRENT_FORMULA != 0 )
      {
         COORDS->set( pt->coordinate_vector() ) ;
         if( !CURRENT_FORMULA->value_can_be_evaluated() )
         {
            PDE_ValuesOnMeshing_ERROR:: n4(
               NAME, CURRENT_FORMULA->undefined_variables() ) ;
         }
         result = CURRENT_FORMULA->to_double_vector() ;
         if( result.size() != NB_COMPS )
         {
            PDE_ValuesOnMeshing_ERROR:: n9( NAME, result.size(), NB_COMPS ) ;
         }
      }
      else
      {
         PDE_ValuesOnMeshing_ERROR:: n5( NAME, pt ) ;
      }
   }
   else if( CURRENT_FORMULA != 0 ) 
   {
      COORDS->set( pt->coordinate_vector() ) ;
      if( !CURRENT_FORMULA->value_can_be_evaluated() )
      {
         PDE_ValuesOnMeshing_ERROR:: n4(
            NAME, CURRENT_FORMULA->undefined_variables() ) ;
      }
      result = CURRENT_FORMULA->to_double_vector()  ;
      if( result.size() != NB_COMPS )
      {
         PDE_ValuesOnMeshing_ERROR:: n9( NAME, result.size(), NB_COMPS ) ;
      }
   }
   else
   {
      PEL_CHECK( TYPE == vertex_defined ) ;
      
      result.re_initialize( VVALS.index_bound(1) ) ;
      if( result.size() != NB_COMPS )
      {
         PDE_ValuesOnMeshing_ERROR:: n9( NAME, result.size(), NB_COMPS ) ;
      }

      GE_Point* pt_ref = GE_Point::create( 0, POLY->dimension() ) ;
      POLY->apply_inverse_mapping( pt, pt_ref ) ;
      for( size_t in=0; in<ELM->nb_nodes() ; ++in )
      {
         double bf = ELM->N_local( in, pt_ref ) ;
         for( size_t ic=0 ; ic<result.size() ; ++ic )
         {
            result(ic) += VVALS(in,ic) * bf ;
         }
      }
      pt_ref->destroy() ;
   }

   PEL_CHECK_POST( result.size() == nb_components() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_ValuesOnMeshing_ERROR:: n0( std::string const& name,
                                 std::string const& m_v,
                                 std::string const& col )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "Initialization of field \"" << name << "\" :" << endl
        << "   there is no " << m_v << " of color \"" << col << "\"" ; 
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_ValuesOnMeshing_ERROR:: n1( std::string const& name,
                                 std::string const& col1,
                                 std::string const& col2 )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "Initialization of field \"" << name << "\" :" << endl
	<< "   initialization ambiguity due the occurence of the" << endl
        << "   two matching colors \"" << col1 << "\" and \"" << col2 << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_ValuesOnMeshing_ERROR:: n2_1( std::string const& name,
                                  GE_Mpolyhedron const* mesh_poly,
                                  GE_Color const* mesh_color )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "Initialization of field \"" << name << "\" :" << endl ;
   mesg << "   when calculating the initialization value," << endl ;
   mesg << "   unable to find a formula for the following mesh" << endl ;
   mesh_poly->print( mesg, 6 ) ;
   mesg << "   whose color is" << endl ;
   mesh_color->print( mesg, 6 ) ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_ValuesOnMeshing_ERROR:: n2_2( std::string const& name,
                                  GE_Point const* pt,
                                  GE_Color const* color )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "Initialization of field \"" << name << "\" :" << endl ;
   mesg << "   when calculating the initialization value," << endl ;
   mesg << "   unable to find a formula for the following vertex" << endl ;
   pt->print( mesg, 6 ) ; mesg << endl ;
   mesg << "   whose color is" << endl ;
   color->print( mesg, 6 ) ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_ValuesOnMeshing_ERROR:: n3( std::string const& name,
                                std::string const& col1, size_t nbc1,
                                std::string const& col2, size_t nbc2 )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "Initialization of field \"" << name << "\" :" << endl
        << "   the value associated color \"" << col1
        << "\" has " << nbc1 << " component(s)" << endl
	<< "   the value associated color \"" << col2
        << "\" has " << nbc2 << " component(s)" << endl ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_ValuesOnMeshing_ERROR:: n4( std::string const& name,
                                stringVector const& undef_vars )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "Initialization of field \"" << name << "\" :" << endl
        << "   invalid initialization value" << endl
        << "   the given expressions cannot be evaluated" << endl ;
   if( undef_vars.size() > 0 )
   {
      mesg << "   undefined variable(s): " << endl ;
      for( size_t i=0 ; i<undef_vars.size() ; ++i )
      {
         mesg << "      - \"" << undef_vars(i) << "\"" << endl ;
      }
   }
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_ValuesOnMeshing_ERROR:: n5( std::string const& name,
                                GE_Point const* pt )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "Initialization of field \"" << name << "\" :" << endl
        << "   node found outside the stored grid" << endl
        << "   without \"out_of_stored_grid_value\" defined" << endl ;
   mesg << "Node coordinates : " ;
   pt->print( mesg, 0 ) ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_ValuesOnMeshing_ERROR:: n6( std::string const& name )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "Initialization of field \"" << name << "\" :" << endl
        << "   type \"defined_by_PDE_ResultReader\" without PDE_ResultReader defined" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_ValuesOnMeshing_ERROR:: n7( std::string const& name )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "Initialization of field \"" << name << "\" :" << endl
        << "   field not found in the PDE_ResultReader" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_ValuesOnMeshing_ERROR:: n8( std::string const& name,
                                size_t current, size_t expected )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "Initialization of field \"" << name << "\" :" << endl
        << "   bad dimension found in the PDE_ResultReader" << endl
        << "      expected : " << expected << endl
        << "      read     : " << current ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PDE_ValuesOnMeshing_ERROR:: n9( std::string const& name,
                                size_t current, size_t expected )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "Initialization of field \"" << name << "\" :" << endl
        << "   bad dimensions for the formula used" << endl
        << "   for initialization" << endl
        << "      expected : " << expected << endl
        << "      read     : " << current ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}
