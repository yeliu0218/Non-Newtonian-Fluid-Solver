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

#include <PDE_LocalFE.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_Vector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QuadratureRule.hh>
#include <GE_QRprovider.hh>

#include <PDE_CrossProcessNodeNumbering.hh>
#include <PDE_BFvalues.hh>
#include <PDE_LinkDOF2Unknown.hh>

#include <ios>
#include <iostream>
#include <iomanip>
#include <sstream>

using std::cout ; using std::endl ;
using std::ios_base ;
using std::setprecision ; using std::setw ;
using std::string ;
using std::ostringstream ;

int const PDE_LocalFE::node = 0 ;
int const PDE_LocalFE::N    = PDE_BFvalues::N   ;
int const PDE_LocalFE::dN   = PDE_BFvalues::dN  ;
int const PDE_LocalFE::d2N  = PDE_BFvalues::d2N ;

//----------------------------------------------------------------------
PDE_LocalFE:: PDE_LocalFE( PEL_Object* a_owner, size_t aNbSpDims )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , NB_FIELDs( 0 )
   , GLOB_2_iF( 0 )
   , F_DERIS( 0 )
   , MAX_NB_COMPs( 0 )
   , NB_SP_DIMS( aNbSpDims )
   , EXCLUDE_COLORS( PEL_Vector::create( this, 0 ) )
   , HESSIAN( false )
{
   PEL_LABEL( "PDE_LocalFE:: PDE_LocalFE" ) ;

   PEL_ASSERT( N   == PDE_BFvalues::N   ) ;
   PEL_ASSERT( dN  == PDE_BFvalues::dN  ) ;
   PEL_ASSERT( d2N == PDE_BFvalues::d2N ) ;

   PEL_CHECK( nb_space_dimensions()>=1 && nb_space_dimensions()<=3 ) ;
}

//----------------------------------------------------------------------
PDE_LocalFE:: ~PDE_LocalFE( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
size_t
PDE_LocalFE:: nb_space_dimensions( void ) const
//----------------------------------------------------------------------
{
   return( NB_SP_DIMS ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFE:: require_field_calculation( PDE_DiscreteField const* ff,
                                         int order )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE:: require_field_calculation" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;
   PEL_CHECK_PRE( !is_valid() ) ;
   PEL_CHECK_PRE( order == PDE_LocalFE::node ||
                  order == PDE_LocalFE::N    ||
                  order == PDE_LocalFE::dN   ||
                  order == PDE_LocalFE::d2N ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( order != PDE_LocalFE::node ) JACOBIAN=true ;
   if( order == PDE_LocalFE::d2N ) HESSIAN=true ;

   size_t const i = ff->id_number() ;

   if( i+1 > GLOB_2_iF.size() )
   {
      GLOB_2_iF.resize( i+1, PEL::bad_index() ) ;
   }
   if( GLOB_2_iF[i] == PEL::bad_index() )
   {
      // New field:
      GLOB_2_iF[i] = NB_FIELDs ;
      F_DERIS.push_back( order ) ;
      FIELDS.push_back( ff ) ;
      size_t const nbc = ff->nb_components() ;
      if( nbc > MAX_NB_COMPs )
      {
         MAX_NB_COMPs = nbc ;
      }
      NB_FIELDs++ ;
   }
   else
   {
      F_DERIS[GLOB_2_iF[i]] |= order ;
   }

   PEL_CHECK( FIELDS[GLOB_2_iF[i]]==ff ) ;


   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( field_is_handled( ff ) ) ;
   PEL_CHECK_POST( field_calculation_is_handled( ff, order ) ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFE:: field_calculation_is_handled( PDE_DiscreteField const* ff,
                                            int order ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE:: field_calculation_is_handled" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;
   PEL_CHECK_PRE( order == PDE_LocalFE::node ||
                  order == PDE_LocalFE::N    ||
                  order == PDE_LocalFE::dN   ||
                  order == PDE_LocalFE::d2N ) ;
   PEL_CHECK_INV( invariant() ) ;

   bool result = field_is_handled( ff ) ;
   if( result && order != PDE_LocalFE::node )
   {
      size_t const i = ff->id_number() ;
      result = ( F_DERIS[GLOB_2_iF[i]] & order ) ;
   }

   PEL_CHECK_POST( IMPLIES( field_calculation_is_handled( ff, order ),
                            field_is_handled( ff ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFE:: field_is_handled( PDE_DiscreteField const* ff ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE:: field_is_handled" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t const i = ff->id_number() ;
   return( i<GLOB_2_iF.size() && GLOB_2_iF[i] != PEL::bad_index() ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFE:: nb_handled_fields( void ) const
//----------------------------------------------------------------------
{
   return( NB_FIELDs ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LocalFE:: field_local_index( PDE_DiscreteField const* ff ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE:: field_local_index" ) ;
   PEL_CHECK( ff != 0 ) ;
   PEL_CHECK( field_is_handled( ff ) ) ;
   PEL_CHECK( GLOB_2_iF[ff->id_number()]<nb_handled_fields() ) ;
   return( GLOB_2_iF[ff->id_number()] ) ;
}

//----------------------------------------------------------------------
PDE_DiscreteField const*
PDE_LocalFE:: handled_field( size_t iF ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE:: handled_field" ) ;
   PEL_CHECK_PRE( iF < nb_handled_fields() ) ;

   PDE_DiscreteField const* result = FIELDS[iF] ;

   if( !field_is_handled( result ) )
   {
      std::cout << nb_handled_fields() << " " << iF << " " << result->id_number() << " " << GLOB_2_iF.size() << " " << GLOB_2_iF[result->id_number()] << std::endl ;
   }
   PEL_CHECK( field_local_index( result ) == iF ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( field_is_handled( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFE:: exclude_color( GE_Color const* a_color )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE:: exclude_color" ) ;
   PEL_CHECK_PRE( !is_valid() ) ;
   PEL_CHECK_PRE( a_color!=0 ) ;
   PEL_CHECK_PRE( !is_excluded( a_color ) ) ;

   EXCLUDE_COLORS->append( const_cast<GE_Color*>( a_color ) ) ;

   PEL_CHECK_POST( is_excluded( a_color ) ) ;
}

//----------------------------------------------------------------------
void
PDE_LocalFE:: include_color( GE_Color const* a_color )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE:: include_color" ) ;
   PEL_CHECK_PRE( !is_valid() ) ;
   PEL_CHECK_PRE( a_color!=0 ) ;
   PEL_CHECK_PRE( is_excluded( a_color ) ) ;

   EXCLUDE_COLORS->remove_at( EXCLUDE_COLORS->index_of( a_color ) ) ;

   PEL_CHECK_POST( !is_excluded( a_color ) ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFE:: is_excluded( GE_Color const* a_color ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE:: is_excluded" ) ;
   PEL_CHECK_PRE( a_color!=0 ) ;

   bool result = EXCLUDE_COLORS->has( a_color ) ;

   return( result ) ;
}


//----------------------------------------------------------------------
void
PDE_LocalFE:: print( std::ostream& os, size_t indent_width  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE:: print" ) ;
   PEL_CHECK_PRE( os ) ;

   string space( indent_width, ' ' ) ;
   os << space << "Handled fields" << endl ;
   print_handled_fields( os, indent_width+3 ) ;
   if( is_valid() )
   {
      os << space << "Mesh-iterator movement" << endl ;
      print_current_mesh( os, indent_width+3 ) ;
      os << space << "Local discretization of fields" << endl ;
      print_local_discretization_of_fields( os, indent_width+3 ) ;
      if( field( row ) != 0 )
      {
         os << space << "Local discrete variational problem" << endl ;
         print_local_discrete_variational_problem( os, indent_width+3 ) ;
      }
      if( valid_IP() )
      {
         os << space << "IP-iterator movement" << endl ;
         print_current_IP( os, indent_width+3 ) ;
         print_values_at_current_IP( os, indent_width+3 ) ;
      }
   }
}

//----------------------------------------------------------------------
void
PDE_LocalFE:: print_handled_fields( std::ostream& os,
                                    size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE:: print_handled_fields" ) ;

   string space( indent_width, ' ' ) ;

   for( size_t iF=0 ; iF<nb_handled_fields() ; ++iF )
   {
      PDE_DiscreteField const* ff = handled_field( iF ) ;
      os << space << "\"" << ff->name() << "\"" ;
      if( field_calculation_is_handled( ff, PDE_LocalFE::N ) )   os << " N" ;
      if( field_calculation_is_handled( ff, PDE_LocalFE::dN ) )  os << " dN" ;
      if( field_calculation_is_handled( ff, PDE_LocalFE::d2N ) ) os << " d2N" ;
      os << endl ;
   }
}

//----------------------------------------------------------------------
void
PDE_LocalFE:: print_local_discretization_of_fields(
                                               std::ostream& os,
                                               size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE:: print_local_discretization_of_fields" ) ;
   PEL_CHECK_PRE( os ) ;
   PEL_CHECK_PRE( print_local_discretization_of_fields_PRE( os, indent_width ) ) ;

   std::string space( indent_width, ' ' ) ;
   std::string more( 3, ' ' ) ;
   
   size_t nb_ranks = PEL_Exec::communicator()->nb_ranks() ;

   for( size_t iF=0 ; iF<nb_handled_fields() ; ++iF )
   {
      PDE_DiscreteField const* ff = handled_field( iF ) ;
      os << space << "\"" << ff->name() << "\"  "
         << nb_local_nodes( ff ) << " local nodes" << endl ;
      for( size_t il=0 ; il<nb_local_nodes(ff) ; ++il )
      {
         os << space << more ;
         size_t nn = global_node( ff, il ) ;
         os << setw( 3 ) << nn ;
         GE_Point const* pt = local_node_location( ff, il ) ;
         pt->print( os, 2 ) ;
         os << " imp : " ;
         for( size_t ic=0 ; ic<ff->nb_components() ; ++ic )
         {
            if( ic == 0 ) os << "(" ;
            else os << "," ;
            if( ff->DOF_has_imposed_value( nn, ic ) )
               os << "y" ;
            else
               os << "n" ;
            if( ic==ff->nb_components()-1 ) os << ")" ;
         }
         if( nb_ranks > 1 )
         {
            PDE_CrossProcessNodeNumbering const* nbr =
                                          ff->cross_process_numbering() ;
            os << " glob:" ;
            os << setw( 3 ) << nbr->global_node_index( nn ) ;
            if( nbr->current_process_handles_node( nn ) )
            {
               os << " owned" ;
            }
         }
         os << endl ;
      }
   }
}

//----------------------------------------------------------------------
void
PDE_LocalFE:: print_local_discrete_variational_problem(
                                                std::ostream& os,
                                                size_t indent_width  ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE:: print_local_discrete_variational_problem" ) ;
   PEL_CHECK_PRE( os ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   PEL_CHECK_PRE( field( row ) != 0 ) ;
   PEL_CHECK_PRE( field( col ) != 0 ) ;

   string space( indent_width, ' ' ) ;

   os << indent_width << "row field : " << field(row)->name() ;
   size_t_vector const& rc = row_field_node_connectivity() ;
   os << indent_width << "connectivity : " ;
   for( size_t i=0 ; i<rc.size() ; ++i )
   {
      os << rc(i) << " " << endl ;
   }
   os << indent_width << "col field : " << field(col)->name() ;
   size_t_vector const& cc = col_field_node_connectivity() ;
   os << indent_width << "connectivity : " ;
   for( size_t i=0 ; i<cc.size() ; ++i )
   {
      os << cc(i) << " " << endl ;
   }
}

//----------------------------------------------------------------------
int
PDE_LocalFE:: handled_field_deri( size_t iF ) const
//----------------------------------------------------------------------
{
   PEL_CHECK( iF < nb_handled_fields() ) ;

   return( F_DERIS[iF] ) ;
}

//------------------------------------------------------------------------
size_t
PDE_LocalFE:: handled_fields_max_nb_comps( void ) const
//------------------------------------------------------------------------
{
   return( MAX_NB_COMPs ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: jacobian_of_mapping_required( void ) const
//------------------------------------------------------------------------
{
   return( JACOBIAN ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: hessian_of_mapping_required( void ) const
//------------------------------------------------------------------------
{
   return( HESSIAN ) ;
}

//------------------------------------------------------------------------
void
PDE_LocalFE:: raise_no_discretization_error( PDE_DiscreteField const* ff )
//------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE:: raise_no_discretization_error" ) ;
   PEL_CHECK( is_valid() ) ;

   ostringstream mesg ;
   mesg << type_name() << " : " << endl ;
   mesg << "   the field \"" << ff->name()
        << "\" is handled" << endl ;
   mesg << "   but its discretization is unknown on the current mesh."
        << endl ;
   polyhedron()->print( mesg, 6 ) ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: go_i_th_PRE( size_t an_id_mesh ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( an_id_mesh == PEL::bad_index() || an_id_mesh<nb_meshes() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: go_i_th_POST( size_t an_id_mesh ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( an_id_mesh == PEL::bad_index(), !is_valid() ) ) ;
   PEL_ASSERT( IMPLIES( is_valid(), mesh_id() == an_id_mesh ) ) ;
   PEL_ASSERT( IMPLIES( is_valid(), field(row)==0 ) ) ;
   PEL_ASSERT( IMPLIES( is_valid(), field(col)==0 ) ) ;
   PEL_ASSERT( IMPLIES( is_valid(), calculation_point()==0 ) ) ;
   PEL_ASSERT( IMPLIES( is_valid(), !valid_IP() ) ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: start_POST( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( is_valid(), field(row)==0 ) ) ;
   PEL_ASSERT( IMPLIES( is_valid(), field(col)==0 ) ) ;
   PEL_ASSERT( IMPLIES( is_valid(), calculation_point()==0 ) ) ;
   PEL_ASSERT( IMPLIES( is_valid(), !valid_IP() ) ) ;
   PEL_ASSERT( calculation_point() == 0 ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: go_next_PRE( void ) const
//-----------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: go_next_POST( void ) const
//-----------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( is_valid(), field(row)==0 ) ) ;
   PEL_ASSERT( IMPLIES( is_valid(), field(col)==0 ) ) ;
   PEL_ASSERT( IMPLIES( is_valid(), calculation_point()==0 ) ) ;
   PEL_ASSERT( IMPLIES( is_valid(), !valid_IP() ) ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: polyhedron_PRE( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: refinement_level_PRE( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: polyhedron_POST( GE_Mpolyhedron const* result ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->nb_space_dimensions() == nb_space_dimensions() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: color_PRE( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: color_POST( GE_Color const* result ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( result!=0 ) ;
   PEL_ASSERT( !is_excluded( result ) ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: mesh_id_PRE( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: mesh_id_POST( size_t result ) const
//------------------------------------------------------------------------
{
//????   PEL_ASSERT( result < nb_meshes() ) ;
   return( true ) ;
}


//------------------------------------------------------------------------
bool
PDE_LocalFE:: nb_local_nodes_PRE( PDE_DiscreteField const* ff ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( ff != 0 ) ;
   PEL_ASSERT( field_is_handled( ff ) ) ;

   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: nb_local_nodes_POST( size_t result,
                                   PDE_DiscreteField const* ff ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( ff==field(row), result==nb_basis_functions(row) ) ) ;
   PEL_ASSERT( IMPLIES( ff==field(col), result==nb_basis_functions(col) ) ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFE:: local_node_is_in_mesh_PRE( PDE_DiscreteField const* ff,
                                         size_t i ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( ff != 0 ) ;
   PEL_ASSERT( field_is_handled( ff ) ) ;
   PEL_ASSERT( i < nb_local_nodes( ff ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFE:: local_node_is_in_mesh_POST( bool result,
                                          PDE_DiscreteField const* ff,
                                          size_t i ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( FORMAL(
      EQUIVALENT( result,
                  polyhedron()->contains( local_node_location( ff, i ) ) ))) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFE:: local_node_location_PRE( PDE_DiscreteField const* ff,
                                       size_t i ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( ff != 0 ) ;
   PEL_ASSERT( field_is_handled( ff ) ) ;
   PEL_ASSERT( i < nb_local_nodes( ff ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFE:: local_node_location_POST( GE_Point const* result,
                                        PDE_DiscreteField const* ff,
                                        size_t i ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->is_under_ownership_of( this ) ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: global_node_PRE( PDE_DiscreteField const* ff,
                               size_t i ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( ff != 0 ) ;
   PEL_ASSERT( field_is_handled( ff ) ) ;
   PEL_ASSERT( i < nb_local_nodes( ff ) ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: global_node_POST( size_t result,
                                PDE_DiscreteField const* ff,
                                size_t i ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( result < ff->nb_nodes() ) ;
   PEL_ASSERT( ff->node_is_active( result ) ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: node_refinement_level_PRE( PDE_DiscreteField const* ff,
                                         size_t i ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( ff != 0 ) ;
   PEL_ASSERT( field_is_handled( ff ) ) ;
   PEL_ASSERT( i < nb_local_nodes( ff ) ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: node_refinement_level_POST( size_t result,
                                          PDE_DiscreteField const * ff,
                                          size_t i ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( result <= refinement_level() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: set_calculation_point_PRE( GE_Point const* pt ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( pt != 0 ) ;
   PEL_ASSERT( pt->nb_coordinates() == nb_space_dimensions() ) ;
   PEL_ASSERT( FORMAL( polyhedron()->contains( pt ) ) ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: set_calculation_point_POST( GE_Point const* pt ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( calculation_point() != 0 ) ;
   PEL_ASSERT( calculation_point() != pt ) ;
   PEL_ASSERT( calculation_point()->nb_coordinates() ==
               pt->nb_coordinates() ) ;
   PEL_ASSERT(
      FORALL( ( size_t i=0 ; i<pt->nb_coordinates() ; ++i ),
         calculation_point()->coordinate( i ) == pt->coordinate( i ) ) ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: calculation_point_POST( GE_Point const* result ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( IMPLIES( result!=0, result->is_under_ownership_of( this ) ) ) ;
   PEL_ASSERT( FORMAL(
               IMPLIES( result!=0, polyhedron()->contains( result ) ) ) ) ;
   PEL_ASSERT( IMPLIES( result!=0,
                        result->nb_coordinates()==nb_space_dimensions() ) ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: value_at_pt_PRE( PDE_DiscreteField const* ff,
                               size_t level,
                               size_t ic ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( calculation_point() != 0 ) ;
   PEL_ASSERT( ff != 0 ) ;
   PEL_ASSERT( field_calculation_is_handled( ff, N ) ) ;
   PEL_ASSERT( nb_local_nodes( ff ) != 0 ) ;
   PEL_ASSERT( level < ff->storage_depth() ) ;
   PEL_ASSERT( ic < ff->nb_components() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: gradient_at_pt_PRE( PDE_DiscreteField const* ff,
                                  size_t level,
                                  size_t a,
                                  size_t ic ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( calculation_point() != 0 ) ;
   PEL_ASSERT( ff != 0 ) ;
   PEL_ASSERT( field_calculation_is_handled( ff, dN ) ) ;
   PEL_ASSERT( nb_local_nodes( ff ) != 0 ) ;
   PEL_ASSERT( level < ff->storage_depth() ) ;
   PEL_ASSERT( ic < ff->nb_components() ) ;
   PEL_ASSERT( a < nb_space_dimensions() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: hessian_at_pt_PRE( PDE_DiscreteField const* ff,
                                 size_t level,
                                 size_t a,
                                 size_t b,
                                 size_t ic ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( calculation_point() != 0 ) ;
   PEL_ASSERT( ff != 0 ) ;
   PEL_ASSERT( field_calculation_is_handled( ff, d2N ) ) ;
   PEL_ASSERT( nb_local_nodes( ff ) != 0 ) ;
   PEL_ASSERT( level < ff->storage_depth() ) ;
   PEL_ASSERT( ic < ff->nb_components() ) ;
   PEL_ASSERT( a < nb_space_dimensions() ) ;
   PEL_ASSERT( b < nb_space_dimensions() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: N_at_pt_PRE( PDE_DiscreteField const* ff,
                           size_t i ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( calculation_point() != 0 ) ;
   PEL_ASSERT( ff != 0 ) ;
   PEL_ASSERT( field_calculation_is_handled( ff, N ) ) ;
   PEL_ASSERT( nb_local_nodes( ff ) != 0 ) ;
   PEL_ASSERT( i < nb_local_nodes( ff )  ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: dN_at_pt_PRE( PDE_DiscreteField const* ff,
                            size_t i,
                            size_t a ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( calculation_point() != 0 ) ;
   PEL_ASSERT( ff != 0 ) ;
   PEL_ASSERT( field_calculation_is_handled( ff, dN ) ) ;
   PEL_ASSERT( nb_local_nodes( ff ) != 0 ) ;
   PEL_ASSERT( i < nb_local_nodes( ff )  ) ;
   PEL_ASSERT( a < nb_space_dimensions() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: d2N_at_pt_PRE( PDE_DiscreteField const* ff,
                             size_t i,
                             size_t a,
                             size_t b ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( calculation_point() != 0 ) ;
   PEL_ASSERT( ff != 0 ) ;
   PEL_ASSERT( field_calculation_is_handled( ff, d2N ) ) ;
   PEL_ASSERT( nb_local_nodes( ff ) != 0 ) ;
   PEL_ASSERT( i < nb_local_nodes( ff )  ) ;
   PEL_ASSERT( a < nb_space_dimensions() ) ;
   PEL_ASSERT( b < nb_space_dimensions() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: set_row_and_col_fields_PRE(
                              PDE_DiscreteField const* row_field,
                              PDE_DiscreteField const* col_field) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( row_field != 0 ) ;
   PEL_ASSERT( col_field != 0 ) ;
   PEL_ASSERT( field_is_handled( row_field ) ) ;
   PEL_ASSERT( nb_local_nodes( row_field ) != 0 ) ;
   PEL_ASSERT( field_is_handled( col_field ) ) ;
   PEL_ASSERT( nb_local_nodes( col_field ) != 0 ) ;
   PEL_ASSERT( !valid_IP() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: set_row_and_col_fields_POST(
                              PDE_DiscreteField const* row_field,
                              PDE_DiscreteField const* col_field) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( field( row ) == row_field ) ;
   PEL_ASSERT( field( col ) == col_field ) ;
   PEL_ASSERT( !valid_IP() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: row_field_node_connectivity_PRE( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: row_field_node_connectivity_POST(
                                    size_t_vector const& result ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( result.size() == nb_basis_functions( row ) ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: col_field_node_connectivity_PRE( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: col_field_node_connectivity_POST(
                                    size_t_vector const& result ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( result.size() == nb_basis_functions( col ) ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: nb_basis_functions_PRE( field_id sf ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( field( row ) != 0 ) ;
   PEL_ASSERT( field( col ) != 0 ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: nb_basis_functions_POST( size_t result, field_id sf ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( result == nb_local_nodes( field( sf ) ) ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: start_IP_iterator_PRE( GE_QRprovider const* qrp ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( qrp != 0 ) ;
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: start_IP_iterator_POST( GE_QRprovider const* qrp ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( valid_IP() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: go_next_IP_PRE( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( valid_IP() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: valid_IP_PRE( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE::  weight_of_IP_PRE( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( valid_IP() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: coordinates_of_IP_PRE( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( valid_IP() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: coordinates_of_IP_POST( GE_Point const* result ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( result!=0 ) ;
   PEL_ASSERT( result->is_under_ownership_of( this ) ) ;
   PEL_ASSERT( result->nb_coordinates()==nb_space_dimensions() ) ;
   PEL_ASSERT( FORMAL( polyhedron()->contains( result ) ) ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: value_at_IP_PRE( PDE_DiscreteField const* ff,
                               size_t level,
                               size_t ic ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( valid_IP() ) ;
   PEL_ASSERT( ff != 0 ) ;
   PEL_ASSERT( field_calculation_is_handled( ff, N ) ) ;
   PEL_ASSERT( nb_local_nodes( ff ) != 0 ) ;
   PEL_ASSERT( level < ff->storage_depth() ) ;
   PEL_ASSERT( ic < ff->nb_components() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: gradient_at_IP_PRE( PDE_DiscreteField const* ff,
                                  size_t level,
                                  size_t a,
                                  size_t ic ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( valid_IP() ) ;
   PEL_ASSERT( ff != 0 ) ;
   PEL_ASSERT( field_calculation_is_handled( ff, dN ) ) ;
   PEL_ASSERT( nb_local_nodes( ff ) != 0 ) ;
   PEL_ASSERT( level < ff->storage_depth() ) ;
   PEL_ASSERT( ic < ff->nb_components() ) ;
   PEL_ASSERT( a < nb_space_dimensions() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: hessian_at_IP_PRE( PDE_DiscreteField const* ff,
                                 size_t level,
                                 size_t a,
                                 size_t b,
                                 size_t ic ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( valid_IP() ) ;
   PEL_ASSERT( ff != 0 ) ;
   PEL_ASSERT( field_calculation_is_handled( ff, d2N ) ) ;
   PEL_ASSERT( nb_local_nodes( ff ) != 0 ) ;
   PEL_ASSERT( level < ff->storage_depth() ) ;
   PEL_ASSERT( ic < ff->nb_components() ) ;
   PEL_ASSERT( a < nb_space_dimensions() ) ;
   PEL_ASSERT( b < nb_space_dimensions() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: N_at_IP_PRE( field_id sf, size_t i ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( valid_IP() ) ;
   PEL_ASSERT( field( row ) != 0 ) ;
   PEL_ASSERT( field( col ) != 0 ) ;
   PEL_ASSERT( field_calculation_is_handled( field( sf ), N ) ) ;
   PEL_ASSERT( i < nb_basis_functions( sf ) ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: dN_at_IP_PRE( field_id sf, size_t i, size_t a ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( valid_IP() ) ;
   PEL_ASSERT( field( row ) != 0 ) ;
   PEL_ASSERT( field( col ) != 0 ) ;
   PEL_ASSERT( field_calculation_is_handled( field( sf ), dN ) ) ;
   PEL_ASSERT( i< nb_basis_functions( sf ) ) ;
   PEL_ASSERT( a < nb_space_dimensions() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: d2N_at_IP_PRE( field_id sf, size_t i, size_t a, size_t b ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( valid_IP() ) ;
   PEL_ASSERT( field( row ) != 0 ) ;
   PEL_ASSERT( field( col ) != 0 ) ;
   PEL_ASSERT( field_calculation_is_handled( field( sf ), d2N ) ) ;
   PEL_ASSERT( i < nb_basis_functions( sf ) ) ;
   PEL_ASSERT( a < nb_space_dimensions() ) ;
   PEL_ASSERT( b < nb_space_dimensions() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFE:: Ns_at_IP_PRE( field_id sf ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( valid_IP() ) ;
   PEL_ASSERT( field( row ) != 0 ) ;
   PEL_ASSERT( field( col ) != 0 ) ;
   PEL_ASSERT( field_calculation_is_handled( field( sf ), N ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFE:: Ns_at_IP_POST( doubleVector const& result, field_id sf ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT(
      FORALL( (size_t i=0 ; i<nb_basis_functions(sf) ; i++ ),
              result(i )== N_at_IP(sf,i) ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFE:: dNs_at_IP_PRE( field_id sf ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( valid_IP() ) ;
   PEL_ASSERT( field( row ) != 0 ) ;
   PEL_ASSERT( field( col ) != 0 ) ;
   PEL_ASSERT( field_calculation_is_handled( field( sf ), dN ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFE:: dNs_at_IP_POST( doubleArray2D const& result, field_id sf ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT(
      FORALL( (size_t i=0 ; i<nb_basis_functions(sf) ; i++ ),
              FORALL( (size_t a=0 ; a<nb_space_dimensions() ; a++ ),
                      result(i,a) == dN_at_IP(sf,i,a) ) ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFE:: d2Ns_at_IP_PRE( field_id sf ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( valid_IP() ) ;
   PEL_ASSERT( field( row ) != 0 ) ;
   PEL_ASSERT( field( col ) != 0 ) ;
   PEL_ASSERT( field_calculation_is_handled( field( sf ), d2N ) ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
PDE_LocalFE:: d2Ns_at_IP_POST( doubleArray3D const& result, field_id sf ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT(
      FORALL( (size_t i=0 ; i<nb_basis_functions(sf) ; i++ ),
              FORALL( (size_t a=0 ; a<nb_space_dimensions() ; a++ ),
                      FORALL( (size_t b=0 ; b<nb_space_dimensions() ; b++ ),
                              result(i,a,b) == d2N_at_IP(sf,i,a,b) ) ) ) ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: mask_value_at_IP_PRE( PDE_DiscreteField const* ff,
                                    size_t level,
                                    double value,
                                    size_t ic ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( valid_IP() ) ;
   PEL_ASSERT( level < ff->storage_depth() ) ;
   PEL_ASSERT( field_calculation_is_handled( ff, N ) ) ;
   PEL_ASSERT( nb_local_nodes( ff ) != 0 ) ;
   PEL_ASSERT( ic < ff->nb_components() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE::  mask_value_at_IP_POST( PDE_DiscreteField const* ff,
                                      size_t level,
                                      double value,
                                      size_t ic ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( value_at_IP(ff,level,ic) == value ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: print_current_mesh_PRE( std::ostream& os,
                                      size_t indent_width ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( os ) ;
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: print_local_discretization_of_fields_PRE(
                                             std::ostream& os,
                                             size_t indent_width ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( os ) ;
   PEL_ASSERT( is_valid() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: print_current_IP_PRE( std::ostream& os,
                                    size_t indent_width ) const

//------------------------------------------------------------------------
{
   PEL_ASSERT( os ) ;
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( valid_IP() ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: print_values_at_current_IP_PRE( std::ostream& os,
                                              size_t indent_width ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( os ) ;
   PEL_ASSERT( is_valid() ) ;
   PEL_ASSERT( valid_IP() ) ;
   PEL_ASSERT( field( row ) != 0 ) ;
   PEL_ASSERT( field( col ) != 0 ) ;
   return( true ) ;
}

//------------------------------------------------------------------------
bool
PDE_LocalFE:: invariant( void ) const
//------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
//   PEL_ASSERT( IMPLIES( is_valid(), mesh_id()<nb_meshes() ) ) ;
   PEL_ASSERT( IMPLIES( valid_IP(), is_valid() ) ) ;
//   PEL_ASSERT( IMPLIES( !is_valid(), calculation_point()==0 ) ) ;
   return( true ) ;
}
