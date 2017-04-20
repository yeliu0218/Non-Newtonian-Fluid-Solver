#include <FE_InterfaceIndicator.hh>

#include <PEL.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Int.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QRprovider.hh>
#include <GE_Vector.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ReferenceElement.hh>
#include <PDE_SetOfBCs.hh>

#include <iostream>
#include <string>

using std::string ;
using std::cout ; using std::endl ;

FE_InterfaceIndicator const*
FE_InterfaceIndicator:: PROTOTYPE = new FE_InterfaceIndicator() ;

//---------------------------------------------------------------------------
FE_InterfaceIndicator:: FE_InterfaceIndicator( void )
//---------------------------------------------------------------------------
   : PDE_AdaptationIndicator( "FE_InterfaceIndicator" )
   , CELL_INTERF( 0 )
{
}

//---------------------------------------------------------------------------
FE_InterfaceIndicator*
FE_InterfaceIndicator:: create_replica( 
                                  PEL_Object* a_owner,
                                  PDE_DomainAndFields const* dom,
                                  PEL_ModuleExplorer const* exp,
                                  size_t a_verbose_level ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_InterfaceIndicator:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp, a_verbose_level ) ) ;

   FE_InterfaceIndicator* result = 
       new FE_InterfaceIndicator( a_owner, dom, exp, a_verbose_level ) ;

   PEL_CHECK( create_replica_POST( result,
                                   a_owner, dom, exp, a_verbose_level ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
FE_InterfaceIndicator:: FE_InterfaceIndicator( 
                                           PEL_Object* a_owner,
                                           PDE_DomainAndFields const* dom,
                                           PEL_ModuleExplorer const* exp,
                                           size_t a_verbose_level  )
//---------------------------------------------------------------------------
   : PDE_AdaptationIndicator( a_owner, a_verbose_level )
   , CTX( 0 )
   , COORDS( 0 )
   , ITER( 0 )
   , PHF( 0 )
   , cFE( dom->create_LocalFEcell( this ) )
   , QRP( GE_QRprovider::object( 
                           exp->string_data( "quadrature_rule_provider" ) ) )
   , CELL_INTERF( 0 )
   , H_INTERF( exp->double_data( "h_for_interface" ) )
   , BF_MIN_REFI( exp->double_data( "bf_min_refinement" ) )
   , BF_MAX_REFI( exp->double_data( "bf_max_refinement" ) )
   , BF_MIN_UNREFI( exp->double_data( "bf_min_unrefinement" ) )
   , BF_MAX_UNREFI( exp->double_data( "bf_max_unrefinement" ) )
   , ICALL( PEL::bad_index() )
{
   PEL_ContextSimple* c = PEL_ContextSimple::create( this ) ;
   COORDS = PEL_DoubleVector::create( c, doubleVector(0) ) ;
   c->extend( PEL_Variable::object( "DV_X" ), COORDS ) ;
   ITER = PEL_Int::create( c, 0 ) ;
   c->extend( PEL_Variable::object( "IS_ITER" ), ITER ) ;
   CTX = c ;
   PHF = exp->abstract_data( this, "phase_field", CTX ) ;
   if( !PHF->value_can_be_evaluated() )
   {
      PEL_Error::object()->raise_not_evaluable(
         exp, "phase_field", PHF->undefined_variables() ) ;
   }
   if( PHF->data_type()!=PEL_Data::Double )
   {
      PEL_Error::object()->raise_bad_data_type(
            exp, "phase_field", PEL_Data::Double ) ;
   }

   cFE->include_color( GE_Color::halo_color() ) ;

   CELL_INTERF.re_initialize( cFE->nb_meshes() ) ;
}

//---------------------------------------------------------------------------
FE_InterfaceIndicator:: ~FE_InterfaceIndicator( void  )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
FE_InterfaceIndicator:: reset( void )
//---------------------------------------------------------------------------
{
   if( ICALL == PEL::bad_index() )
   {
      ICALL = 0 ;
   }
   else
   {
      ++ICALL ;
   }
}

//---------------------------------------------------------------------------
void
FE_InterfaceIndicator:: build( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_InterfaceIndicator:: build" ) ;

   ITER->set( ICALL ) ;

   CELL_INTERF.set( 0.0 ) ;

   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      size_t id = cFE->mesh_id() ;
      if( id >= CELL_INTERF.size() ) CELL_INTERF.resize( id+1 ) ;
      cFE->start_IP_iterator( QRP ) ;
      double ss = 0.0 ;
      for( ; cFE->valid_IP() ; cFE->go_next_IP() )
      {
         COORDS->set( cFE->coordinates_of_IP()->coordinate_vector() ) ;
         ss += PHF->to_double() * cFE->weight_of_IP() ;
      }
      CELL_INTERF( id ) = ss / cFE->polyhedron()->measure() ;
   }
}

//---------------------------------------------------------------------------
bool
FE_InterfaceIndicator:: to_be_refined( double bf_indicator,
                                       GE_Mpolyhedron const* poly,
                                       PDE_ReferenceElement const* elm,
                                       size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_InterfaceIndicator:: to_be_refined" ) ;
   PEL_CHECK_PRE( to_be_refined_PRE( bf_indicator, poly, elm, local_node ) ) ;

   bool result = false ;

   if( (bf_indicator > BF_MIN_REFI) && (bf_indicator < BF_MAX_REFI) )
   {
      double h = poly->inter_vertices_maximum_distance() ;
      if( h > H_INTERF ) result = true ;
   }

   return( result ) ;
}

//---------------------------------------------------------------------------
bool
FE_InterfaceIndicator:: to_be_unrefined( double bf_indicator,
                                         GE_Mpolyhedron const* poly,
                                         PDE_ReferenceElement const* elm,
                                         size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_InterfaceIndicator:: to_be_unrefined" ) ;
   PEL_CHECK_PRE( to_be_unrefined_PRE( bf_indicator, poly, elm, local_node ) );

   bool result = false ;

   if( (bf_indicator < BF_MIN_UNREFI) || (bf_indicator > BF_MAX_UNREFI) )
   {
      result = true ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------------
double
FE_InterfaceIndicator:: cell_indicator( size_t cell_id ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_InterfaceIndicator:: cell_indicator" ) ;

   cFE->go_i_th( cell_id ) ;
   PEL_ASSERT( cFE->is_valid() ) ;

   double result = CELL_INTERF( cell_id ) ;
   return( result ) ;
}
