#include <PDE_LocalFE_DEMO1.hh>

#include <PEL_assertions.hh>
#include <PEL_ModuleExplorer.hh>
#include <size_t_vector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_QRprovider.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <fstream>
#include <ios>
#include <iostream>

using std::endl ;
using std::ios_base ;
using std::ofstream ;

PDE_LocalFE_DEMO1 const* 
PDE_LocalFE_DEMO1::PROTOTYPE = new PDE_LocalFE_DEMO1() ;

//----------------------------------------------------------------------------
PDE_LocalFE_DEMO1:: PDE_LocalFE_DEMO1( void )
//----------------------------------------------------------------------------
   : PEL_Application( "PDE_LocalFE_DEMO1" )
{
}

//----------------------------------------------------------------------------
PDE_LocalFE_DEMO1*
PDE_LocalFE_DEMO1:: create_replica( PEL_Object* a_owner, 
                                 PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE_DEMO1:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PDE_LocalFE_DEMO1* result = new PDE_LocalFE_DEMO1( a_owner, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PDE_LocalFE_DEMO1:: PDE_LocalFE_DEMO1( PEL_Object* a_owner, 
                                       PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
   : PEL_Application( a_owner, exp )
   , DOM( 0 )
   , OUT_NAME( exp->string_data( "trace_file" ) )
{
   PEL_ModuleExplorer* e = 
                    exp->create_subexplorer( 0, "PDE_DomainAndFields" ) ;
   DOM = PDE_DomainAndFields::create( this, e ) ;
   e->destroy() ;
}

//----------------------------------------------------------------------------
PDE_LocalFE_DEMO1:: ~PDE_LocalFE_DEMO1( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
PDE_LocalFE_DEMO1:: run( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LocalFE_DEMO1:: run" ) ;

   ofstream ofs( OUT_NAME.c_str() ) ;

   ofs << "do_1()-----------------------" << endl ;
   do_1( ofs ) ;
   ofs << "do_2()-----------------------" << endl ;
   do_2( ofs ) ;
   ofs << "do_3()-----------------------" << endl ;
   do_3( ofs ) ;
   ofs << "do_4()-----------------------" << endl ;
   do_4( ofs ) ;
   ofs << "do_5()-----------------------" << endl ;
   do_5( ofs ) ;
   ofs << "do_6()-----------------------" << endl ;
   do_6( ofs ) ;
   ofs << "do_7()-----------------------" << endl ;
   do_7( ofs ) ;
   ofs << "do_8()-----------------------" << endl ;
   do_8( ofs ) ;

   ofs.close() ;
}
 
//----------------------------------------------------------------------------
void
PDE_LocalFE_DEMO1:: do_1( std::ostream& os )
//----------------------------------------------------------------------------
{
   PDE_LocalFEcell* fe = DOM->create_LocalFEcell( 0 ) ;
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      os << "cell " << fe->mesh_id() << "   " ;
      GE_Mpolyhedron const* poly = fe->polyhedron() ;
      poly->print( os, 0 ) ;
   }
   fe->destroy() ;
}

//----------------------------------------------------------------------------
void
PDE_LocalFE_DEMO1:: do_2( std::ostream& os )
//----------------------------------------------------------------------------
{
   PDE_LocalFEbound* fe = DOM->create_LocalFEbound( 0 ) ;
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      os << "bound " << fe->mesh_id() << "   " ;
      GE_Mpolyhedron const* poly = fe->polyhedron() ;
      poly->print( os, 0 ) ;
   }
   fe->destroy() ;
}

//----------------------------------------------------------------------------
void
PDE_LocalFE_DEMO1:: do_3( std::ostream& os )
//----------------------------------------------------------------------------
{
   PDE_SetOfDiscreteFields const* sdf = DOM->set_of_discrete_fields() ;
   PDE_DiscreteField const* uu = sdf->item( "temp" ) ;
   PDE_LocalFEcell* fe = DOM->create_LocalFEcell( 0 ) ;
   fe->require_field_calculation( uu, PDE_LocalFE::N ) ;
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      fe->set_row_and_col_fields( uu, uu ) ;
      
      size_t_vector const& rc = fe->row_field_node_connectivity() ;
      os << "cell " << fe->mesh_id() << " : nodes " ;
      for( size_t i=0 ; i<rc.size() ; ++i )
      {
         os << "  " << rc(i) ;
      }
      os << endl ;
   }
   fe->destroy() ;
}

//----------------------------------------------------------------------------
void
PDE_LocalFE_DEMO1:: do_4( std::ostream& os )
//----------------------------------------------------------------------------
{
   PDE_SetOfDiscreteFields const* sdf = DOM->set_of_discrete_fields() ;
   PDE_DiscreteField const* uu = sdf->item( "temp" ) ;
   PDE_LocalFEbound* fe = DOM->create_LocalFEbound( 0 ) ;
   fe->require_field_calculation( uu, PDE_LocalFE::N ) ;
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      fe->set_row_and_col_fields( uu, uu ) ;

      size_t_vector const& rc = fe->row_field_node_connectivity() ;
      os << "bound " << fe->mesh_id() << " : nodes " ;
      for( size_t i=0 ; i<rc.size() ; ++i )
      {
         os << "  " << rc(i) ;
      }
      os << endl ;
   }
   fe->destroy() ;
}

//----------------------------------------------------------------------------
void
PDE_LocalFE_DEMO1:: do_5( std::ostream& os )
//----------------------------------------------------------------------------
{
   os.setf( ios_base::fixed ) ; os.precision( 2 ) ;

   GE_QRprovider const* qrp = GE_QRprovider::object( "GE_QRprovider_3" ) ;
   PDE_SetOfDiscreteFields const* sdf = DOM->set_of_discrete_fields() ;
   PDE_DiscreteField const* uu = sdf->item( "temp" ) ;
   PDE_LocalFEcell* fe = DOM->create_LocalFEcell( 0 ) ;
   fe->require_field_calculation( uu, PDE_LocalFE::N ) ;
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      fe->set_row_and_col_fields( uu, uu ) ;

      os << "cell " << fe->mesh_id() << endl ;
      fe->start_IP_iterator( qrp ) ;
      for( ; fe->valid_IP() ; fe->go_next_IP() )
      {
         GE_Point const* pt = fe->coordinates_of_IP() ;
         os << pt->coordinate( 0 ) << " "
              << pt->coordinate( 1 ) << endl ;
      }
   }
   fe->destroy() ;
}

//----------------------------------------------------------------------------
void
PDE_LocalFE_DEMO1:: do_6( std::ostream& os )
//----------------------------------------------------------------------------
{
   os.setf( ios_base::fixed ) ; os.precision( 2 ) ;

   GE_QRprovider const* qrp = GE_QRprovider::object( "GE_QRprovider_3" ) ;
   PDE_SetOfDiscreteFields const* sdf = DOM->set_of_discrete_fields() ;
   PDE_DiscreteField const* uu = sdf->item( "temp" ) ;
   PDE_LocalFEbound* fe = DOM->create_LocalFEbound( 0 ) ;
   fe->require_field_calculation( uu, PDE_LocalFE::N ) ;
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      fe->set_row_and_col_fields( uu, uu ) ;

      os << "bound " << fe->mesh_id() << endl ;
      fe->start_IP_iterator( qrp ) ;
      for( ; fe->valid_IP() ; fe->go_next_IP() )
      {
         GE_Point const* pt = fe->coordinates_of_IP() ;
         os << pt->coordinate( 0 ) << " "
              << pt->coordinate( 1 ) << endl ;
      }
   }
   fe->destroy() ;
}

//----------------------------------------------------------------------------
void
PDE_LocalFE_DEMO1:: do_7( std::ostream& os )
//----------------------------------------------------------------------------
{
   os.setf( ios_base::fixed ) ; os.precision( 2 ) ;

   PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
   GE_QRprovider const* qrp = GE_QRprovider::object( "GE_QRprovider_3" ) ;
   PDE_SetOfDiscreteFields const* sdf = DOM->set_of_discrete_fields() ;
   PDE_DiscreteField const* uu = sdf->item( "temp" ) ;
   PDE_LocalFEcell* fe = DOM->create_LocalFEcell( 0 ) ;
   fe->require_field_calculation( uu, PDE_LocalFE::N ) ;
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      fe->set_row_and_col_fields( uu, uu ) ;
      size_t nb_bfs = fe->nb_basis_functions( row ) ;
      os << "cell " << fe->mesh_id() << endl ;
      fe->start_IP_iterator( qrp ) ;
      for( ; fe->valid_IP() ; fe->go_next_IP() )
      {
         for( size_t i=0 ; i<nb_bfs ; ++i )
            os << fe->N_at_IP( row, i ) << " " ;
         os << endl ;
      }
   }
   fe->destroy() ;
}

//----------------------------------------------------------------------------
void
PDE_LocalFE_DEMO1:: do_8( std::ostream& os )
//----------------------------------------------------------------------------
{
   os.setf( ios_base::fixed ) ; os.precision( 2 ) ;

   PDE_LocalFE::field_id const row = PDE_LocalFE::row ;
   GE_QRprovider const* qrp = GE_QRprovider::object( "GE_QRprovider_3" ) ;
   PDE_SetOfDiscreteFields const* sdf = DOM->set_of_discrete_fields() ;
   PDE_DiscreteField const* uu = sdf->item( "temp" ) ;
   PDE_LocalFEbound* fe = DOM->create_LocalFEbound( 0 ) ;
   fe->require_field_calculation( uu, PDE_LocalFE::N ) ;
   for( fe->start() ; fe->is_valid() ; fe->go_next() )
   {
      fe->set_row_and_col_fields( uu, uu ) ;
      size_t nb_bfs = fe->nb_basis_functions( row ) ;
      os << "bound " << fe->mesh_id() << endl ;
      fe->start_IP_iterator( qrp ) ;
      for( ; fe->valid_IP() ; fe->go_next_IP() )
      {
         for( size_t i=0 ; i<nb_bfs ; ++i )
            os << fe->N_at_IP( row, i ) << " " ;
         os << endl ;
      }
   }
   fe->destroy() ;
}

