#include <MY_MUSCL_DataStructure.hh>

#include <FE.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_Vector.hh>

#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>

PEL_Vector*
MY_MUSCL_DataStructure::OBJECTS = PEL_Vector::create( PEL_Root::object(), 0 ) ;

class MY_OneSideStreamCellNodes : public PEL_Object
{
   public :
      
      MY_OneSideStreamCellNodes( PEL_Object* a_owner,
				 PDE_DomainAndFields const* dom,
				 PDE_DiscreteField const* f,
				 PDE_CursorFEside const* fe ) ;
     ~MY_OneSideStreamCellNodes( void ) ;

      void one_direction_determination( PDE_DomainAndFields const* dom,
                                        PDE_CursorFEside const* fe,
					PDE_LocalFEcell const* cFE,
					PDE_DiscreteField const* field,
                                        size_t idir,
                                        double ini_dist,
					bool& has_node,
					size_t& node,
                                        double& dist ) ;

      size_t UPS_N ;
      double UPS_N_D ;
      size_t  DOWNS_N ;
      double DOWNS_N_D ;
      bool HAS_UP_UPS_N ;
      size_t UP_UPS_N ;
      double UP_UPS_N_D ;
      bool HAS_DOWN_DOWNS_N ;
      size_t DOWN_DOWNS_N ;
      double DOWN_DOWNS_N_D ;
} ;

class MY_OneBoundStreamCellNodes : public PEL_Object
{
   public :
      
      MY_OneBoundStreamCellNodes( PEL_Object* a_owner,
				  PDE_DomainAndFields const* dom,
				  PDE_DiscreteField const* f,
				  PDE_LocalFEbound const* fe ) ;
     ~MY_OneBoundStreamCellNodes( void ) ;

      void one_direction_determination( PDE_DomainAndFields const* dom,
                                        PDE_LocalFEbound const* fe,
					PDE_LocalFEcell const* cFE,
					PDE_DiscreteField const* field,
                                        size_t idir,
                                        double ini_dist,
					bool& has_node,
					size_t& node,
                                        double& dist ) ;

      size_t UPS_N ;
      double UPS_N_D ;
      bool HAS_UP_UPS_N ;
      size_t UP_UPS_N ;
      double UP_UPS_N_D ;
} ;

//---------------------------------------------------------------------------
MY_MUSCL_DataStructure const*
MY_MUSCL_DataStructure:: object( PDE_DomainAndFields const* dom )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: object" ) ;
   PEL_CHECK_PRE( dom != 0 ) ;

   MY_MUSCL_DataStructure const* result = 0 ;

   for( size_t i=0 ; ( result == 0 ) && i<OBJECTS->index_limit() ; ++i )
   {
      PEL_Object const* obj = OBJECTS->at(i) ;
      if( obj != 0 && obj->owner() == dom )
      {
         result = static_cast<MY_MUSCL_DataStructure const*>( obj ) ;
      }
   }
   if( result == 0 )
   {
      MY_MUSCL_DataStructure* new_item = new MY_MUSCL_DataStructure( dom ) ;
      result = new_item ;
      OBJECTS->append( new_item ) ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == dom ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
MY_MUSCL_DataStructure:: MY_MUSCL_DataStructure(
                                             PDE_DomainAndFields const* dom )
//---------------------------------------------------------------------------
   : PEL_Object( const_cast<PDE_DomainAndFields*>( dom ) )
   , DOM( dom )
   , CONNECTIVITY_SIDES( PEL_Vector::create( this, 0 ) )
   , CONNECTIVITY_BOUNDS( PEL_Vector::create( this, 0 ) )
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: MY_MUSCL_DataStructure" ) ;

   if( !structured_mesh( dom, 1.E-06 ) )
   {
      PEL_Error::object()->raise_plain(
         "*** MY_MUSCL_DataStructure error:\n"
         "    structured meshing is expected. " ) ;
   }
}

//---------------------------------------------------------------------------
MY_MUSCL_DataStructure:: ~MY_MUSCL_DataStructure( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: ~MY_MUSCL_DataStructure" ) ;
   OBJECTS->remove( this ) ;
}

//---------------------------------------------------------------------------
size_t 
MY_MUSCL_DataStructure:: upstream_node(
          PDE_CursorFEside const* fe, PDE_DiscreteField const* field ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: upstream_node" ) ;
   PEL_CHECK_PRE( fe!=0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( field != 0 ) ;

   size_t result = connectivity( fe, field )->UPS_N ;

   PEL_CHECK_POST( result < field->nb_nodes() ) ;
   PEL_CHECK_POST( field->node_is_active( result ) ) ; 
   return( result ) ;
}
 
//---------------------------------------------------------------------------
double
MY_MUSCL_DataStructure:: upstream_node_dist(
           PDE_CursorFEside const* fe, PDE_DiscreteField const* field ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: upstream_node_dist" ) ;
   PEL_CHECK_PRE( fe!=0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( field != 0 ) ;

   double result = connectivity( fe, field )->UPS_N_D ;

   PEL_CHECK_POST( result > 0. ) ; 
   return( result ) ;
}
 
//---------------------------------------------------------------------------
size_t 
MY_MUSCL_DataStructure:: upstream_node(
          PDE_LocalFEbound const* fe, PDE_DiscreteField const* field ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: upstream_node" ) ;
   PEL_CHECK_PRE( fe!=0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( field != 0 ) ;

   size_t result = connectivity( fe, field )->UPS_N ;

   PEL_CHECK_POST( result < field->nb_nodes() ) ;
   PEL_CHECK_POST( field->node_is_active( result ) ) ; 
   return( result ) ;
}
 
//---------------------------------------------------------------------------
double
MY_MUSCL_DataStructure:: upstream_node_dist(
           PDE_LocalFEbound const* fe, PDE_DiscreteField const* field ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: upstream_node_dist" ) ;
   PEL_CHECK_PRE( fe!=0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( field != 0 ) ;

   double result = connectivity( fe, field )->UPS_N_D ;

   PEL_CHECK_POST( result > 0. ) ; 
   return( result ) ;
}
 
//---------------------------------------------------------------------------
size_t 
MY_MUSCL_DataStructure:: downstream_node(
           PDE_CursorFEside const* fe, PDE_DiscreteField const* field ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: downstream_node" ) ;
   PEL_CHECK_PRE( fe!=0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( field != 0 ) ;

   size_t result = connectivity( fe, field )->DOWNS_N ;

   PEL_CHECK_POST( field->node_is_active( result ) ) ; 
   return( result ) ;
}

//---------------------------------------------------------------------------
double
MY_MUSCL_DataStructure:: downstream_node_dist(
           PDE_CursorFEside const* fe, PDE_DiscreteField const* field ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: downstream_node_dist" ) ;
   PEL_CHECK_PRE( fe!=0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( field != 0 ) ;

   double result = connectivity( fe, field )->DOWNS_N_D ;

   PEL_CHECK_POST( result < field->nb_nodes() ) ;
   PEL_CHECK_POST( result > 0. ) ; 
   return( result ) ;
}
 
//---------------------------------------------------------------------------
bool
MY_MUSCL_DataStructure:: has_up_upstream_node(
           PDE_CursorFEside const* fe, PDE_DiscreteField const* field ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: has_up_upstream_node" ) ;
   PEL_CHECK_PRE( fe!=0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( field != 0 ) ;

   return( connectivity( fe, field )->HAS_UP_UPS_N ) ;
}

//---------------------------------------------------------------------------
size_t 
MY_MUSCL_DataStructure:: up_upstream_node(
           PDE_CursorFEside const* fe, PDE_DiscreteField const* field ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: up_upstream_node" ) ;
   PEL_CHECK_PRE( fe!=0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( field != 0 ) ;
   PEL_CHECK_PRE( has_up_upstream_node( fe, field ) ) ;

   size_t result = connectivity( fe, field )->UP_UPS_N ;

   PEL_CHECK_POST( result < field->nb_nodes() ) ;
   PEL_CHECK_POST( field->node_is_active( result ) ) ; 
   return( result ) ;
}
 
//---------------------------------------------------------------------------
double
MY_MUSCL_DataStructure:: up_upstream_node_dist(
           PDE_CursorFEside const* fe, PDE_DiscreteField const* field ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: up_upstream_node_dist" ) ;
   PEL_CHECK_PRE( fe!=0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( field != 0 ) ;
   PEL_CHECK_PRE( has_up_upstream_node( fe, field ) ) ;

   double result = connectivity( fe, field )->UP_UPS_N_D ;

   PEL_CHECK_POST( result > 0. ) ; 
   return( result ) ;
}
 
//---------------------------------------------------------------------------
bool
MY_MUSCL_DataStructure:: has_up_upstream_node(
           PDE_LocalFEbound const* fe, PDE_DiscreteField const* field ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: has_up_upstream_node" ) ;
   PEL_CHECK_PRE( fe!=0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( field != 0 ) ;

   return( connectivity( fe, field )->HAS_UP_UPS_N ) ;
}

//---------------------------------------------------------------------------
size_t 
MY_MUSCL_DataStructure:: up_upstream_node(
           PDE_LocalFEbound const* fe, PDE_DiscreteField const* field ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: up_upstream_node" ) ;
   PEL_CHECK_PRE( fe!=0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( field != 0 ) ;
   PEL_CHECK_PRE( has_up_upstream_node( fe, field ) ) ;

   size_t result = connectivity( fe, field )->UP_UPS_N ;

   PEL_CHECK_POST( result < field->nb_nodes() ) ;
   PEL_CHECK_POST( field->node_is_active( result ) ) ; 
   return( result ) ;
}
 
//---------------------------------------------------------------------------
double
MY_MUSCL_DataStructure:: up_upstream_node_dist(
           PDE_LocalFEbound const* fe, PDE_DiscreteField const* field ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: up_upstream_node_dist" ) ;
   PEL_CHECK_PRE( fe!=0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( field != 0 ) ;
   PEL_CHECK_PRE( has_up_upstream_node( fe, field ) ) ;

   double result = connectivity( fe, field )->UP_UPS_N_D ;

   PEL_CHECK_POST( result > 0. ) ; 
   return( result ) ;
}
 
//---------------------------------------------------------------------------
bool
MY_MUSCL_DataStructure:: has_down_downstream_node(
           PDE_CursorFEside const* fe, PDE_DiscreteField const* field ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: has_down_downstream_node" ) ;
   PEL_CHECK_PRE( fe!=0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( field != 0 ) ;

   return( connectivity( fe, field )->HAS_DOWN_DOWNS_N ) ;
}

//---------------------------------------------------------------------------
size_t 
MY_MUSCL_DataStructure:: down_downstream_node(
          PDE_CursorFEside const* fe,  PDE_DiscreteField const* field ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: down_downstream_node" ) ;
   PEL_CHECK_PRE( fe!=0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( field != 0 ) ;
   PEL_CHECK_PRE( has_down_downstream_node( fe, field ) ) ;
 
   size_t result = connectivity( fe, field )->DOWN_DOWNS_N ;

   PEL_CHECK_POST( result < field->nb_nodes() ) ;
   PEL_CHECK_POST( field->node_is_active( result ) ) ; 
   return( result ) ;
}

//---------------------------------------------------------------------------
double
MY_MUSCL_DataStructure:: down_downstream_node_dist(
           PDE_CursorFEside const* fe, PDE_DiscreteField const* field ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: down_downstream_node_dist" ) ;
   PEL_CHECK_PRE( fe!=0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( field != 0 ) ;
   PEL_CHECK_PRE( has_down_downstream_node( fe, field ) ) ;

   double result = connectivity( fe, field )->DOWN_DOWNS_N_D ;

   PEL_CHECK_POST( result > 0. ) ; 
   return( result ) ;
}
 
//---------------------------------------------------------------------------
bool
MY_MUSCL_DataStructure:: structured_mesh( PDE_DomainAndFields const* dom,
                                          double eps ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: strutured_mesh" ) ;
   PEL_CHECK_PRE( dom!=0 ) ;

   bool result = true ;
   // Check sides
   PDE_CursorFEside* sFE = dom->create_CursorFEside(0) ;
   for( sFE->start(); result && sFE->is_valid() ; sFE->go_next() )
   {
      result = ( MY_MUSCL_DataStructure::side_normal_direction( sFE, eps ) != PEL::bad_index() ) ;
   }
   sFE->destroy() ; sFE = 0 ;

   // Check bounds 
   PDE_LocalFEbound* bFE = dom->create_LocalFEbound(0) ;
   for( bFE->start(); result && bFE->is_valid() ; bFE->go_next() )
   {
      result = ( MY_MUSCL_DataStructure::bound_normal_direction( bFE, eps ) != PEL::bad_index() ) ;
   }
   bFE->destroy() ; bFE = 0 ;
   
   return( result ) ;
}
 
//----------------------------------------------------------------------
size_t
MY_MUSCL_DataStructure:: side_normal_direction( PDE_CursorFEside const* fe,
                                                double eps )
//----------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: side_normal_direction" ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( eps > 0. ) ;

   GE_Vector const* n = fe->normal() ;
   size_t result = PEL::bad_index() ;
   for( size_t i=0 ;
        result == PEL::bad_index() && i<n->nb_components() ; ++i )
   {
      if( PEL::abs( n->component(i) )>1.-eps )
      {
         result = i ;
      }
   }
 
   PEL_CHECK_POST( result<fe->nb_space_dimensions() ) ;
   PEL_CHECK_POST( 
      IMPLIES( result != PEL::bad_index(),
               PEL::abs( fe->normal()->component(result) )>1.-eps ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
MY_MUSCL_DataStructure:: bound_normal_direction( PDE_LocalFEbound const* fe,
                                                 double eps )
//----------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: bound_normal_direction" ) ;
   PEL_CHECK_PRE( fe != 0 ) ;
   PEL_CHECK_PRE( fe->is_valid() ) ;
   PEL_CHECK_PRE( eps > 0. ) ;

   GE_Vector const* n = fe->outward_normal() ;
   size_t result = PEL::bad_index() ;
   for( size_t i=0 ;
        result == PEL::bad_index() && i<n->nb_components() ; ++i )
   {
      if( PEL::abs( n->component(i) )>1.-eps )
      {
         result = i ;
      }
   }
            
   PEL_CHECK_POST( result<fe->nb_space_dimensions() ) ;
   PEL_CHECK_POST( 
      IMPLIES( result != PEL::bad_index(),
               PEL::abs( fe->outward_normal()->component(result) )>1.-eps ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
MY_OneSideStreamCellNodes* 
MY_MUSCL_DataStructure::connectivity(
           PDE_CursorFEside const* fe, PDE_DiscreteField const* field ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure::connectivity" ) ;
   PEL_CHECK( fe != 0 ) ;
   PEL_CHECK( fe->is_valid() ) ;
   PEL_CHECK( field != 0 ) ;

   size_t id_mesh = fe->mesh_id() ;
   if( id_mesh>=CONNECTIVITY_SIDES->index_limit() )
   {
      CONNECTIVITY_SIDES->resize( id_mesh+1 ) ;
   }
   if( CONNECTIVITY_SIDES->at( id_mesh ) == 0 )
   {
      MY_OneSideStreamCellNodes*  ofscn
              = new MY_OneSideStreamCellNodes( CONNECTIVITY_SIDES, DOM, 
                                               field, fe ) ;
      CONNECTIVITY_SIDES->set_at( id_mesh, ofscn ) ;
   }
   
   MY_OneSideStreamCellNodes* result =
     static_cast<MY_OneSideStreamCellNodes*>( CONNECTIVITY_SIDES->at( id_mesh ) ) ;

   PEL_CHECK( result != 0 ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
MY_OneBoundStreamCellNodes* 
MY_MUSCL_DataStructure::connectivity(
           PDE_LocalFEbound const* fe, PDE_DiscreteField const* field ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure::connectivity" ) ;
   PEL_CHECK( fe != 0 ) ;
   PEL_CHECK( fe->is_valid() ) ;
   PEL_CHECK( field != 0 ) ;

   size_t id_mesh = fe->mesh_id() ;
   if( id_mesh>=CONNECTIVITY_BOUNDS->index_limit() )
   {
      CONNECTIVITY_BOUNDS->resize( id_mesh+1 ) ;
   }
   if( CONNECTIVITY_BOUNDS->at( id_mesh ) == 0 )
   {
      MY_OneBoundStreamCellNodes*  obscn
              = new MY_OneBoundStreamCellNodes( CONNECTIVITY_BOUNDS, DOM, 
                                               field, fe ) ;
      CONNECTIVITY_BOUNDS->set_at( id_mesh, obscn ) ;
   }
   
   MY_OneBoundStreamCellNodes* result =
     static_cast<MY_OneBoundStreamCellNodes*>( CONNECTIVITY_BOUNDS->at( id_mesh ) ) ;

   PEL_CHECK( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
MY_MUSCL_DataStructure:: display_connectivities(
                                     PDE_DiscreteField const* field,
                                     std::ostream& os,
                                     size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "MY_MUSCL_DataStructure:: display_connectivities_for_field" ) ;
   PEL_CHECK_PRE( field != 0 ) ;
   
   std::string const s( indent_width, ' ' ) ;

   PDE_CursorFEside* sFE = DOM->create_CursorFEside( 0 ) ;   
   for( sFE->start() ; sFE->is_valid() ; sFE->go_next() )
   {
      os << s << "MODULE face#" << sFE->mesh_id() << std::endl ;
      os << s << "   face_id = " << sFE->mesh_id() << std::endl ;
      os << s << "   MODULE upstream" << std::endl ;
      os << s << "      node = " << upstream_node( sFE, field ) << std::endl ;
      os << s << "      distance = " << upstream_node_dist( sFE, field ) << std::endl ;
      os << s << "   END MODULE upstream" << std::endl ;
      os << s << "   MODULE downstream" << std::endl ;
      os << s << "      node = " << downstream_node( sFE, field ) << std::endl ;
      os << s << "      distance = " << downstream_node_dist( sFE, field ) << std::endl ;
      os << s << "   END MODULE downstream" << std::endl ;
      if( has_up_upstream_node( sFE, field ) )
      {
         os << s << "   MODULE up_upstream" << std::endl ;
     	 os << s << "      node = " << up_upstream_node( sFE, field ) << std::endl ;
	 os << s << "      distance = " << up_upstream_node_dist( sFE, field ) << std::endl ;
         os << s << "   END MODULE up_upstream" << std::endl ;
      }
      if( has_down_downstream_node( sFE, field ) )
      {
         os << s << "   MODULE down_downstream" << std::endl ;
	 os << s << "      node = " << down_downstream_node( sFE, field ) << std::endl ;
	 os << s << "      distance = " << down_downstream_node_dist( sFE, field ) << std::endl ;
         os << s << "   END MODULE down_downstream" << std::endl ;
      }
      os << s << "END MODULE face#" <<sFE->mesh_id() << std::endl ;
   }
   sFE->destroy() ; sFE = 0 ;

   PDE_LocalFEbound* bFE = DOM->create_LocalFEbound( 0 ) ;   
   for( bFE->start() ; bFE->is_valid() ; bFE->go_next() )
   {
      os << s << "MODULE bound#" << bFE->mesh_id() << std::endl ;
      os << s << "   bound_id = " << bFE->mesh_id() << std::endl ;
      os << s << "   MODULE upstream" << std::endl ;
      os << s << "      node = " << upstream_node( bFE, field ) << std::endl ;
      os << s << "      distance = " << upstream_node_dist( bFE, field ) << std::endl ;
      os << s << "   END MODULE upstream" << std::endl ;
      if( has_up_upstream_node( bFE, field ) )
      {
         os << s << "   MODULE up_upstream" << std::endl ;
     	 os << s << "      node = " << up_upstream_node( bFE, field ) << std::endl ;
	 os << s << "      distance = " << up_upstream_node_dist( bFE, field ) << std::endl ;
         os << s << "   END MODULE up_upstream" << std::endl ;
      }
      os << s << "END MODULE bound#" <<bFE->mesh_id() << std::endl ;
   }
   bFE->destroy() ; bFE = 0 ;
}

//internal-------------------------------------------------------------------
MY_OneSideStreamCellNodes:: MY_OneSideStreamCellNodes(
                                     PEL_Object* a_owner,
                                     PDE_DomainAndFields const* dom,
                                     PDE_DiscreteField const* field,
                                     PDE_CursorFEside const* fe )
//internal-------------------------------------------------------------------
   : PEL_Object( a_owner )
   , UPS_N( PEL::bad_index() )
   , UPS_N_D( PEL::bad_double() )
   , DOWNS_N( PEL::bad_index() )
   , DOWNS_N_D( PEL::bad_double() )
   , HAS_UP_UPS_N( false )
   , UP_UPS_N( PEL::bad_index() )
   , UP_UPS_N_D( PEL::bad_double() )
   , HAS_DOWN_DOWNS_N( false )
   , DOWN_DOWNS_N( PEL::bad_index() )
   , DOWN_DOWNS_N_D( PEL::bad_double() )
{
   PDE_CursorFEside* sFE = dom->create_CursorFEside( 0 ) ;
   sFE->include_color( GE_Color::halo_color() ) ;
   sFE->require_field_calculation( field, PDE_LocalFE::node ) ;
   sFE->go_i_th( fe->mesh_id() ) ;

   size_t idir = MY_MUSCL_DataStructure::side_normal_direction( fe ) ;

   // Upstream and up-upstream nodes
   PDE_LocalFEcell const* c0 = sFE->adjacent_localFEcell( 0 ) ;
   PEL_CHECK( c0->nb_local_nodes( field )==1 ) ;
   UPS_N = c0->global_node( field, 0 ) ;
   UPS_N_D = 0.5*c0->polyhedron()->inter_vertices_maximum_distance( idir ) ;
   one_direction_determination( dom, sFE, c0, field, idir, UPS_N_D,
                                HAS_UP_UPS_N, UP_UPS_N, UP_UPS_N_D ) ;

   // Downstream and down_downstream nodes
   PDE_LocalFEcell const* c1 = sFE->adjacent_localFEcell( 1 ) ;
   PEL_CHECK( c1->nb_local_nodes( field )==1 ) ;
   DOWNS_N = c1->global_node( field, 0 ) ;
   DOWNS_N_D = 0.5*c1->polyhedron()->inter_vertices_maximum_distance( idir ) ;
   one_direction_determination( dom, sFE, c1, field, idir, DOWNS_N_D, 
                                HAS_DOWN_DOWNS_N, DOWN_DOWNS_N, DOWN_DOWNS_N_D ) ;

   sFE->destroy() ; sFE = 0 ;

   PEL_CHECK( field->node_is_active( UPS_N ) ) ;
   PEL_CHECK( field->node_is_active( DOWNS_N ) ) ;
   PEL_CHECK( IMPLIES( HAS_UP_UPS_N, 
		       field->node_is_active( UP_UPS_N ) ) ) ;
   PEL_CHECK( IMPLIES( HAS_DOWN_DOWNS_N, 
		       field->node_is_active( DOWN_DOWNS_N ) ) ) ;
}

//internal-------------------------------------------------------------------
MY_OneSideStreamCellNodes::~MY_OneSideStreamCellNodes( void )
//internal-------------------------------------------------------------------
{
}

//internal-------------------------------------------------------------------
void
MY_OneSideStreamCellNodes::one_direction_determination(
                        PDE_DomainAndFields const* dom,
                        PDE_CursorFEside const* fe,
                        PDE_LocalFEcell const* cFE,
                        PDE_DiscreteField const* field,
                        size_t idir,
                        double ini_dist,
                        bool& has_node,
                        size_t & node,
                        double& dist ) 
//internal-------------------------------------------------------------------
{   
   PDE_CursorFEside* sFE = dom->create_CursorFEside( 0 ) ;
   sFE->include_color( GE_Color::halo_color() ) ;
   sFE->require_field_calculation( field, PDE_LocalFE::node ) ;

   size_t_vector const& asi = cFE->adjacent_side_ids() ;
   has_node = false ;
   for( size_t is=0 ; !has_node && is<asi.size(); ++is )
   {
      sFE->go_i_th( asi( is ) ) ;
      if( sFE->mesh_id()!=fe->mesh_id() )
      {
	 has_node = ( MY_MUSCL_DataStructure::side_normal_direction( sFE )==idir ) ;
      }
   }
   
   if( has_node )
   {
      if( sFE->adjacent_localFEcell( 0 )->mesh_id()!=cFE->mesh_id() )
      {
         PDE_LocalFEcell const* c0 = sFE->adjacent_localFEcell( 0 ) ;
         PEL_ASSERT( c0->field_is_handled( field ) ) ;
         node = c0->global_node( field, 0 ) ;
         dist = ini_dist 
              + 0.5*c0->polyhedron()->inter_vertices_maximum_distance( idir ) ;
      }
      else
      {
         PDE_LocalFEcell const* c1 = sFE->adjacent_localFEcell( 1 ) ;
         node = c1->global_node( field, 0 ) ;
         dist = ini_dist 
              + 0.5*c1->polyhedron()->inter_vertices_maximum_distance( idir ) ;
      }
   }
   else // Check
   {
      PDE_LocalFEbound* bFE = dom->create_LocalFEbound( 0 ) ;
      bFE->include_color( GE_Color::halo_color() ) ;
      bFE->require_field_calculation( field, PDE_LocalFE::node ) ;
      
      bool found = false ;
      size_t_vector const& abi = cFE->adjacent_bound_ids() ;
      for( size_t ib=0; !found && ib<abi.size(); ++ib )
      {
	 bFE->go_i_th( abi( ib ) ) ;
	 found = ( MY_MUSCL_DataStructure::bound_normal_direction( bFE )== idir ) ;
      }
      if( !found )
      {
	 PEL_Error::object()->raise_internal(
            "*** MY_MUSCL_DataStructure error:\n"
            "    Error with side cell neighbours determination." ) ;
      }

      bFE->destroy() ; bFE = 0 ;
   }

   sFE->destroy() ; sFE = 0 ;
}

//internal-------------------------------------------------------------------
MY_OneBoundStreamCellNodes:: MY_OneBoundStreamCellNodes(
                                     PEL_Object* a_owner,
                                     PDE_DomainAndFields const* dom,
                                     PDE_DiscreteField const* field,
                                     PDE_LocalFEbound const* fe )
//internal-------------------------------------------------------------------
   : PEL_Object( a_owner )
   , UPS_N( PEL::bad_index() )
   , UPS_N_D( PEL::bad_double() )
   , HAS_UP_UPS_N( false )
   , UP_UPS_N( PEL::bad_index() )
   , UP_UPS_N_D( PEL::bad_double() )
{
   PDE_LocalFEbound* bFE = dom->create_LocalFEbound( 0 ) ;
   bFE->include_color( GE_Color::halo_color() ) ;
   bFE->require_field_calculation( field, PDE_LocalFE::node ) ;
   bFE->go_i_th( fe->mesh_id() ) ;

   size_t idir = MY_MUSCL_DataStructure::bound_normal_direction( fe ) ;

   // Upstream and up-upstream nodes
   PDE_LocalFEcell* c = dom->create_LocalFEcell( 0 ) ;
   c->include_color( GE_Color::halo_color() ) ;
   c->require_field_calculation( field, PDE_LocalFE::node ) ;
   c->go_i_th( fe->adjacent_cell_id() ) ;
   PEL_CHECK( c->nb_local_nodes( field )==1 ) ;
   UPS_N = c->global_node( field, 0 ) ;
   UPS_N_D = 0.5*c->polyhedron()->inter_vertices_maximum_distance( idir ) ;
   one_direction_determination( dom, bFE, c, field, idir, UPS_N_D,
                                HAS_UP_UPS_N, UP_UPS_N, UP_UPS_N_D ) ;

   c->destroy() ; c = 0 ;
   bFE->destroy() ; bFE = 0 ;

   PEL_CHECK( field->node_is_active( UPS_N ) ) ;
   PEL_CHECK( IMPLIES( HAS_UP_UPS_N, 
		       field->node_is_active( UP_UPS_N ) ) ) ;
}

//internal-------------------------------------------------------------------
MY_OneBoundStreamCellNodes::~MY_OneBoundStreamCellNodes( void )
//internal-------------------------------------------------------------------
{
}

//internal-------------------------------------------------------------------
void
MY_OneBoundStreamCellNodes::one_direction_determination(
                        PDE_DomainAndFields const* dom,
                        PDE_LocalFEbound const* fe,
                        PDE_LocalFEcell const* cFE,
                        PDE_DiscreteField const* field,
                        size_t idir,
                        double ini_dist,
                        bool& has_node,
                        size_t & node,
                        double& dist ) 
//internal-------------------------------------------------------------------
{   
   PDE_CursorFEside* sFE = dom->create_CursorFEside( 0 ) ;
   sFE->include_color( GE_Color::halo_color() ) ;
   sFE->require_field_calculation( field, PDE_LocalFE::node ) ;
   size_t_vector const& asi = cFE->adjacent_side_ids() ;
   has_node = false ;
   for( size_t is=0 ; !has_node && is<asi.size(); ++is )
   {
      sFE->go_i_th( asi( is ) ) ;
      has_node = ( MY_MUSCL_DataStructure::side_normal_direction( sFE )==idir ) ;
   }
   
   if( has_node )
   {
      if( sFE->adjacent_localFEcell( 0 )->mesh_id()!=cFE->mesh_id() )
      {
         PDE_LocalFEcell const* c0 = sFE->adjacent_localFEcell( 0 ) ;
         PEL_ASSERT( c0->field_is_handled( field ) ) ;
         node = c0->global_node( field, 0 ) ;
         dist = ini_dist 
              + 0.5*c0->polyhedron()->inter_vertices_maximum_distance( idir ) ;
      }
      else
      {
         PDE_LocalFEcell const* c1 = sFE->adjacent_localFEcell( 1 ) ;
         node = c1->global_node( field, 0 ) ;
         dist = ini_dist 
              + 0.5*c1->polyhedron()->inter_vertices_maximum_distance( idir ) ;
      }
   }
   else // Check
   {
      PDE_LocalFEbound* bFE = dom->create_LocalFEbound( 0 ) ;
      bFE->include_color( GE_Color::halo_color() ) ;
      bFE->require_field_calculation( field, PDE_LocalFE::node ) ;
      
      bool found = false ;
      size_t_vector const& abi = cFE->adjacent_bound_ids() ;
      for( size_t ib=0; !found && ib<abi.size(); ++ib )
      {
	 bFE->go_i_th( abi( ib ) ) ;
	 found = ( MY_MUSCL_DataStructure::bound_normal_direction( bFE )== idir ) ;
      }
      if( !found )
      {
	 PEL_Error::object()->raise_internal(
            "*** MY_MUSCL_DataStructure error:\n"
            "    Error with side cell neighbours determination." ) ;
      }

      bFE->destroy() ; bFE = 0 ;
   }

   sFE->destroy() ; sFE = 0 ;
}
