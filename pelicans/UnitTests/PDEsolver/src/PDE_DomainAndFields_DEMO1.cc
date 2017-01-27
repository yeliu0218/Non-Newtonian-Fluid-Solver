#include <PDE_DomainAndFields_DEMO1.hh>

#include <GE_Color.hh>

#include <PEL_assertions.hh>
#include <PEL_ModuleExplorer.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_SetOfBCs.hh>
#include <PDE_SetOfDiscreteFields.hh>

#include <fstream>
#include <iostream>

using std::endl ;
using std::ofstream ;
using std::string ;

PDE_DomainAndFields_DEMO1 const* 
PDE_DomainAndFields_DEMO1::PROTOTYPE = new PDE_DomainAndFields_DEMO1() ;

//----------------------------------------------------------------------------
PDE_DomainAndFields_DEMO1:: PDE_DomainAndFields_DEMO1( void )
//----------------------------------------------------------------------------
   : PEL_Application( "PDE_DomainAndFields_DEMO1" )
{
}

//----------------------------------------------------------------------------
PDE_DomainAndFields_DEMO1*
PDE_DomainAndFields_DEMO1:: create_replica( PEL_Object* a_owner, 
                                 PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields_DEMO1:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   PDE_DomainAndFields_DEMO1* result = new PDE_DomainAndFields_DEMO1( a_owner,
                                                                      exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PDE_DomainAndFields_DEMO1:: PDE_DomainAndFields_DEMO1( PEL_Object* a_owner, 
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
PDE_DomainAndFields_DEMO1:: ~PDE_DomainAndFields_DEMO1( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
PDE_DomainAndFields_DEMO1:: run( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DomainAndFields_DEMO1:: run" ) ;
   
   ofstream ofs( OUT_NAME.c_str() ) ;

   PDE_SetOfDiscreteFields const* sdf = DOM->set_of_discrete_fields() ;

   ofs << "number of fields: " << sdf->nb_fields() << endl ;

   if( sdf->has( "my_nice_field" ) )
   {
      PDE_DiscreteField const* ff = sdf->item( "my_nice_field" ) ;
      ofs << ff->name() << endl ; // "my_nice_field"
   }

   for( sdf->start() ; sdf->is_valid() ; sdf->go_next() )
   {
      PDE_DiscreteField const* ff = sdf->item() ;
      ofs << ff->name() << endl ;
   }

   PDE_SetOfBCs const* bcs = DOM->set_of_boundary_conditions() ;
   GE_Color const* col = GE_Color::object( "turquoise" ) ;
   PDE_DiscreteField const* ff = sdf->item( "my_nice_field" ) ;
   if( bcs->has_BC( col, ff ) )
   {
      PEL_ModuleExplorer const* ee = bcs->BC_explorer( col, ff ) ;
      string const& bc_type = ee->string_data( "type" ) ;
      if( bc_type == "my_user_bc" )
      {
         double xx = ee->double_data( "the_coef" ) ;
         ofs << "coef of bc: " << xx << endl ;
         string const& ss = ee->string_data( "the_string" ) ;
         ofs << "string of bc: " << ss << endl ;
      }
   }

   ofs.close() ;
}

