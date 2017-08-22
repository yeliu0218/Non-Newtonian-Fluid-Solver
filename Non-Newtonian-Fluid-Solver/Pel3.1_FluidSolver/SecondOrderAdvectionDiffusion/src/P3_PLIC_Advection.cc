//This is the implementation of the Piecewise Linear Interface Reconstruction (PLIC) scheme with a correction scheme to resolve the moving contact line
//Notice that it only works with a structural rectangular domain, with rectangular elements, while the mesh size is allowed to be non-uniform


#include <P3_PLIC_Advection.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_ContextSimple.hh>
#include <PEL_DataWithContext.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Double.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_MemoryTracer.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>

#include <GE_Color.hh>
#include <GE_QRprovider.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>
#include <GE_Vector.hh>

#include <PDE.hh>
#include <PDE_CursorFEside.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEbound.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_ResultSaver.hh>
#include <PDE_SetOfBCs.hh>
#include <PDE_SystemNumbering.hh>

#include <FE.hh>
#include <FE_Parameter.hh>
#include <FE_SetOfParameters.hh>
#include <FE_TimeIterator.hh>
#include <FE_LocalTimeIteratorAdapter.hh>

// KWB check if we need all
#include <LA_Matrix.hh>	
#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>
#include <LA_Solver.hh>

#include <LA_Vector.hh> // KWB do we need it?

#include <P3_PLIC_Advection.hh>

#include <cmath>
#include <iostream>


using std::string ;
using std::cout ;
using std::endl ;

P3_PLIC_Advection const*
P3_PLIC_Advection::PROTOTYPE = new P3_PLIC_Advection() ;

//---------------------------------------------------------------------------
P3_PLIC_Advection:: P3_PLIC_Advection( void )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( "P3_PLIC_Advection" )
   , CC( 0 )
   , UU( 0 )
   , cFE( 0 )
   , sFE( 0 )
   , bFE( 0 )
   , EPSILON( 0. )
   , Cmax( 0. )
   , useCOURANT( false )
   , FIELDS( NULL )
   , LEVEL0( PEL::bad_index() )
   , FIELDS_TABLE( 0 )
   , COURANT( 0 )
   , A_OWNER( NULL )
   , DOM( NULL )
   , PRMS ( NULL )
   , EXP ( NULL )
   , NX ( 0 )
   , NY ( 0 )
   , Xvector ( 0 )
   , Yvector ( 0 )
   , Xcoord( 0 )
{
}

//---------------------------------------------------------------------------
P3_PLIC_Advection*
P3_PLIC_Advection:: create_replica( PEL_Object* a_owner,
		           PDE_DomainAndFields const* dom,
		           FE_SetOfParameters const* prms,
		           PEL_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "P3_PLIC_Advection:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   P3_PLIC_Advection* result = new P3_PLIC_Advection( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
P3_PLIC_Advection:: P3_PLIC_Advection( PEL_Object* a_owner,
		     PDE_DomainAndFields const* dom,
		     FE_SetOfParameters const* prms,
		     PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , CC( dom->set_of_discrete_fields()->item( exp->string_data("AD_unknown_field") ) )
   , UU( dom->set_of_discrete_fields()->item( exp->string_data("AD_param_advective_velocity") ))
   , L_UPDATE(0)    
   , cFE( dom->create_LocalFEcell( this ) )
   , sFE( dom->create_CursorFEside( this ) )
   , bFE( dom->create_LocalFEbound( this ) )
   , EPSILON( 1.E-03 )
   , useCOURANT( false )
   , FIELDS( dom->set_of_discrete_fields() )
   , LEVEL0( PEL::bad_index() )
   , FIELDS_TABLE( 0 )
   , COURANT( 0 )
   , A_OWNER( a_owner )
   , DOM( dom )
   , PRMS ( prms )
   , EXP ( exp )  
   , Xvector( doubleVector( 0 ) )
   , Yvector( doubleVector( 0 ) )
   , Xcoord( doubleVector( 0 ) )
{

   PEL_LABEL( "P3_PLIC_Advection:: P3_PLIC_Advection" ) ;

   cFE->require_field_calculation( CC, PDE_LocalFE::N ) ;
   sFE->require_field_calculation( CC, PDE_LocalFE::node ) ;
   cFE->require_field_calculation( CC, PDE_LocalFE::dN ) ;
   bFE->require_field_calculation( CC, PDE_LocalFE::node ) ;
   cFE->require_field_calculation( UU, PDE_LocalFE::N ) ;
   sFE->require_field_calculation( UU, PDE_LocalFE::node ) ;

   //GRAD_C->transfer_cell_calculation_requirements( cFE, FE_Parameter::Val ) ;

   if (exp->has_module( "CourantNumber" ) )
   {
                PEL_ModuleExplorer const* se = 0 ;
		se = exp->create_subexplorer( 0, "CourantNumber" ) ;
		useCOURANT = true;
		COURANT = se->double_data( "Courant" );
		LEVEL0 = se->int_data( "active_level" ) ;  
                doubleVector const& cleft = se->doubleVector_data( "vertices_coordinate_0" ) ;
                doubleVector const& cbottom = se->doubleVector_data( "vertices_coordinate_1" ) ; 
                NX = cleft.size()-1;
                NY = cbottom.size()-1;  
                Xvector.resize( NX );
                Yvector.resize( NY );
                Xcoord.resize( NX );
                for( size_t i = 0 ; i<NX ; i++ )
                {   
                   Xvector( i ) = cleft(i+1)-cleft(i) ;
                   Xcoord( i ) = 0.5*(cleft(i+1)+cleft(i));
                }
                for( size_t i = 0 ; i<NY ; i++ )
                   Yvector( i ) = cbottom(i+1)-cbottom(i) ;    
                Cmax = se->double_data("Cmax");
		FIELDS_TABLE = se->stringVector_data( "discrete_fields" ) ;
		t_global = FE_TimeIterator::create( A_OWNER, EXP );
		se->destroy() ; se = 0 ;
   }
   
}

//---------------------------------------------------------------------------
P3_PLIC_Advection:: ~P3_PLIC_Advection( void )
//---------------------------------------------------------------------------
{
}
//-----------------------------------------------------------------------------
void
P3_PLIC_Advection:: do_before_inner_iterations_stage( FE_TimeIterator const* t_it )
//-----------------------------------------------------------------------------
{
   PEL_LABEL( "P3_PLIC_Advection:: do_before_inner_iterations_stage" ) ;
   PEL_CHECK_PRE( do_before_inner_iterations_stage_PRE( t_it ) ) ;

   start_total_timer( "P3_PLIC_Advection:: do_before_inner_iterations_stage" ) ;
   stop_total_timer() ;

   PEL_CHECK_POST( do_before_inner_iterations_stage_POST( t_it ) ) ;
}
//---------------------------------------------------------------------------
void P3_PLIC_Advection:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "P3_PLIC_Advection:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "P3_PLIC_Advection:: do_one_inner_iteration" ) ;
   // --------------
    
   PEL_MemoryTracer::object()->start_event(
      "P3_PLIC_Advection::do_one_inner_iteration \""+CC->name()+"\"" ) ;
  
   reset_discrete_problem( t_it );
   
   // Assembling
   start_assembling_timer() ; 

   PEL_MemoryTracer::object()->start_event( "assemble" ) ;
   loop_on_cellsXX( t_it ) ;
   loop_on_cellsZZ( t_it ) ;
   loop_on_cellsZZ( t_it ) ;
   loop_on_cellsXX( t_it ) ;
   PEL_MemoryTracer::object()->stop_event() ;
   
   stop_assembling_timer() ;

   stop_total_timer() ;
   PEL_MemoryTracer::object()->stop_event() ;
}

//------------------------------------------------------------------------
void
P3_PLIC_Advection:: print_mesh( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "P3_PLIC_Advection:: print_mesh" ) ;
   for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
       size_t n_U = cFE->global_node( UU, 0 ) ;
       PEL::out() << indent() << " global_node UU " << n_U << std::endl; 
       cFE->print_current_mesh( PEL::out(), 3 ) ;
   }
}

//Move the interface along the X-axis
//------------------------------------------------------------------------
void
P3_PLIC_Advection:: loop_on_cellsXX( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "P3_PLIC_Advection:: loop_on_cellsXX" ) ;
   double const dt = 0.5*t_it->time_step() ;
   double s1,s2,s3,s4,mm1,mm2,V1,V2,mx,mz,alfa;
   size_t n_C, n_U, inv;
   double temp;
   for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
       n_C = cFE->global_node( CC, 0 ) ;
       temp = CC->DOF_value(0,n_C);
       CC->set_DOF_value(4, n_C, temp );
       if (temp > 1.0)
       CC->set_DOF_value(0, n_C, 2.0-temp);
   }
   for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
       //Compute the velocity
       n_U = cFE->global_node( UU, 0 ) ;
       if (n_U>2*NY)
       {
          s1= ( UU->DOF_value( 0, n_U, 0 )+ UU->DOF_value( 0, n_U+1, 0 ) )*0.5 ;
          s2= ( UU->DOF_value( 0, n_U+NY+1, 0 ) + UU->DOF_value( 0, n_U+NY+2, 0 ) )*0.5 ;
       }
       else if (n_U%2==0 && n_U>0)
       {
          s1= ( UU->DOF_value( 0, n_U, 0 )+ UU->DOF_value( 0, n_U+2, 0 ) )*0.5 ;
          s2= ( UU->DOF_value( 0, n_U/2+2*NY+2, 0 )+ UU->DOF_value( 0, n_U/2+2*NY+3, 0 ) )*0.5 ;
       }
       else if (n_U%2==1 && n_U>1)
       {
          s1= ( UU->DOF_value( 0, n_U, 0 )+ UU->DOF_value( 0, n_U+2, 0 ) )*0.5 ;
          s2= ( UU->DOF_value( 0, n_U-1, 0 )+ UU->DOF_value( 0, n_U+1, 0 ) )*0.5 ;
       }
       else if (n_U==1)
       {
          s1= ( UU->DOF_value( 0, 1, 0 )+ UU->DOF_value( 0, 2, 0 ) )*0.5 ;
          s2= ( UU->DOF_value( 0, 2*NY+2, 0 )+ UU->DOF_value( 0, 2*NY+3, 0 ) )*0.5 ;
       }
       else
       {
          s1= ( UU->DOF_value( 0, 0, 0 )+ UU->DOF_value( 0, 3, 0 ) )*0.5 ;
          s2= ( UU->DOF_value( 0, 1, 0 )+ UU->DOF_value( 0, 2, 0 ) )*0.5 ;
       }
       
      //Compute the flux
       n_C = cFE->global_node( CC, 0 ) ;
       s1=s1*dt/Xvector(n_C/NY);
       s2=s2*dt/Xvector(n_C/NY);
       if ( CC->DOF_value(0,n_C)==0.0 )
       {
          CC->set_DOF_value(1, n_C, 0.0);
          CC->set_DOF_value(2, n_C, 0.0); 
          CC->set_DOF_value(3, n_C, 0.0);  
       }
       else if ( CC->DOF_value(0,n_C)==1.0 )
       {
          s3=s1*Xvector(n_C/NY)/(n_C/NY==0 ? Xvector(n_C/NY) : Xvector(n_C/NY-1));
          s4=s2*Xvector(n_C/NY)/(n_C/NY==NX-1 ? Xvector(n_C/NY) : Xvector(n_C/NY+1));
          CC->set_DOF_value(1, n_C, (-s3 > 0.0) ? -s3 : 0.0 );
          CC->set_DOF_value(2, n_C, 1.0 - (s1 > 0.0 ? s1 : 0.0) + ( s2 < 0.0 ? s2 : 0.0 )); 
          CC->set_DOF_value(3, n_C, s4 > 0.0 ? s4 : 0.0);
       }
       else
       {
          if (n_C<NY-1 && n_C>0)
          {
             mm1=CC->DOF_value(0,n_C-1)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C+1);
             mm2=CC->DOF_value(0,n_C+NY-1)+2.0*CC->DOF_value(0,n_C+NY)+CC->DOF_value(0,n_C+NY+1);
             mx=(mm1-mm2)*2.0;
             mm1=CC->DOF_value(0,n_C-1)+2.0*CC->DOF_value(0,n_C-1)+CC->DOF_value(0,n_C+NY-1);
             mm2=CC->DOF_value(0,n_C+1)+2.0*CC->DOF_value(0,n_C+1)+CC->DOF_value(0,n_C+NY+1);
             mz=mm1-mm2;
          }
          else if (n_C>(NX-1)*NY && n_C<NX*NY-1)
          {
             mm1=CC->DOF_value(0,n_C-NY-1)+2.0*CC->DOF_value(0,n_C-NY)+CC->DOF_value(0,n_C-NY+1);
             mm2=CC->DOF_value(0,n_C-1)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C+1);
             mx=(mm1-mm2)*2.0;
             mm1=CC->DOF_value(0,n_C-NY-1)+2.0*CC->DOF_value(0,n_C-1)+CC->DOF_value(0,n_C-1);
             mm2=CC->DOF_value(0,n_C-NY+1)+2.0*CC->DOF_value(0,n_C+1)+CC->DOF_value(0,n_C+1);
             mz=mm1-mm2;
          }
          else if (n_C%NY==0 && n_C>0 && n_C<(NX-1)*NY)
          {
             mm1=CC->DOF_value(0,n_C-NY)+2.0*CC->DOF_value(0,n_C-NY)+CC->DOF_value(0,n_C-NY+1);
             mm2=CC->DOF_value(0,n_C+NY)+2.0*CC->DOF_value(0,n_C+NY)+CC->DOF_value(0,n_C+NY+1);
             mx=mm1-mm2;
             mm1=CC->DOF_value(0,n_C-NY)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C+NY);
             mm2=CC->DOF_value(0,n_C-NY+1)+2.0*CC->DOF_value(0,n_C+1)+CC->DOF_value(0,n_C+NY+1);
             mz=(mm1-mm2)*2.0;
          }
          else if (n_C%NY==NY-1 && n_C>NY-1 && n_C<NX*NY-1)
          {
             mm1=CC->DOF_value(0,n_C-NY-1)+2.0*CC->DOF_value(0,n_C-NY)+CC->DOF_value(0,n_C-NY);
             mm2=CC->DOF_value(0,n_C+NY-1)+2.0*CC->DOF_value(0,n_C+NY)+CC->DOF_value(0,n_C+NY);
             mx=mm1-mm2;
             mm1=CC->DOF_value(0,n_C-NY-1)+2.0*CC->DOF_value(0,n_C-1)+CC->DOF_value(0,n_C+NY-1);
             mm2=CC->DOF_value(0,n_C-NY)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C+NY);
             mz=(mm1-mm2)*2.0;
          }
          else if (n_C==NY-1)
          {
             mm1=CC->DOF_value(0,n_C-1)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C);
             mm2=CC->DOF_value(0,n_C+NY-1)+2.0*CC->DOF_value(0,n_C+NY)+CC->DOF_value(0,n_C+NY);
             mx=mm1-mm2;
             mm1=CC->DOF_value(0,n_C-1)+2.0*CC->DOF_value(0,n_C-1)+CC->DOF_value(0,n_C+NY-1);
             mm2=CC->DOF_value(0,n_C)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C+NY);
             mz=mm1-mm2;
          }
          else if (n_C==0)
          {
             mm1=CC->DOF_value(0,n_C)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C+1);
             mm2=CC->DOF_value(0,n_C+NY)+2.0*CC->DOF_value(0,n_C+NY)+CC->DOF_value(0,n_C+NY+1);
             mx=mm1-mm2;
             mm1=CC->DOF_value(0,n_C)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C+NY);
             mm2=CC->DOF_value(0,n_C+1)+2.0*CC->DOF_value(0,n_C+1)+CC->DOF_value(0,n_C+NY+1);
             mz=mm1-mm2;
          }
          else if (n_C==(NX-1)*NY)
          {
             mm1=CC->DOF_value(0,n_C-NY)+2.0*CC->DOF_value(0,n_C-NY)+CC->DOF_value(0,n_C-NY+1);
             mm2=CC->DOF_value(0,n_C)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C+1);
             mx=mm1-mm2;
             mm1=CC->DOF_value(0,n_C-NY)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C);
             mm2=CC->DOF_value(0,n_C-NY+1)+2.0*CC->DOF_value(0,n_C+1)+CC->DOF_value(0,n_C+1);
             mz=mm1-mm2;
          }
          else if (n_C==NX*NY-1)
          {
             mm1=CC->DOF_value(0,n_C-NY-1)+2.0*CC->DOF_value(0,n_C-NY)+CC->DOF_value(0,n_C-NY);
             mm2=CC->DOF_value(0,n_C-1)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C);
             mx=mm1-mm2;
             mm1=CC->DOF_value(0,n_C-NY-1)+2.0*CC->DOF_value(0,n_C-1)+CC->DOF_value(0,n_C+NY-1);
             mm2=CC->DOF_value(0,n_C-NY)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C);
             mz=mm1-mm2;
          }
          else
          {
             mm1=CC->DOF_value(0,n_C-NY-1)+2.0*CC->DOF_value(0,n_C-NY)+CC->DOF_value(0,n_C-NY+1);
             mm2=CC->DOF_value(0,n_C+NY-1)+2.0*CC->DOF_value(0,n_C+NY)+CC->DOF_value(0,n_C+NY+1);
             mx=mm1-mm2;
             mm1=CC->DOF_value(0,n_C-NY-1)+2.0*CC->DOF_value(0,n_C-1)+CC->DOF_value(0,n_C+NY-1);
             mm2=CC->DOF_value(0,n_C-NY+1)+2.0*CC->DOF_value(0,n_C+1)+CC->DOF_value(0,n_C+NY+1);
             mz=mm1-mm2;
          }          
   
          inv = 0;
          if (mx<0.0)
          {
             mm1=-s1;
             s1=-s2;
             s2=mm1;
             mx=-mx;
             inv = 1;
          }
          mx=mx+1.E-50;
          mz=PEL::abs(mz)+1.E-50;
          mm2= mx > mz ? mx : mz;
          mx=mx/mm2;
          mz=mz/mm2;

          mm1= mx < mz ? mx : mz;
          V1=0.5*mm1;
          V2=1.0-V1;
          if (CC->DOF_value(0,n_C)<=V1)
             alfa=sqrt(2.0*CC->DOF_value(0,n_C)*mm1);
          else if (CC->DOF_value(0,n_C)<=V2)
             alfa=CC->DOF_value(0,n_C)+0.5*mm1;
          else 
             alfa=mm1+1.0-sqrt(2.0*(1.0-CC->DOF_value(0,n_C))*mm1);

          mx=mx/(1.0-s1+s2);
          alfa=alfa+mx*s1;
          
          mm1=(-s1 > 0.0) ? -s1 : 0.0;
          mm2=alfa+mx*mm1;
          CC->set_DOF_value(1, n_C, vol2(mx,mz,mm2,mm1));

          mm1=1.0-(s1>0.0?s1:0.0)+(s2<0.0?s2:0.0);
          mm2=alfa-mx*(s1>0.0?s1:0.0);
          CC->set_DOF_value(2, n_C, vol2(mx,mz,mm2,mm1) );
          
          mm1=s2>0.0?s2:0.0;
          mm2=alfa-mx;
          CC->set_DOF_value(3, n_C, vol2(mx,mz,mm2,mm1) );
          
          // it is very likely that this is the problem!!
          if (inv==1)
          {
             mm1=CC->DOF_value(1,n_C);
             CC->set_DOF_value(1, n_C, CC->DOF_value(3,n_C) );
             CC->set_DOF_value(3, n_C, mm1 );
          }

          CC->set_DOF_value(1,n_C,CC->DOF_value(1,n_C)*Xvector(n_C/NY)/(n_C/NY==0 ? Xvector(n_C/NY) : Xvector(n_C/NY-1)));
          CC->set_DOF_value(3,n_C,CC->DOF_value(3,n_C)*Xvector(n_C/NY)/(n_C/NY==NX-1 ? Xvector(n_C/NY) : Xvector(n_C/NY+1)));
      }
   }

   //update C
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      n_C = cFE->global_node( CC, 0 ) ;
      double l,r;
      r = (n_C+NY<NX*NY-1) ? CC->DOF_value(1,n_C+NY) : 0.0;
      l = (n_C>NY) ? CC->DOF_value(3,n_C-NY) : 0.0;
      CC->set_DOF_value(0, n_C, l+CC->DOF_value(2,n_C)+r );
   }
 
   //regularize C 
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      n_C = cFE->global_node( CC, 0 ) ;
      mm1=0.0>CC->DOF_value(0,n_C)?0.0:CC->DOF_value(0,n_C);
      mm2=1.0<mm1?1.0:mm1;
      CC->set_DOF_value(0, n_C, mm2);
   }
   for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
       n_C = cFE->global_node( CC, 0 ) ;
       temp = CC->DOF_value(0,n_C);
       if (CC->DOF_value(4, n_C)>1.0 && temp < 1.0)   //previous C2 or interface; current not C1 so it has to be C2 or interface
          CC->set_DOF_value(0, n_C, 2.0-temp);        
       //if (CC->DOF_value(4, n_C)>1.0 && temp == 1.0)   //previous C2; current C1: do nothing!
          //CC->set_DOF_value(0, n_C, 2.0-temp);
       //if (CC->DOF_value(4, n_C)==1.0 && temp == 1.0)  //previous C1; current C1: do nothing!
          //CC->set_DOF_value(0, n_C, 2.0-temp);
       if (CC->DOF_value(4, n_C)==1.0 && temp < 1.0)   //previous C1; current interface(C0 or C2)
       { 
          //check the left cell(closer to symmetry axis) C=0: IC0; 0<C<1: IC0; C=1: IC2; 1<C<2: IC2; C=2:unlikely
          if (CC->DOF_value(4, n_C-NY)>=1.0)
             CC->set_DOF_value(0, n_C, 2.0-temp); 
       }
       if (n_C%NY==0)
          CC->set_DOF_value(0, n_C, CC->DOF_value(4, n_C)); 
   }
}


//Move the interface along the Z-axis
//------------------------------------------------------------------------
void
P3_PLIC_Advection:: loop_on_cellsZZ( FE_TimeIterator const* t_it )
//------------------------------------------------------------------------
{
   PEL_LABEL( "P3_PLIC_Advection:: loop_on_cellsZZ" ) ;
   double const dt = 0.5*t_it->time_step() ;
   double s1,s2,s3,s4,mm1,mm2,V1,V2,mx,mz,alfa;
   size_t n_C, n_U, inv;
   double temp;
   for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
       n_C = cFE->global_node( CC, 0 ) ;
       temp = CC->DOF_value(0,n_C);
       CC->set_DOF_value(4, n_C, temp );
       if (temp > 1.0)
       CC->set_DOF_value(0, n_C, 2.0-temp);
   }
   for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
       n_U = cFE->global_node( UU, 0 ) ;
       if (n_U>2*NY)
       {
          s1= ( UU->DOF_value( 0, n_U, 1 )+ UU->DOF_value( 0, n_U+NY+1, 1 ) )*0.5 ;
          s2= ( UU->DOF_value( 0, n_U+1, 1 ) + UU->DOF_value( 0, n_U+NY+2, 1 ) )*0.5 ;
       }
       else if (n_U%2==0 && n_U>0)
       {
          s1= ( UU->DOF_value( 0, n_U, 1 )+ UU->DOF_value( 0, n_U/2+2*NY+2, 1 ) )*0.5 ;
          s2= ( UU->DOF_value( 0, n_U+2, 1 )+ UU->DOF_value( 0, n_U/2+2*NY+3, 1 ) )*0.5 ;
       }
       else if (n_U%2==1 && n_U>1)
       {
          s1= ( UU->DOF_value( 0, n_U, 1 )+ UU->DOF_value( 0, n_U-1, 1 ) )*0.5 ;
          s2= ( UU->DOF_value( 0, n_U+2, 1 )+ UU->DOF_value( 0, n_U+1, 1 ) )*0.5 ;
       }
       else if (n_U==1)
       {
          s1= ( UU->DOF_value( 0, 1, 1 )+ UU->DOF_value( 0, 2*NY+2, 1 ) )*0.5 ;
          s2= ( UU->DOF_value( 0, 2, 1 )+ UU->DOF_value( 0, 2*NY+3, 1 ) )*0.5 ;
       }    
       else
       {
          s1= ( UU->DOF_value( 0, 0, 1 )+ UU->DOF_value( 0, 1, 1 ) )*0.5 ;
          s2= ( UU->DOF_value( 0, 3, 1 )+ UU->DOF_value( 0, 2, 1 ) )*0.5 ;
       }
       n_C = cFE->global_node( CC, 0 ) ;
       s1=s1*dt/Yvector(n_C%NY);
       s2=s2*dt/Yvector(n_C%NY);
       if ( CC->DOF_value(0,n_C)==0.0 )
       {
          CC->set_DOF_value(1, n_C, 0.0);
          CC->set_DOF_value(2, n_C, 0.0); 
          CC->set_DOF_value(3, n_C, 0.0);  
       }
       else if ( CC->DOF_value(0,n_C)==1.0 )
       {
          s3=s1*Yvector(n_C%NY)/(n_C%NY==0 ? Yvector(n_C%NY) : Yvector(n_C%NY-1));
          s4=s2*Yvector(n_C%NY)/(n_C%NY==NY-1 ? Yvector(n_C%NY) : Yvector(n_C%NY+1));
          CC->set_DOF_value(1, n_C, (-s3 > 0.0) ? -s3 : 0.0 );
          CC->set_DOF_value(2, n_C, 1.0 - (s1 > 0.0 ? s1 : 0.0) + ( s2 < 0.0 ? s2 : 0.0 )); 
          CC->set_DOF_value(3, n_C, s4 > 0.0 ? s4 : 0.0);
       }
       else
       {
          if (n_C<NY-1 && n_C>0)
          {
             mm1=CC->DOF_value(0,n_C-1)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C+1);
             mm2=CC->DOF_value(0,n_C+NY-1)+2.0*CC->DOF_value(0,n_C+NY)+CC->DOF_value(0,n_C+NY+1);
             mx=(mm1-mm2)*2.0;
             mm1=CC->DOF_value(0,n_C-1)+2.0*CC->DOF_value(0,n_C-1)+CC->DOF_value(0,n_C+NY-1);
             mm2=CC->DOF_value(0,n_C+1)+2.0*CC->DOF_value(0,n_C+1)+CC->DOF_value(0,n_C+NY+1);
             mz=mm1-mm2;
          }
          else if (n_C>(NX-1)*NY && n_C<NX*NY-1)
          {
             mm1=CC->DOF_value(0,n_C-NY-1)+2.0*CC->DOF_value(0,n_C-NY)+CC->DOF_value(0,n_C-NY+1);
             mm2=CC->DOF_value(0,n_C-1)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C+1);
             mx=(mm1-mm2)*2.0;
             mm1=CC->DOF_value(0,n_C-NY-1)+2.0*CC->DOF_value(0,n_C-1)+CC->DOF_value(0,n_C-1);
             mm2=CC->DOF_value(0,n_C-NY+1)+2.0*CC->DOF_value(0,n_C+1)+CC->DOF_value(0,n_C+1);
             mz=mm1-mm2;
          }
          else if (n_C%NY==0 && n_C>0 && n_C<(NX-1)*NY)
          {
             mm1=CC->DOF_value(0,n_C-NY)+2.0*CC->DOF_value(0,n_C-NY)+CC->DOF_value(0,n_C-NY+1);
             mm2=CC->DOF_value(0,n_C+NY)+2.0*CC->DOF_value(0,n_C+NY)+CC->DOF_value(0,n_C+NY+1);
             mx=mm1-mm2;
             mm1=CC->DOF_value(0,n_C-NY)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C+NY);
             mm2=CC->DOF_value(0,n_C-NY+1)+2.0*CC->DOF_value(0,n_C+1)+CC->DOF_value(0,n_C+NY+1);
             mz=(mm1-mm2)*2.0;
          }
          else if (n_C%NY==NY-1 && n_C>NY-1 && n_C<NX*NY-1)
          {
             mm1=CC->DOF_value(0,n_C-NY-1)+2.0*CC->DOF_value(0,n_C-NY)+CC->DOF_value(0,n_C-NY);
             mm2=CC->DOF_value(0,n_C+NY-1)+2.0*CC->DOF_value(0,n_C+NY)+CC->DOF_value(0,n_C+NY);
             mx=mm1-mm2;
             mm1=CC->DOF_value(0,n_C-NY-1)+2.0*CC->DOF_value(0,n_C-1)+CC->DOF_value(0,n_C+NY-1);
             mm2=CC->DOF_value(0,n_C-NY)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C+NY);
             mz=(mm1-mm2)*2.0;
          }
          else if (n_C==NY-1)
          {
             mm1=CC->DOF_value(0,n_C-1)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C);
             mm2=CC->DOF_value(0,n_C+NY-1)+2.0*CC->DOF_value(0,n_C+NY)+CC->DOF_value(0,n_C+NY);
             mx=mm1-mm2;
             mm1=CC->DOF_value(0,n_C-1)+2.0*CC->DOF_value(0,n_C-1)+CC->DOF_value(0,n_C+NY-1);
             mm2=CC->DOF_value(0,n_C)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C+NY);
             mz=mm1-mm2;
          }
          else if (n_C==0)
          {
             mm1=CC->DOF_value(0,n_C)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C+1);
             mm2=CC->DOF_value(0,n_C+NY)+2.0*CC->DOF_value(0,n_C+NY)+CC->DOF_value(0,n_C+NY+1);
             mx=mm1-mm2;
             mm1=CC->DOF_value(0,n_C)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C+NY);
             mm2=CC->DOF_value(0,n_C+1)+2.0*CC->DOF_value(0,n_C+1)+CC->DOF_value(0,n_C+NY+1);
             mz=mm1-mm2;
          }
          else if (n_C==(NX-1)*NY)
          {
             mm1=CC->DOF_value(0,n_C-NY)+2.0*CC->DOF_value(0,n_C-NY)+CC->DOF_value(0,n_C-NY+1);
             mm2=CC->DOF_value(0,n_C)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C+1);
             mx=mm1-mm2;
             mm1=CC->DOF_value(0,n_C-NY)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C);
             mm2=CC->DOF_value(0,n_C-NY+1)+2.0*CC->DOF_value(0,n_C+1)+CC->DOF_value(0,n_C+1);
             mz=mm1-mm2;
          }
          else if (n_C==NX*NY-1)
          {
             mm1=CC->DOF_value(0,n_C-NY-1)+2.0*CC->DOF_value(0,n_C-NY)+CC->DOF_value(0,n_C-NY);
             mm2=CC->DOF_value(0,n_C-1)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C);
             mx=mm1-mm2;
             mm1=CC->DOF_value(0,n_C-NY-1)+2.0*CC->DOF_value(0,n_C-1)+CC->DOF_value(0,n_C+NY-1);
             mm2=CC->DOF_value(0,n_C-NY)+2.0*CC->DOF_value(0,n_C)+CC->DOF_value(0,n_C);
             mz=mm1-mm2;
          }
          else
          {
             mm1=CC->DOF_value(0,n_C-NY-1)+2.0*CC->DOF_value(0,n_C-NY)+CC->DOF_value(0,n_C-NY+1);
             mm2=CC->DOF_value(0,n_C+NY-1)+2.0*CC->DOF_value(0,n_C+NY)+CC->DOF_value(0,n_C+NY+1);
             mx=mm1-mm2;
             mm1=CC->DOF_value(0,n_C-NY-1)+2.0*CC->DOF_value(0,n_C-1)+CC->DOF_value(0,n_C+NY-1);
             mm2=CC->DOF_value(0,n_C-NY+1)+2.0*CC->DOF_value(0,n_C+1)+CC->DOF_value(0,n_C+NY+1);
             mz=mm1-mm2;
          }          
          inv = 0;
          if (mz<0.0)
          {
             mm1=-s1;
             s1=-s2;
             s2=mm1;
             mz=-mz;
             inv = 1;
          }
          mz=mz+1.E-50;
          mx=PEL::abs(mx)+1.E-50;
          mm2= mx > mz ? mx : mz;
          mx=mx/mm2;
          mz=mz/mm2;

          mm1= mx < mz ? mx : mz;
          V1=0.5*mm1;
          V2=1.0-V1;
          if (CC->DOF_value(0,n_C)<=V1)
             alfa=sqrt(2.0*CC->DOF_value(0,n_C)*mm1);
          else if (CC->DOF_value(0,n_C)<=V2)
             alfa=CC->DOF_value(0,n_C)+0.5*mm1;
          else 
             alfa=mm1+1.0-sqrt(2.0*(1.0-CC->DOF_value(0,n_C))*mm1);
          mz=mz/(1.0-s1+s2);
          alfa=alfa+mz*s1;
          
          mm1=(-s1 > 0.0) ? -s1 : 0.0;
          mm2=alfa+mz*mm1;
          CC->set_DOF_value(1, n_C, vol2(mz,mx,mm2,mm1) );

          mm1=1.0-(s1>0.0?s1:0.0)+(s2<0.0?s2:0.0);
          mm2=alfa-mz*(s1>0.0?s1:0.0);
          CC->set_DOF_value(2, n_C, vol2(mz,mx,mm2,mm1) );
          
          mm1=s2>0.0?s2:0.0;
          mm2=alfa-mz;
          CC->set_DOF_value(3, n_C, vol2(mz,mx,mm2,mm1) );

          if(inv==1)
          {
             mm1=CC->DOF_value(1,n_C);
             CC->set_DOF_value(1, n_C, CC->DOF_value(3,n_C) );
             CC->set_DOF_value(3, n_C, mm1 );
          }

          CC->set_DOF_value(1,n_C,CC->DOF_value(1,n_C)*Yvector(n_C%NY)/(n_C%NY==0 ? Yvector(n_C%NY) : Yvector(n_C%NY-1)));
          CC->set_DOF_value(3,n_C,CC->DOF_value(3,n_C)*Yvector(n_C%NY)/(n_C%NY==NY-1 ? Yvector(n_C%NY) : Yvector(n_C%NY+1)));
      }
   }
   //update C
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      n_C = cFE->global_node( CC, 0 ) ;
      double u,d;
      u = (n_C%NY==0) ? 0.0 : CC->DOF_value(3,n_C-1);
      d = (n_C%NY==NY-1) ? 0.0 : CC->DOF_value(1,n_C+1);
      CC->set_DOF_value(0, n_C, u+CC->DOF_value(2,n_C)+d );
   }
   //regularize C 
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      n_C = cFE->global_node( CC, 0 ) ;
      mm1=0.0>CC->DOF_value(0,n_C)?0.0:CC->DOF_value(0,n_C);
      mm2=1.0<mm1?1.0:mm1;
      CC->set_DOF_value(0, n_C, mm2);
   }
   for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
       n_C = cFE->global_node( CC, 0 ) ;
       temp = CC->DOF_value(0,n_C);
       if (CC->DOF_value(4, n_C)>1.0 && temp < 1.0)   //previous C2 or interface; current not C1 so it has to be C2 or interface
          CC->set_DOF_value(0, n_C, 2.0-temp);        
       //if (CC->DOF_value(4, n_C)>1.0 && temp == 1.0)   //previous C2; current C1: do nothing!
          //CC->set_DOF_value(0, n_C, 2.0-temp);
       //if (CC->DOF_value(4, n_C)==1.0 && temp == 1.0)  //previous C1; current C1: do nothing!
          //CC->set_DOF_value(0, n_C, 2.0-temp);
       if (CC->DOF_value(4, n_C)==1.0 && temp < 1.0)   //previous C1; current interface(C0 or C2)
       { 
          //check the left cell(closer to symmetry axis) C=0: IC0; 0<C<1: IC0; C=1: IC2; 1<C<2: IC2; C=2:unlikely
          if (CC->DOF_value(4, n_C-NY)>=1.0)
             CC->set_DOF_value(0, n_C, 2.0-temp); 
       }
       if (n_C%NY==0)
          CC->set_DOF_value(0, n_C, CC->DOF_value(4, n_C)); 
   }
}

//-------------------------------------------------------------------------
void
P3_PLIC_Advection::do_inner_iterations_stage( FE_TimeIterator const* t_it )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "P3_PLIC_Advection:: do_inner_iterations_stage" ) ;

   if( useCOURANT ){
	double epsilon= -1.e-6;
	double dt_new = t_it->time_step();
      t_local =  FE_LocalTimeIteratorAdapter::make_default(t_global, DOM, PRMS);
      t_global->start();
      t_global->set_time_offset( t_it->time()-t_it->time_step() );
	
      dt_new = time_step_by_Courant( dt_new, t_it);
	  t_local->restart_iteration_with_new_time_step(dt_new);
      do {
	    PEL::out() << indent() << "- local time step = " << t_global->time() << std::endl;
	    t_local->initialize_time_step();
	    t_local->propose_next_time_step( dt_new );
	    t_local->adapt_time_iterator();

            do_one_inner_iteration( t_global ) ;

		t_global->go_next_time( t_it->time() + epsilon); //go_next_time(double tend)
      }
      while ( !AD_inner_iterations_are_completed( t_global, t_it, epsilon ) );

      delete t_local;
   }
   else
       do_one_inner_iteration( t_it ) ;
}


//-------------------------------------------------------------------------
bool
P3_PLIC_Advection::AD_inner_iterations_are_completed( FE_TimeIterator * t_glo, FE_TimeIterator const* t_it, double epsilon ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "P3_PLIC_Advection:: AD_inner_iterations_are_completed( t_global, t_it )" ) ;

   if ( useCOURANT )
	{
	   if( t_glo->time() <= t_it->time()+epsilon )
        {
           return false ;
        }
	  else
	  {
	     t_glo->finish_iterations();
	     return true;
	  }
	}
   else
	return true;

}

//----------------------------------------------------------------------
void
P3_PLIC_Advection:: reset_discrete_problem( FE_TimeIterator const* t_it )
//----------------------------------------------------------------------
{
   PEL_LABEL( "P3_PLIC_Advection:: reset_discrete_problem" ) ;
}



//-------------------------------------------------------------------------
void
P3_PLIC_Advection::update_fields( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "P3_PLIC_Advection:: update_fields" ) ;
}

//-------------------------------------------------------------------------
double
P3_PLIC_Advection::smallest_vertices( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "P3_PLIC_Advection:: smallest_vertices" ) ;

   PEL_Communicator const* com = PEL_Exec::communicator() ;

   double min_distance = PEL::bad_double() ;;
   for ( sFE->start() ; sFE->is_valid(); sFE->go_next() )
   {
	   double d_KL = sFE->distance_to_adjacent_finite_volume_center( 0 ) +
	                 sFE->distance_to_adjacent_finite_volume_center( 1 );

       if( PEL::abs(d_KL) < min_distance ) min_distance = PEL::abs(d_KL) ;
   }

   min_distance = com->min(min_distance);
   return (min_distance);
}


//-------------------------------------------------------------------------
double
P3_PLIC_Advection::time_step_by_Courant( double dt_new, FE_TimeIterator const* t_it )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "P3_PLIC_Advection:: time_step_by_Courant" ) ;

   double min_distance = smallest_vertices();
   PEL_Communicator const* com = PEL_Exec::communicator() ;

   for( size_t i=0 ; i<FIELDS_TABLE.size() ; ++i )
   {
      PDE_DiscreteField const* f = FIELDS->item( FIELDS_TABLE(i) ) ;
      doubleVector v_max( f->nb_components() ) ;
      double cfl_tmp=0.;
      for( size_t ic=0 ; ic<f->nb_components() ; ++ic )
      {
         v_max(ic) = -PEL::bad_double() ;
         for( size_t i_node=0 ; i_node<f->nb_nodes() ; ++i_node )
         {
             double const v_cur = f->DOF_value( LEVEL0, i_node, ic ) ;
             if( PEL::abs(v_cur) > v_max(ic) ) v_max(ic) = PEL::abs(v_cur) ;
         }
         v_max(ic) = com->max( v_max(ic) ) ;
         cfl_tmp += v_max(ic) ;
      }

	  bool change_dt = false;

	  double cfl = cfl_tmp * dt_new/min_distance; // initial value
      while ( cfl > COURANT ) {
         dt_new = 0.9 * dt_new ;
         cfl = cfl_tmp * dt_new/min_distance;
         change_dt = true;
      }

	  if(change_dt) {
	     PEL::out() << indent() << "Given Courant number = " << COURANT << ", (grid) min distance = " << min_distance << std::endl;
	     PEL::out() << indent() << "Time step by Courant number = " <<  dt_new << std::endl;
         if (dt_new < 1.E-06) {
                PEL::out() << "*****************************************" << endl;                
                PEL::out() << "**  New time step too small (<1.E-6)!  **" << endl;
                PEL::out() << "*****************************************" << endl;   
                PEL_Error::exit();
         }
	     int tmp = int(t_it->time_step()/dt_new)+1 ;
	     dt_new = t_it->time_step()/tmp ;
	     PEL::out() << indent() << "Number of local time steps " << tmp << std::endl;
	     PEL::out() << indent() << "Adjusted new time step to fit time interval = " <<  dt_new << std::endl;
	  }
   }
   return( dt_new );
}



//Given mx*x+mz*z=alfa and flux, compute the volume change in the cell
double
P3_PLIC_Advection::vol2( double mx, double mz, double alfa, double b)
//----------------------------------------------------------------------
{
   PEL_LABEL( "P3_PLIC_Advection::vol2" ) ;
   double alfaa,mb;
   mb=mx*b;
   if (alfa<=0.0)
      return 0.0;
   else if (alfa>=(mb+mz) || b==0.0)
      return b;
   else
   {
      if (mb<=mz)
      {
         if (alfa<mb)
            return 0.5*alfa*alfa/(mx*mz);
         else if (alfa>=mz)
         {
            alfaa=mb+mz-alfa;
            return b-0.5*alfaa*alfaa/(mx*mz);
         }
         else
            return b*(alfa-0.5*mb)/mz;
      }
      else
      {
         if (alfa<mz)
            return 0.5*alfa*alfa/(mx*mz);
         else if (alfa>mb)
         {
            alfaa=mb+mz-alfa;
            return b-0.5*alfaa*alfaa/(mx*mz);
         }
         else
            return (alfa-0.5*mz)/mx;
      }
    }           
}


void P3_PLIC_Advection::update_interface(double diff)
//----------------------------------------------------------------------
{
   PEL_LABEL( "P3_PLIC_Advection::update_interface" ) ;
   double step,ratio,C,volume=0.0;
   size_t n_C;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      n_C = cFE->global_node( CC, 0 ) ;
      C = CC->DOF_value(0,n_C) ;
      step = (C<Cmax ? 1.0 : 0.0);
      if ( FE::geometry() == FE::axisymmetrical )
         volume+=C*step*Xcoord(n_C/NY);
      else
         volume+=C*step;
      //num+=step;
   }
   if (volume==0.0)
      return;
   ratio = 1.0-diff/volume;
   PEL::out() << indent() << "interface shrink ratio = " <<  ratio << std::endl;
   for( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
   {
      n_C = cFE->global_node( CC, 0 ) ;
      C = CC->DOF_value(0,n_C) ;
      if (C<Cmax)  //one choice of interface
      CC->set_DOF_value(0, n_C, CC->DOF_value(0,n_C)*ratio);
   }
}
