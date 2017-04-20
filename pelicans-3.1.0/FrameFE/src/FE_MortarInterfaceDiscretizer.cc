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

#include <FE_MortarInterfaceDiscretizer.hh>

#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <LA_Matrix.hh>
#include <LA_Vector.hh>

#include <GE_QRprovider.hh>

#include <PDE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_InterfaceAndFields.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_LocalEquation.hh>
#include <PDE_LocalFEmortarSide.hh>
#include <PDE_SetOfDomains.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PDE_SystemNumbering.hh>

#include <FE.hh>
#include <FE_OneStepIterationOpen.hh>

#include <iostream>
#include <sstream>

using std::cout ;
using std::endl ;
using std::ostringstream ;
using std::string ;

struct FE_MortarInterfaceDiscretizer_ERROR 
{
   static void n0( std::string const& fname, 
                   PEL_ModuleExplorer const* exp ) ;
   static void n1( size_t n_interf, size_t n_dom0, size_t n_dom1,
                   PEL_ModuleExplorer const* exp ) ; 
} ;

//----------------------------------------------------------------------------
FE_MortarInterfaceDiscretizer*
FE_MortarInterfaceDiscretizer:: create( PEL_Object* a_owner,
                                        PDE_InterfaceAndFields const* interf,
                                        PEL_Vector const* dom_discs,
                                        PEL_ModuleExplorer* exp )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_MortarInterfaceDiscretizer:: create" ) ;
   PEL_CHECK_PRE( interf != 0 ) ;
   PEL_CHECK_PRE( dom_discs != 0 ) ;
   PEL_CHECK_PRE( dom_discs->index_limit() >= 2 ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   PEL_CHECK_PRE( 
      FORALL( ( size_t i=0 ; i<dom_discs->index_limit() ; ++i ),
        dynamic_cast< FE_OneStepIterationOpen* >( dom_discs->at( i ) ) != 0 ));
   
   FE_MortarInterfaceDiscretizer* result = 
      new FE_MortarInterfaceDiscretizer( a_owner, interf, dom_discs, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
FE_MortarInterfaceDiscretizer:: FE_MortarInterfaceDiscretizer( 
                                         PEL_Object* a_owner,
                                         PDE_InterfaceAndFields const* interf,
                                         PEL_Vector const* dom_discs,
                                         PEL_ModuleExplorer* exp )
//----------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , INTERF( interf )
   , ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , T_ELEMENT_EQ( PDE_LocalEquation::create( this ) )
   , msFE( interf->create_LocalFEmortarSide( this ) )
   , QRP( GE_QRprovider::object( 
                         exp->string_data( "quadrature_rule_provider" ) ) )
{
   PEL_LABEL( "FE_MortarInterfaceDiscretizer:: FE_MortarInterfaceDiscretizer" ) ;
   
   PDE_DomainAndFields const* dom0 = interf->adjacent_domain( 0 ) ;
   PDE_DomainAndFields const* dom1 = interf->adjacent_domain( 1 ) ;
   
   FE_OneStepIterationOpen* disc0 = 0 ;
   FE_OneStepIterationOpen* disc1 = 0 ;
   for( size_t i=0 ; i<dom_discs->index_limit() ; ++i )
   {
      FE_OneStepIterationOpen* disc = 
         static_cast< FE_OneStepIterationOpen* >( dom_discs->at( i ) ) ;
      if( disc->domain() == dom0 )
      {
         PEL_ASSERT( disc0 == 0 ) ;
         disc0 = disc ;
      }
      else if( disc->domain() == dom1 )
      {
         PEL_ASSERT( disc1 == 0 ) ;
         disc1 = disc ;
      }
   }
   PEL_ASSERT( disc0 != 0 && disc1 != 0  ) ;
   
   PEL_ModuleExplorer* e = exp->create_subexplorer( 0, "unknowns" ) ;
   e->start_module_iterator() ;
   for( ; e->is_valid_module() ; e->go_next_module() )
   {
      PEL_ModuleExplorer* ee = e->create_subexplorer( 0 ) ;
      PDE_DiscreteField* ff = INTERF->set_of_discrete_fields()->item( 
                                    ee->string_data( "interface_field" ) ) ;
      msFE->require_field_calculation( ff, PDE_LocalFE::N ) ;
      LLs.push_back( ff ) ;
      LL_links.push_back( PDE_LinkDOF2Unknown::create( this, ff,
                                  "sequence_of_the_components",true ) ) ;

      {
         string const& nn = ee->string_data( "field_of_adjacent_domain_0" ) ;
         size_t idx = disc0->index_of_field( nn ) ;
         if( idx == PEL::bad_index() )
            FE_MortarInterfaceDiscretizer_ERROR::n0( nn, ee ) ;
         IDs_0.push_back( idx ) ;
         PDE_DiscreteField const* dom_ff = disc0->field( idx ) ;
         check_field_consistency( ff, dom_ff, ee ) ;
         msFE->require_field_calculation( dom_ff, PDE_LocalFE::N ) ;
         UUs_0.push_back(  dom_ff ) ;
      }
      {
         string const& nn = ee->string_data( "field_of_adjacent_domain_1" ) ;
         size_t idx = disc1->index_of_field( nn ) ;
         if( idx == PEL::bad_index() )
            FE_MortarInterfaceDiscretizer_ERROR::n0( nn, ee ) ;
         IDs_1.push_back( idx ) ;   
         PDE_DiscreteField const* dom_ff = disc1->field( idx ) ;
         check_field_consistency( ff, dom_ff, ee ) ;
         msFE->require_field_calculation( dom_ff, PDE_LocalFE::N ) ;
         UUs_1.push_back( dom_ff ) ;
      }
      ee->destroy() ; ee = 0 ;
   }
   e->destroy() ; e = 0 ;
   
   if( nb_unknowns() != disc0->nb_unknowns() ||
       nb_unknowns() != disc1->nb_unknowns() )
   { 
      FE_MortarInterfaceDiscretizer_ERROR::n1( nb_unknowns(),
                                               disc0->nb_unknowns(),
                                               disc1->nb_unknowns(), exp ) ;
   }
}

//----------------------------------------------------------------------------
FE_MortarInterfaceDiscretizer:: ~FE_MortarInterfaceDiscretizer( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
PDE_InterfaceAndFields const*
FE_MortarInterfaceDiscretizer:: interface( void ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_MortarInterfaceDiscretizer:: interface" ) ;
   
   return( INTERF ) ;
}

//----------------------------------------------------------------------------
size_t
FE_MortarInterfaceDiscretizer:: nb_unknowns( void ) const
//----------------------------------------------------------------------------
{
   return( LLs.size() ) ;
}

//----------------------------------------------------------------------------
PDE_DiscreteField*
FE_MortarInterfaceDiscretizer:: field( size_t i_unk ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_MortarInterfaceDiscretizer:: field" ) ;
   PEL_CHECK_PRE( i_unk < nb_unknowns() ) ;
   
   PDE_DiscreteField* result = LLs[i_unk] ;
   
   PEL_CHECK_POST( result != 0 ) ;
//   PEL_CHECK_POST( result == link_DOF_2_unknown( i_unk )->field() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PDE_LinkDOF2Unknown* 
FE_MortarInterfaceDiscretizer:: create_link_DOF_2_unknown( PEL_Object* a_owner,
                                                           size_t i_unk ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "FE_MortarInterfaceDiscretizer:: create_link_DOF_2_unknown" ) ;
   PEL_CHECK_PRE( i_unk < nb_unknowns() ) ;
   
   PDE_LinkDOF2Unknown* result = 
         PDE_LinkDOF2Unknown::create( a_owner, 
                                      LLs[i_unk],
                                      "sequence_of_the_components", 
                                      true ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->field( ) == field( i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
FE_MortarInterfaceDiscretizer:: assemble_contribution(
                                         LA_Matrix* matrix,
                                         LA_Vector* vector,
                                         PDE_SystemNumbering const* nmb,
                                         size_t interf_shift, 
                                         size_t domain_0_shift, 
                                         size_t domain_1_shift ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_MortarInterfaceDiscretizer:: assemble_contribution" ) ;
   PEL_CHECK_PRE( matrix != 0 ) ;
   PEL_CHECK_PRE( vector != 0 ) ;
   PEL_CHECK_PRE( nmb != 0 ) ;
   //??? many many other preconditions

   for( size_t i=0 ; i<LLs.size() ; ++i )
   {
      size_t i_eq = interf_shift + i ;
      size_t i_unk_0 = domain_0_shift + IDs_0[i] ;
      size_t i_unk_1 = domain_1_shift + IDs_1[i] ;
      
      for( msFE->start() ; msFE->is_valid() ; msFE->go_next() )
      {
         msFE->set_row_and_col_fields( LLs[i], UUs_0[i] ) ; 
         ELEMENT_EQ->initialize( msFE->row_field_node_connectivity(), 1,
                                 msFE->col_field_node_connectivity(), 1 ) ;
         msFE->start_IP_iterator( QRP ) ;
         for( ; msFE->valid_IP() ; msFE->go_next_IP() )
         {
            FE::add_row_col_NS( ELEMENT_EQ, msFE, +1.0 ) ;
         }
         PDE::assemble_in_matrix_vector_0( matrix, vector, ELEMENT_EQ, 
                                           nmb, i_eq, i_unk_0 ) ;
      
         T_ELEMENT_EQ->set_as_transpose( ELEMENT_EQ ) ;
         PDE::assemble_in_matrix_vector_0( matrix, vector, T_ELEMENT_EQ, 
                                           nmb, i_unk_0, i_eq ) ;

         msFE->set_row_and_col_fields( LLs[i], UUs_1[i] ) ;
         ELEMENT_EQ->initialize( msFE->row_field_node_connectivity(), 1,
                                 msFE->col_field_node_connectivity(), 1 ) ;
         msFE->start_IP_iterator( QRP ) ;
         for( ; msFE->valid_IP() ; msFE->go_next_IP() )
         {
            FE::add_row_col_NS( ELEMENT_EQ, msFE, -1.0 ) ;
         }
         PDE::assemble_in_matrix_vector_0( matrix, vector, ELEMENT_EQ, 
                                           nmb, i_eq, i_unk_1 ) ;
      
         T_ELEMENT_EQ->set_as_transpose( ELEMENT_EQ ) ;
         PDE::assemble_in_matrix_vector_0( matrix, vector, T_ELEMENT_EQ, 
                                           nmb, i_unk_1, i_eq ) ;
      }
   }
}

//---------------------------------------------------------------------------
void
FE_MortarInterfaceDiscretizer:: check_field_consistency( 
                                      PDE_DiscreteField const* ff,
                                      PDE_DiscreteField const* dom_ff,
                                      PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
{
   if( ff->nb_components() != dom_ff->nb_components() )
   {
      ostringstream mesg ;
      mesg << "In module:" << endl ;
      mesg << exp->absolute_path_name() << endl ;
      mesg << "the fields of name \"" << ff->name() 
           << "\" and \"" << dom_ff->name() << "\"" << endl ;
      mesg << "are incompatible" ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
      
   }
}

//internal-------------------------------------------------------------------
void
FE_MortarInterfaceDiscretizer_ERROR:: n0( std::string const& fname,
                                          PEL_ModuleExplorer const* exp )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "In module:" << endl ;
   mesg << exp->absolute_path_name() << endl ;
   mesg << "the field of name \"" << fname << "\" is unknown" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
FE_MortarInterfaceDiscretizer_ERROR:: n1( size_t n_interf, 
                                          size_t n_dom0, 
                                          size_t n_dom1,
                                          PEL_ModuleExplorer const* exp )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "In module:" << endl ;
   mesg << exp->absolute_path_name() << endl ;
   mesg << "there should be the same number of unknowns associated" << endl ;
   mesg << "with the interface and with the adjacent domains" << endl ;
   mesg << "Number of unknowns declared" << endl ;
   mesg << "   * in the interface: " << n_interf << endl ;
   mesg << "   * in the first  adjacent domain: " << n_dom0 << endl ;
   mesg << "   * in the second adjacent domain: " << n_dom1 ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
   
}


