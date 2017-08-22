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

#include <PDE_SetOfBCs.hh>

#include <GE_Color.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_SetOfDiscreteFields.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_KeywordDataIterator.hh>
#include <PEL_ModuleIterator.hh>
#include <PEL_Root.hh>
#include <PEL_Vector.hh>
#include <stringVector.hh>

size_t const
PDE_SetOfBCs:: all_components = (size_t)-1 ;

class PDE_BC : public PEL_Object
{
   public:

      PDE_BC( PEL_Object* a_owner,
              PDE_DiscreteField const* f,
              GE_Color const* col,
              size_t comp,
              PEL_ModuleExplorer const* exp = 0 ) ;
     ~PDE_BC( void ) ;

      void re_initialize(
              PDE_DiscreteField const* f,
              GE_Color const* col,
              size_t comp,
              PEL_ModuleExplorer const* exp = 0 ) ;

      PEL_ModuleExplorer const* explorer( void ) const ;

      bool is_equal( PEL_Object const* obj ) const ;
      int three_way_comparison( PEL_Object const* obj ) const ;
      size_t hash_code( void ) const ;

   private:
      
      PDE_DiscreteField const* FIELD  ;
      GE_Color const * COLOR ;
      size_t COMP ;
      PEL_ModuleExplorer const* EXP ;
} ;

static PDE_BC* DUMMY_BC =
    new PDE_BC( PEL_Root::object(),  0, 0, PDE_SetOfBCs::all_components ) ;
      
//-------------------------------------------------------------------------
PDE_SetOfBCs const*
PDE_SetOfBCs:: create( PEL_Object* a_owner,
                       PEL_ModuleExplorer const* exp,
                       PDE_SetOfDiscreteFields const* fields )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBCs:: create" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;
   PEL_CHECK_PRE( fields != 0 ) ;
   
   PDE_SetOfBCs const* result = new PDE_SetOfBCs( a_owner, exp, fields ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PDE_SetOfBCs:: PDE_SetOfBCs( PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp,
                             PDE_SetOfDiscreteFields const* fields )
//-------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , BCS( PEL_List::create( this ) )
   , MACRO_COLS( PEL_Vector::create( this, 0 ) )
   , MACRO_TYPE( 0 )
{
   if( exp->has_module( "boundary_conditions" ) &&
       exp->has_module( "macro_boundary_conditions" ) )
   {
      PEL_Error::object()->raise_module_error(
         exp,
         "modules \"boundary_conditions\" and\n"
         "        \"macro_boundary_conditions\"\n"
         "can not be defined together" ) ;
   }
   if( exp->has_module( "boundary_conditions" ) )
   {
      PEL_ModuleExplorer* bcTree =
         exp->create_subexplorer( 0, "boundary_conditions" ) ;
      explore( bcTree, fields ) ;
      bcTree->destroy() ; bcTree = 0 ;
   }
   else if( exp->has_module( "macro_boundary_conditions" ) )
   {
      PEL_ModuleExplorer* macro_bcTree =
         exp->create_subexplorer( 0, "macro_boundary_conditions" ) ;
      for( macro_bcTree->start_module_iterator();
           macro_bcTree->is_valid_module() ;
           macro_bcTree->go_next_module() )
      {
         PEL_ModuleExplorer* bcTree = macro_bcTree->create_subexplorer( 0 ) ;
         GE_Color const* col =
                         GE_Color::object( bcTree->string_data( "color" ) ) ;
         std::string const& type = bcTree->string_data( "type" ) ;
         MACRO_COLS->append( const_cast<GE_Color*>( col ) ) ;
         MACRO_TYPE.append( type ) ;
         explore( bcTree, fields, col ) ;
         bcTree->destroy() ; bcTree = 0 ;
      }
      macro_bcTree->destroy() ; macro_bcTree = 0 ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
}

//-------------------------------------------------------------------------
PDE_SetOfBCs:: ~PDE_SetOfBCs( void )
//-------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
}

//-------------------------------------------------------------------------
bool
PDE_SetOfBCs:: has_BC( GE_Color const* color,
                       PDE_DiscreteField const* field,
                       size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBCs:: has_BC" ) ;
   PEL_CHECK_PRE( color != 0 ) ;
   PEL_CHECK_PRE( field != 0 ) ;
   PEL_CHECK_PRE( ic==all_components || ic<field->nb_components() ) ;
   PEL_CHECK_INV( invariant() ) ;

   DUMMY_BC->re_initialize( field, color, ic ) ;
   return( BCS->has( DUMMY_BC ) ) ;
}

//-------------------------------------------------------------------------
PEL_ModuleExplorer const*
PDE_SetOfBCs:: BC_explorer( GE_Color const* color,
                            PDE_DiscreteField const* field,
                            size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBCs:: BC_explorer" ) ;
   PEL_CHECK_PRE( color != 0 ) ;
   PEL_CHECK_PRE( field != 0 ) ;
   PEL_CHECK_PRE( ic==all_components || ic<field->nb_components() ) ;
   PEL_CHECK_PRE( has_BC( color, field, ic ) ) ;
   PEL_CHECK_INV( invariant() ) ;

   DUMMY_BC->re_initialize( field, color, ic ) ;
   PEL_ModuleExplorer const* result =
         static_cast<PDE_BC const*>( BCS->item( DUMMY_BC ) )->explorer() ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;   
}

//-------------------------------------------------------------------------
bool
PDE_SetOfBCs:: has_macro_BC( GE_Color const* color ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBCs:: has_macro_BC" ) ;
   PEL_CHECK_PRE( color != 0 ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( macro_BC_index( color ) != PEL::bad_index() ) ;
}

//-------------------------------------------------------------------------
std::string const&
PDE_SetOfBCs:: macro_BC( GE_Color const* color ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBCs:: macro_BC" ) ;
   PEL_CHECK_PRE( color != 0 ) ;
   PEL_CHECK_PRE( has_macro_BC( color ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( MACRO_TYPE( macro_BC_index( color ) ) ) ;
}

//-------------------------------------------------------------------------
void
PDE_SetOfBCs:: explore( PEL_ModuleExplorer* bcTree,
                        PDE_SetOfDiscreteFields const* fields,
                        GE_Color const* bcTree_color ) 
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBCs:: explore" ) ;
   
   PEL_CHECK_INV( invariant() ) ;
      
   for( bcTree->start_module_iterator();
        bcTree->is_valid_module() ;
        bcTree->go_next_module() )
   {
      PEL_ModuleExplorer const* exp = bcTree->create_subexplorer( this ) ;
      std::string const& f_name = exp->string_data( "field" ) ;
      if( !fields->has( f_name ) )
      {
         PEL_Error::object()->raise_data_error(
            exp, "field",
            "field of name \""+f_name+"\" does not exist" ) ;
      }
      PDE_DiscreteField const* field = fields->item( f_name ) ;

      GE_Color const* color = 0 ;
      int component = all_components ;
      if( bcTree_color != 0 )
      {
         color = bcTree_color ;
      }
      else
      {
         if( exp->has_entry( "color" ) )
         {
            color = GE_Color::object( exp->string_data( "color" ) ) ;
         }
      }
      if( exp->has_entry( "component" ) )
      {
         component = exp->int_data( "component" ) ;
         if( component<0 || component>=(int) field->nb_components() )
         {
            PEL_Error::object()->raise_data_error(
               exp, "component",
               "bad value (see number of components of field of name \""+f_name+"\")" ) ;
            
         }
      }
      PDE_BC* bc = new PDE_BC( BCS, field, color, component, exp ) ;
      if( BCS->has( bc ) )
      {
         PEL_Error::object()->raise_module_error(
            exp, "Duplicate boundary condition" ) ;
      }
      BCS->append( bc ) ; 
   }
   
   PEL_CHECK_INV( invariant() ) ;
}

//-------------------------------------------------------------------------
size_t
PDE_SetOfBCs:: macro_BC_index( GE_Color const* color ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBCs:: macro_BC_index" ) ;

   size_t result = PEL::bad_index() ;
   bool found = false ;
   for( size_t i=0 ; !found && i<MACRO_COLS->index_limit() ; ++i )
   {
      GE_Color const* col = static_cast<GE_Color const*>( MACRO_COLS->at(i) ) ;
      found = color->is_overlapping( col ) ;
      if( found ) result = i ;
   }
   
   return( result ) ;
}

//-------------------------------------------------------------------------
bool
PDE_SetOfBCs:: invariant( void ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   PEL_ASSERT( BCS!=0 ) ;
   PEL_ASSERT( MACRO_TYPE.size() == MACRO_COLS->index_limit() ) ;
   return( true ) ;
}

//internal-----------------------------------------------------------------
PDE_BC:: PDE_BC( PEL_Object* a_owner,
                 PDE_DiscreteField const* f,
                 GE_Color const* col,
                 size_t comp,
                 PEL_ModuleExplorer const* exp )
//internal-----------------------------------------------------------------
   : PEL_Object(a_owner)
   , FIELD( f )
   , COLOR( col )
   , COMP( comp )
   , EXP( exp )
{
}

//internal-----------------------------------------------------------------
PDE_BC:: ~PDE_BC( void )

//internal-----------------------------------------------------------------
{
}

//internal-----------------------------------------------------------------
void PDE_BC:: re_initialize( PDE_DiscreteField const* f,
                             GE_Color const* col,
                             size_t comp,
                             PEL_ModuleExplorer const* exp )
//internal-----------------------------------------------------------------
{
   FIELD = f ;
   COLOR = col ;
   COMP = comp ;
   EXP = exp ;
}

//internal-----------------------------------------------------------------
PEL_ModuleExplorer const*
PDE_BC:: explorer( void ) const
//internal-----------------------------------------------------------------
{
   return( EXP ) ;
}

//internal-----------------------------------------------------------------
bool
PDE_BC:: is_equal( PEL_Object const* obj ) const
//internal-----------------------------------------------------------------
{
   PDE_BC const* bc = static_cast<PDE_BC const* >(obj) ;
   
   bool result = ( FIELD->id_number() == bc->FIELD->id_number() ) ;
   if( result && COLOR!=0 )
   {
      result = ( bc->COLOR == 0 || COLOR->is_overlapping( bc->COLOR ) ) ;
   } 
   if( result && COMP != PDE_SetOfBCs::all_components )
   {
      result = ( COMP == bc->COMP ) ;
   }
   return( result ) ;
}

//internal-----------------------------------------------------------------
int
PDE_BC:: three_way_comparison( PEL_Object const* ) const
//internal-----------------------------------------------------------------
{
   PEL_Error::object()->raise_not_implemented( this, "three_way_comparison" ) ;
   return( 0 ) ;
}

//internal-----------------------------------------------------------------
size_t
PDE_BC:: hash_code( void ) const
//internal-----------------------------------------------------------------
{
   PEL_Error::object()->raise_not_implemented( this, "hash_code" ) ;
   return( 0 ) ;
}






