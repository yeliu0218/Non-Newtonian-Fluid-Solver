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

#include <PDE_LinkDOF2Unknown.hh>

#include <PEL.hh>
#include <PEL_Error.hh>

#include <LA_SeqVector.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_CrossProcessUnknownNumbering.hh>

#include <iostream>

#ifdef OUTLINE
   #define inline
   #include <PDE_LinkDOF2Unknown.icc>
   #undef inline
#endif

using std::string ;

//----------------------------------------------------------------------
PDE_LinkDOF2Unknown*
PDE_LinkDOF2Unknown:: create( PEL_Object* a_owner,
                              PDE_DiscreteField const* ff,
                              std::string const& ordering,
                              bool imposed_out )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LinkDOF2Unknown:: create" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;
   PEL_CHECK_PRE( ordering == "sequence_of_the_components" || 
                  ordering == "sequence_of_the_nodes" ) ;

   PDE_LinkDOF2Unknown* result = new PDE_LinkDOF2Unknown( a_owner, 
                                                          ff, 
                                                          ordering,
                                                          imposed_out ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->field() == ff ) ;
   PEL_CHECK_POST( result->nb_field_nodes() == ff->nb_nodes() ) ;
   PEL_CHECK_POST( result->components_table().size() == ff->nb_components() ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<ff->nb_components() ; i++ ),
                           result->components_table()(i)==i ) ) ;
   PEL_CHECK_POST( result->DOFs_with_imposed_value_are_dropped() == 
                                                               imposed_out ) ;
   PEL_CHECK_POST( result->DOFs_ordering_in_unknown() == ordering ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_LinkDOF2Unknown:: PDE_LinkDOF2Unknown( PEL_Object* a_owner,
                                           PDE_DiscreteField const* ff,
                                           std::string const& ordering, 
                                           bool imposed_out )
//----------------------------------------------------------------------
  : PEL_Object( a_owner )
  , FIELD( ff )
  , NB_COMPS( ff->nb_components() )
  , COMPS( ff->nb_components() )
  , INDEX_in_COMPS( 0 )
  , DROP_IMPOSED( imposed_out )
  , NB_DOFs( 0 )
  , DOF_2_UNKNOWN( 0 )
  , DOF_IN_UNKNOWNS( 0 )
  , DIS_LINK( 0 )
{
   set_ordering_option( ordering  ) ;

   for ( size_t iDOFc=0 ; iDOFc<NB_COMPS ; iDOFc++ )
   {
      COMPS(iDOFc) = iDOFc ;
   }
   INDEX_in_COMPS = COMPS ;
   reset() ;
}

//----------------------------------------------------------------------
PDE_LinkDOF2Unknown*
PDE_LinkDOF2Unknown:: create( PEL_Object* a_owner,
                              PDE_DiscreteField const* ff,
                              size_t ic,
                              bool imposed_out )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LinkDOF2Unknown:: create" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;
   PEL_CHECK_PRE( ic < ff->nb_components() ) ;

   PDE_LinkDOF2Unknown* result =
       new PDE_LinkDOF2Unknown( a_owner, ff, ic, imposed_out ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->field() == ff ) ;
   PEL_CHECK_POST( result->nb_field_nodes() == ff->nb_nodes() ) ;
   PEL_CHECK_POST( result->components_table().size() == 1 ) ;
   PEL_CHECK_POST( result->components_table()(0) == ic ) ;
   PEL_CHECK_POST( result->DOFs_with_imposed_value_are_dropped() == 
                                                               imposed_out ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_LinkDOF2Unknown*
PDE_LinkDOF2Unknown:: create( PEL_Object* a_owner,
                              PDE_DiscreteField const* ff,
                              bool imposed_out )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LinkDOF2Unknown:: create" ) ;
   PEL_CHECK_PRE( ff != 0 ) ;
   PEL_CHECK_PRE( ff->nb_components() == 1 ) ;

   PDE_LinkDOF2Unknown* result =
         new PDE_LinkDOF2Unknown( a_owner, ff, 0, imposed_out ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->field() == ff ) ;
   PEL_CHECK_POST( result->nb_field_nodes() == ff->nb_nodes() ) ;
   PEL_CHECK_POST( result->components_table().size() == 1 ) ;
   PEL_CHECK_POST( result->components_table()(0) == 0 ) ;
   PEL_CHECK_POST( result->DOFs_with_imposed_value_are_dropped() == 
                                                               imposed_out ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_LinkDOF2Unknown:: PDE_LinkDOF2Unknown( PEL_Object* a_owner,
                                           PDE_DiscreteField const* ff,
                                           size_t ic,
                                           bool imposed_out )
//----------------------------------------------------------------------
  : PEL_Object( a_owner )
  , FIELD( ff )
  , NB_COMPS( 1 )
  , COMPS( 1 )
  , INDEX_in_COMPS( ic+1 )
  , DROP_IMPOSED( imposed_out )
  , ORDERING( sequenceOfTheNodes )
  , NB_DOFs( 0 )
  , DOF_2_UNKNOWN( 0 )
  , DOF_IN_UNKNOWNS( 0 )
  , DIS_LINK( 0 )
{
   COMPS( 0 ) = ic ;
   INDEX_in_COMPS( ic ) = 0 ;
   reset() ;
}

//----------------------------------------------------------------------
PDE_LinkDOF2Unknown:: ~PDE_LinkDOF2Unknown( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
PDE_LinkDOF2Unknown*
PDE_LinkDOF2Unknown:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LinkDOF2Unknown:: create_clone" ) ;

   PDE_LinkDOF2Unknown* result = new PDE_LinkDOF2Unknown( a_owner, this ) ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_LinkDOF2Unknown:: PDE_LinkDOF2Unknown( PEL_Object* a_owner,
                                           PDE_LinkDOF2Unknown const* other )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , FIELD( other->FIELD )
   , NB_NODES( other->NB_NODES )
   , NB_COMPS( other->NB_COMPS )
   , COMPS( other->COMPS )
   , INDEX_in_COMPS( other->INDEX_in_COMPS )
   , DROP_IMPOSED( other->DROP_IMPOSED )
   , ORDERING( other->ORDERING )
   , NB_DOFs( other->NB_DOFs )
   , DOF_2_UNKNOWN( other->DOF_2_UNKNOWN )
   , DOF_IN_UNKNOWNS( other->DOF_IN_UNKNOWNS )
{
   if( FIELD->is_distributed() )
   {
      DIS_LINK = PDE_CrossProcessUnknownNumbering::create( this, this ) ;
   }
}

//----------------------------------------------------------------------
void 
PDE_LinkDOF2Unknown:: reset( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LinkDOF2Unknown:: reset" ) ;

   NB_NODES = FIELD->nb_nodes() ;
   NB_DOFs = NB_NODES*NB_COMPS  ;

   //??? on pourrait utiliser NOT_DOF_IN_UNKNOWNS, ce qui éviterait
   //??? l'affectation à true
   DOF_IN_UNKNOWNS.re_initialize( NB_DOFs ) ;
   DOF_IN_UNKNOWNS.set( true ) ;

   DOF_2_UNKNOWN.re_initialize( NB_DOFs ) ;
   for( size_t n=0 ; n<NB_NODES ; n++ )
   {
      for( size_t iC=0 ; iC<NB_COMPS ; iC++ )
      {
         if( !FIELD->node_is_active( n ) ||
             ( DROP_IMPOSED && FIELD->DOF_has_imposed_value(n,COMPS(iC)) ) ||
             ( FIELD->DOF_is_constrained( n, COMPS(iC) ) ) )
         {
            DOF_IN_UNKNOWNS( local_index_of_DOF( n, iC ) ) = false ;
         }
      }
   }
   NB_DOFs = 0 ;
   for( size_t idx=0 ; idx<DOF_2_UNKNOWN.size() ; idx++ )
   {
      if( DOF_IN_UNKNOWNS( idx ) )
      {
         DOF_2_UNKNOWN(idx) = NB_DOFs++ ;
      }
   }
   if( FIELD->is_distributed() )
   {
      if( DIS_LINK==0 )
         DIS_LINK = PDE_CrossProcessUnknownNumbering::create( this, this ) ;
      else
         DIS_LINK->reset() ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( nb_field_nodes() == field()->nb_nodes() ) ;
}

//----------------------------------------------------------------------
void 
PDE_LinkDOF2Unknown:: reset( boolVector const& observed_nodes )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LinkDOF2Unknown:: reset" ) ;
   PEL_CHECK_PRE( observed_nodes.size() == field()->nb_nodes() ) ;

   NB_NODES = FIELD->nb_nodes() ;
   NB_DOFs = NB_NODES*NB_COMPS  ;

   //??? on pourrait utiliser NOT_DOF_IN_UNKNOWNS, ce qui éviterait
   //??? l'affectation à true
   DOF_IN_UNKNOWNS.re_initialize( NB_DOFs ) ;
   DOF_IN_UNKNOWNS.set( true ) ;

   DOF_2_UNKNOWN.re_initialize( NB_DOFs ) ;
   for( size_t n=0 ; n<NB_NODES ; n++ )
   {
      for( size_t iC=0 ; iC<NB_COMPS ; iC++ )
      {
         //??? seule difference par rapport à la fonction reset()
         if( !observed_nodes( n ) ||
             ( DROP_IMPOSED && FIELD->DOF_has_imposed_value(n,COMPS(iC)) ) ||
             ( FIELD->DOF_is_constrained( n, COMPS(iC) ) ) )
         {
            DOF_IN_UNKNOWNS( local_index_of_DOF( n, iC ) ) = false ;
         }
      }
   }
   NB_DOFs = 0 ;
   for( size_t idx=0 ; idx<DOF_2_UNKNOWN.size() ; idx++ )
   {
      if( DOF_IN_UNKNOWNS( idx ) )
      {
         DOF_2_UNKNOWN(idx) = NB_DOFs++ ;
      }
   }
   if( FIELD->is_distributed() )
   {
      if( DIS_LINK==0 )
         DIS_LINK = PDE_CrossProcessUnknownNumbering::create( this, this ) ;
      else
         DIS_LINK->reset() ;
   }

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( nb_field_nodes() == field()->nb_nodes() ) ;
}

//----------------------------------------------------------------------
PDE_DiscreteField const*
PDE_LinkDOF2Unknown:: field( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LinkDOF2Unknown:: field" ) ;
   PEL_CHECK_INV( invariant() ) ;

   return( FIELD ) ;
}

//----------------------------------------------------------------------
size_t_vector const&
PDE_LinkDOF2Unknown:: components_table( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LinkDOF2Unknown:: components_table" ) ;
   PEL_CHECK_INV( invariant() ) ;

   size_t_vector const& result = COMPS ;

   PEL_CHECK_POST( result.size() == 1 || 
                   result.size() == field()->nb_components() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LinkDOF2Unknown:: nb_field_nodes( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LinkDOF2Unknown:: nb_field_nodes" ) ;

   return( NB_NODES ) ;
}

//----------------------------------------------------------------------
bool
PDE_LinkDOF2Unknown:: DOFs_with_imposed_value_are_dropped( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LinkDOF2Unknown:: DOFs_with_imposed_value_are_dropped" ) ;

   return( DROP_IMPOSED ) ;
}

//----------------------------------------------------------------------
std::string
PDE_LinkDOF2Unknown:: DOFs_ordering_in_unknown( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LinkDOF2Unknown:: DOFs_ordering_in_unknown" ) ;

   std::string result ;
   if( ORDERING == sequenceOfTheComponents )
   {
      result = "sequence_of_the_components" ;
   }
   else if( ORDERING == sequenceOfTheNodes )
   {
      result = "sequence_of_the_nodes" ;
   }

   PEL_CHECK_POST( result == "sequence_of_the_components" || 
                   result == "sequence_of_the_nodes" )  ;

   return( result ) ;
}

//-----------------------------------------------------------------------
PDE_CrossProcessUnknownNumbering const*
PDE_LinkDOF2Unknown:: cross_process_numbering( void ) const 
//-----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LinkDOF2Unknown:: cross_process_numbering" ) ;
   PEL_CHECK_PRE( field()->is_distributed() ) ;

   PDE_CrossProcessUnknownNumbering const* result = DIS_LINK ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->link() == this ) ;   
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LinkDOF2Unknown:: unknown_vector_size( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LinkDOF2Unknown:: unknown_vector_size" ) ;

   size_t result = NB_DOFs ;

   PEL_CHECK_POST( result <= components_table().size()*nb_field_nodes() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_LinkDOF2Unknown:: DOF_is_unknown( size_t n, size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LinkDOF2Unknown:: DOF_is_unknown" ) ;
   PEL_CHECK_PRE( n < nb_field_nodes() ) ;
   PEL_CHECK_PRE( ic < field()->nb_components() ) ;
   PEL_CHECK_PRE( components_table().has( ic ) ) ;

   return( DOF_IN_UNKNOWNS( local_index_of_DOF( n, INDEX_in_COMPS( ic ) ) ) ) ;
}

//----------------------------------------------------------------------
size_t
PDE_LinkDOF2Unknown:: unknown_linked_to_DOF( size_t n, size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LinkDOF2Unknown:: unknown_linked_to_DOF" ) ;
   PEL_CHECK_PRE( n < nb_field_nodes() ) ;
   PEL_CHECK_PRE( ic < field()->nb_components() ) ;
   PEL_CHECK_PRE( components_table().has( ic ) ) ;
   PEL_CHECK_PRE( DOF_is_unknown( n, ic ) ) ;

   size_t result = 
          DOF_2_UNKNOWN( local_index_of_DOF( n, INDEX_in_COMPS( ic ) ) ) ;
   
   PEL_CHECK_POST( result < unknown_vector_size() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_LinkDOF2Unknown:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_LinkDOF2Unknown:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( DIS_LINK!=0 )
   {
      DIS_LINK->print( os, indent_width ) ;
   }
   else
   {
      std::string s( indent_width, ' ' ) ;
      os << s << "Number of unknowns for " << FIELD->name()
         << " is " <<  NB_DOFs << std::endl ;
   }   
}

//-----------------------------------------------------------------------
void
PDE_LinkDOF2Unknown:: set_ordering_option( std::string const& option ) 
//-----------------------------------------------------------------------
{
   if( option == "sequence_of_the_components" )
   {
      ORDERING = sequenceOfTheComponents ;
   }
   else if( option == "sequence_of_the_nodes" )
   {
      ORDERING = sequenceOfTheNodes ;
   }  
   else
   {
      PEL_Error::object()->raise_plain( option + ": bad value " ) ;
   }
}

//----------------------------------------------------------------------
bool
PDE_LinkDOF2Unknown:: invariant( void ) const
//----------------------------------------------------------------------
{
   PEL_ASSERT( FIELD != 0 ) ;
   PEL_ASSERT( NB_COMPS==1 || NB_COMPS==FIELD->nb_components() ) ;
   return( true ) ;
}


