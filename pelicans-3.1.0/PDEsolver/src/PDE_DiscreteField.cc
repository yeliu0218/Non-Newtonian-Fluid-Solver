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

#include <PDE_DiscreteField.hh>

#include <PEL.hh>
#include <PEL_Bool.hh>
#include <PEL_BoolVector.hh>
#include <PEL_Int.hh>
#include <PEL_IntArray2D.hh>
#include <PEL_Communicator.hh>
#include <PEL_Data.hh>
#include <PEL_DoubleArray3D.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Module.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectReader.hh>
#include <PEL_ObjectWriter.hh>
#include <PEL_String.hh>
#include <PEL_Vector.hh>
#include <stringVector.hh>

#include <LA_SeqVector.hh>

#include <PDE_DOFconstraintsIterator.hh>
#include <PDE_CrossProcessNodeNumbering.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_ReferenceElement.hh>

#include <iomanip>
#include <iostream>
#include <sstream>

using std::cout ; using std::endl ;
using std::string ; using std::ostringstream ;
using std::vector ;

size_t
PDE_DiscreteField:: NB_INSTANCES = 0 ;

//----------------------------------------------------------------------
PDE_DiscreteField*
PDE_DiscreteField:: create( PEL_Object* a_owner,
                            std::string const& a_name,
                            size_t a_nb_components,
                            size_t a_depth )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: create" ) ;
   PEL_CHECK_PRE( !a_name.empty() ) ;
   PEL_CHECK_PRE( a_nb_components>0 ) ;
   
   PDE_DiscreteField* result = new PDE_DiscreteField( a_owner,
                                                      a_name,
                                                      a_nb_components,
                                                      a_depth ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   PEL_CHECK_POST( result->storage_depth() == a_depth ) ;
   PEL_CHECK_POST( result->nb_nodes() == 0 ) ;
   PEL_CHECK_POST( result->nb_components() == a_nb_components ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_DiscreteField:: PDE_DiscreteField( PEL_Object* a_owner,
                                       std::string const& a_name,
                                       size_t a_nb_components,
                                       size_t a_depth )
//----------------------------------------------------------------------
    : PEL_Object( a_owner )
    , FNAME( a_name )
    , ID( NB_INSTANCES++ )
    , NB_COMPS( a_nb_components )
    , NB_NODES( 0 )
    , STO_DEPTH( a_depth )
    , VALUES( 0, 0, 0, 0. )
    , DBC_IDX( 0, 0 )
    , DBC_FLAG( 0, 0 )
    , DBC_VALUES( 0 )
    , ACTIVE_NODES( 0 )
    , DF( 0 )
    , CSTR( 0 )
    , CSTR_IT( 0 )
{
}

//----------------------------------------------------------------------
PDE_DiscreteField*
PDE_DiscreteField:: create_duplication( PEL_Object* a_owner,
                                        std::string const& a_name ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: create_duplication" ) ;
   PEL_CHECK_PRE( !a_name.empty()  ) ;

   //??? here constraints are ignored and should not be
   PEL_ASSERT( CSTR == 0 ) ;
   
   PDE_DiscreteField* result = new PDE_DiscreteField( a_owner, this, a_name ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   PEL_CHECK_POST( result->storage_depth() == storage_depth() ) ;
   PEL_CHECK_POST( result->nb_nodes() == nb_nodes() ) ;
   PEL_CHECK_POST( result->nb_components() == nb_components() ) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t s=0 ; s<result->storage_depth() ; ++s ),
         FORALL( ( size_t n=0 ; n<result->nb_nodes() ; ++n ),
       FORALL( ( size_t ic=0 ; ic<result->nb_components() ; ++ic ),
           result->DOF_value( s, n, ic ) == DOF_value( s, n, ic ) ) ) ) ) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t n=0 ; n<result->nb_nodes() ; ++n ),
         FORALL( ( size_t ic=0 ; ic<result->nb_components() ; ++ic ),
       result->DOF_has_imposed_value( n, ic ) == DOF_has_imposed_value( n, ic ) ))) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t n=0 ; n<result->nb_nodes() ; ++n ),
         FORALL( ( size_t ic=0 ; ic<result->nb_components() ; ++ic ),
       IMPLIES( result->DOF_has_imposed_value( n, ic ),
                     result->DOF_imposed_value(n,ic) == DOF_imposed_value(n,ic) )))) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t n=0 ; n<result->nb_nodes() ; ++n ),
         EQUIVALENT( result->node_is_active( n ), node_is_active( n ) ) ) ) ;
   PEL_CHECK_POST( !result->is_distributed() ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
PDE_DiscreteField*
PDE_DiscreteField:: create_duplication( PEL_Object* a_owner,
                                        std::string const& a_name,
                                        size_t a_nb_components,
                                        size_t a_depth ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: create_duplication" ) ;
   PEL_CHECK_PRE( !a_name.empty()  ) ;

   //??? here constraints are ignored and should not be
   PEL_ASSERT( CSTR == 0 ) ;
   
   PDE_DiscreteField* result = new PDE_DiscreteField( a_owner, 
                                                      this, 
                                                      a_name,
                                                      a_nb_components,
                                                      a_depth ) ;

   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   PEL_CHECK_POST( result->storage_depth() == a_depth ) ;
   PEL_CHECK_POST( result->nb_nodes() == nb_nodes() ) ;
   PEL_CHECK_POST( result->nb_components() == a_nb_components ) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t s=0 ; s<result->storage_depth() ; ++s ),
         FORALL( ( size_t n=0 ; n<result->nb_nodes() ; ++n ),
       FORALL( ( size_t ic=0 ; ic<result->nb_components() ; ++ic ),
           result->DOF_value( s, n, ic ) == DOF_value( s, n, ic ) ) ) ) ) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t n=0 ; n<result->nb_nodes() ; ++n ),
         FORALL( ( size_t ic=0 ; ic<result->nb_components() ; ++ic ),
       result->DOF_has_imposed_value( n, ic ) == DOF_has_imposed_value( n, ic ) ))) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t n=0 ; n<result->nb_nodes() ; ++n ),
         FORALL( ( size_t ic=0 ; ic<result->nb_components() ; ++ic ),
       IMPLIES( result->DOF_has_imposed_value( n, ic ),
                     result->DOF_imposed_value(n,ic) == DOF_imposed_value(n,ic) )))) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t n=0 ; n<result->nb_nodes() ; ++n ),
         EQUIVALENT( result->node_is_active( n ), node_is_active( n ) ) ) ) ;
   PEL_CHECK_POST( !result->is_distributed() ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
PDE_DiscreteField:: PDE_DiscreteField( PEL_Object* a_owner,
                                       PDE_DiscreteField const* other,
                                       std::string const& a_name )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , FNAME( a_name )
   , ID( NB_INSTANCES++ )
   , NB_COMPS( other->NB_COMPS )
   , NB_NODES( other->NB_NODES )
   , STO_DEPTH( other->STO_DEPTH )
   , VALUES( other->VALUES )
   , DBC_IDX( other->DBC_IDX )
   , DBC_FLAG( other->DBC_FLAG )
   , DBC_VALUES( other->DBC_VALUES )
   , ACTIVE_NODES( other->ACTIVE_NODES )
   , DF( 0 )
   , CSTR( 0 )
   , CSTR_IT( 0 )
{
}

//----------------------------------------------------------------------
PDE_DiscreteField:: PDE_DiscreteField( PEL_Object* a_owner,
                                       PDE_DiscreteField const* other,
                                       std::string const& a_name,
                                       size_t a_nb_components,
                                       size_t a_depth )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , FNAME( a_name )
   , ID( NB_INSTANCES++ )
   , NB_COMPS( a_nb_components )
   , NB_NODES( other->NB_NODES )
   , STO_DEPTH( a_depth )
   , VALUES( other->VALUES )
   , DBC_IDX( other->DBC_IDX )
   , DBC_FLAG( other->DBC_FLAG )
   , DBC_VALUES( other->DBC_VALUES )
   , ACTIVE_NODES( other->ACTIVE_NODES )
   , DF( 0 )
   , CSTR( 0 )
   , CSTR_IT( 0 )
{
}

//----------------------------------------------------------------------
PDE_DiscreteField:: ~PDE_DiscreteField( void ) 
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
size_t
PDE_DiscreteField:: nb_objects( void )
//----------------------------------------------------------------------
{
   return( NB_INSTANCES ) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: set_nb_nodes( size_t a_nb_nodes )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: set_nb_nodes" ) ;
   PEL_CHECK_PRE( nb_nodes() == 0 ) ;
   
   NB_NODES = a_nb_nodes ;
   VALUES.re_initialize( a_nb_nodes, NB_COMPS, STO_DEPTH, 0. ) ;
   
   PEL_CHECK_POST( nb_nodes() == a_nb_nodes ) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t s=0 ; s<storage_depth() ; ++s ),
         FORALL( ( size_t n=0 ; n<nb_nodes() ; ++n ),
            FORALL( ( size_t ic=0 ; ic<nb_components() ; ++ic ),
               DOF_value( s, n, ic ) == 0. ) ) ) ) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t n=0 ; n<nb_nodes() ; ++n ),
         FORALL( ( size_t ic=0 ; ic<nb_components() ; ++ic ),
            DOF_has_imposed_value( n, ic ) == false ))) ;
   PEL_CHECK_POST( FORALL( ( size_t n=0 ; n<nb_nodes() ; ++n ),
                      node_is_active( n ) ) ) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: add_nodes( size_t nb_supp_nodes )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: add_nodes" ) ;
   PEL_CHECK_PRE( nb_supp_nodes != 0 ) ;
   PEL_SAVEOLD( size_t, nb_nodes, nb_nodes() ) ;
   
   if( ACTIVE_NODES.size() == 0 )
   {
      ACTIVE_NODES.re_initialize( NB_NODES, true ) ;
   }
   NB_NODES += nb_supp_nodes ;  
   VALUES.raise_first_index_bound( NB_NODES, 0. ) ;
   ACTIVE_NODES.resize( NB_NODES, false ) ;
   if( DBC_IDX.index_bound(0) != 0 )
   {
      DBC_IDX.raise_first_index_bound( NB_NODES, PEL::bad_index() ) ;
   }
   
   if( CSTR != 0 ) CSTR->raise_nb_nodes( NB_NODES ) ;
   
  
   PEL_CHECK_POST( nb_nodes() == OLD( nb_nodes ) + nb_supp_nodes ) ;
   PEL_CHECK_POST( FORALL( ( size_t n=OLD( nb_nodes ) ; n<nb_nodes() ; ++n ),
                      !node_is_active( n ) ) ) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t s=0 ; s<storage_depth() ; ++s ),
         FORALL( ( size_t n=OLD( nb_nodes ) ; n<nb_nodes() ; ++n ),
	    FORALL( ( size_t ic=0 ; ic<nb_components() ; ++ic ),
	        DOF_value( s, n, ic ) == 0.0 ) ) ) ) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t n=OLD( nb_nodes ) ; n<nb_nodes() ; ++n ),
         FORALL( ( size_t ic=0 ; ic<nb_components() ; ++ic ),
	    DOF_has_imposed_value( n, ic ) == false ))) ;
}

//----------------------------------------------------------------------
size_t
PDE_DiscreteField:: id_number( void ) const
//----------------------------------------------------------------------
{
   size_t result = ID ;

   PEL_CHECK_POST( result < nb_objects() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
string const& 
PDE_DiscreteField:: name( void ) const
//----------------------------------------------------------------------
{
   return( FNAME ) ;
}

//----------------------------------------------------------------------
bool
PDE_DiscreteField:: is_distributed( void ) const
//----------------------------------------------------------------------
{
   return( DF!=0 ) ;   
}

//----------------------------------------------------------------------
PDE_CrossProcessNodeNumbering* 
PDE_DiscreteField:: cross_process_numbering( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: cross_process_numbering" ) ;
   PEL_CHECK_PRE( is_distributed() ) ;

   PDE_CrossProcessNodeNumbering* result = DF ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==this ) ;
   PEL_CHECK_POST( result->field()==this ) ;
   return( result ) ;   
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: build_cross_process_globalization( 
                                PDE_CrossProcessNodeNumbering* numbering ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: build_cross_process_globalization" ) ;
   PEL_CHECK_PRE( !is_distributed() ) ;
   PEL_CHECK_PRE( numbering != 0 ) ;
   PEL_CHECK_PRE( numbering->owner() == this ) ;

   DF = numbering ;
   DF->attach_field( this ) ;

   PEL_CHECK_POST( is_distributed() ) ;
   PEL_CHECK_POST( cross_process_numbering() == numbering ) ;
}

//----------------------------------------------------------------------
size_t
PDE_DiscreteField:: nb_nodes( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: nb_nodes" ) ;
   return( NB_NODES ) ;
}

//----------------------------------------------------------------------
size_t
PDE_DiscreteField:: nb_components( void ) const
//----------------------------------------------------------------------
{
   return( NB_COMPS ) ;
}

//----------------------------------------------------------------------
bool
PDE_DiscreteField:: node_is_active( size_t n ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: node_is_active" ) ;
   PEL_CHECK_PRE( n < nb_nodes() ) ;

   return( ACTIVE_NODES.size()==0 || ACTIVE_NODES( n ) ) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: set_node_active( size_t n )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: set_node_active" ) ;
   PEL_CHECK_PRE( n < nb_nodes() ) ;

   if( ACTIVE_NODES.size() != 0 ) ACTIVE_NODES(n) = true ;

   PEL_CHECK_POST( node_is_active( n ) ) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: set_node_inactive( size_t n )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: set_node_inactive" ) ; 
   PEL_CHECK_PRE( n < nb_nodes() ) ;
   
   if( ACTIVE_NODES.size() == 0 )
   {
      ACTIVE_NODES.re_initialize( NB_NODES, true ) ;
   }
   ACTIVE_NODES(n) = false ;

   PEL_CHECK_POST( !node_is_active( n ) ) ;
}

//----------------------------------------------------------------------
bool
PDE_DiscreteField:: DOF_is_fixed( size_t n, size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: DOF_is_fixed" ) ;
   PEL_CHECK_PRE( n < nb_nodes() ) ;
   PEL_CHECK_PRE( ic < nb_components() ) ;

   bool result = !node_is_active( n ) || DOF_has_imposed_value(n,ic) ;

   PEL_CHECK_POST( 
      EQUIVALENT( result,
                  !node_is_active( n ) || DOF_has_imposed_value(n,ic) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PDE_DiscreteField:: DOF_is_free( size_t n, size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: DOF_is_free" ) ;
   PEL_CHECK_PRE( n < nb_nodes() ) ;
   PEL_CHECK_PRE( ic < nb_components() ) ;

   bool result = !DOF_is_fixed( n, ic) ;

   PEL_CHECK_POST( EQUIVALENT( result, !DOF_is_fixed(n,ic) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
size_t
PDE_DiscreteField:: storage_depth( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: storage_depth" ) ;

   return( STO_DEPTH ) ;
}

//----------------------------------------------------------------------
bool
PDE_DiscreteField:: DOF_has_imposed_value( size_t n, size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: DOF_has_imposed_value" ) ;
   PEL_CHECK_PRE( n < nb_nodes() ) ;
   PEL_CHECK_PRE( ic< nb_components() ) ;

   return( DBC_IDX.index_bound( 0 ) != 0 &&
           DBC_IDX( n, ic ) != PEL::bad_index() &&
           DBC_FLAG( DBC_IDX( n, ic ) ) ) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: set( PDE_DiscreteField const* other )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: set" ) ;
   PEL_CHECK_PRE( other->nb_nodes() == nb_nodes() ) ;
   PEL_CHECK_PRE( other->nb_components() == nb_components() ) ;
   PEL_CHECK_PRE( other->storage_depth() == storage_depth() ) ;

   //??? here constraints are ignored and should not be
   PEL_ASSERT( CSTR == 0 ) ;
   
   VALUES = other->VALUES ;
   DBC_IDX = other->DBC_IDX ;
   DBC_FLAG = other->DBC_FLAG ;
   DBC_VALUES = other->DBC_VALUES ;
   ACTIVE_NODES = other->ACTIVE_NODES ;

   if( is_distributed() )
   {
      DF->synchronize_imposed_values() ;
      DF->synchronize_valid_nodes() ;
   }
      
   PEL_CHECK_POST( 
      FORALL( ( size_t s=0 ; s<storage_depth() ; ++s ),
         FORALL( ( size_t n=0 ; n<nb_nodes() ; ++n ),
	    FORALL( ( size_t ic=0 ; ic<nb_components() ; ++ic ),
	        DOF_value( s, n, ic ) == other->DOF_value( s, n, ic ) ) ) ) ) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t n=0 ; n<nb_nodes() ; ++n ),
         FORALL( ( size_t ic=0 ; ic<nb_components() ; ++ic ),
	    DOF_has_imposed_value( n, ic ) == other->DOF_has_imposed_value( n, ic ) ))) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t n=0 ; n<nb_nodes() ; ++n ),
         FORALL( ( size_t ic=0 ; ic<nb_components() ; ++ic ),
	    IMPLIES( DOF_has_imposed_value( n, ic ),
                     DOF_imposed_value(n,ic) == other->DOF_imposed_value(n,ic) )))) ;
   PEL_CHECK_POST( 
      FORALL( ( size_t n=0 ; n<nb_nodes() ; ++n ),
         EQUIVALENT( node_is_active( n ), other->node_is_active( n ) ) ) ) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: modify_DOF( bool imposed,
                                double value, size_t n, size_t ic )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: modify_DOF" ) ;
   PEL_CHECK_PRE( n < nb_nodes() ) ;
   PEL_CHECK_PRE( ic < nb_components() ) ;

   if( imposed )
   {
      equip_DOF_with_imposed_value( n, value, ic ) ;
   }
   else if( DOF_has_imposed_value( n, ic ) )
   {
      DBC_FLAG( DBC_IDX( n, ic ) ) = false ;
   }
      
   for( size_t level=0 ; level<STO_DEPTH ; level++ )
   {
      VALUES( n, ic, level ) = value ;
   }

   PEL_CHECK_POST( DOF_has_imposed_value( n, ic ) == imposed ) ;
   PEL_CHECK_POST( IMPLIES( imposed, DOF_imposed_value( n, ic ) == value ) ) ;
   PEL_CHECK_POST( FORALL( ( size_t s=0 ; s<storage_depth() ; s++),
                           DOF_value( s, n, ic ) == value ) ) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: equip_DOF_with_imposed_value(
                                           size_t n, double x, size_t ic ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: equip_DOF_with_imposed_value" ) ;
   PEL_CHECK_PRE( n < nb_nodes()  ) ;
   PEL_CHECK_PRE( ic < nb_components() ) ;

   if( DBC_IDX.index_bound(0) == 0 )
   {
      DBC_IDX.re_initialize( NB_NODES, NB_COMPS, PEL::bad_index() ) ;
   }
   if( DBC_IDX( n, ic ) == PEL::bad_index() )
   {
      size_t const idx = DBC_VALUES.size() ;
      DBC_IDX( n, ic ) = idx ;
      DBC_FLAG.resize( idx+1, true ) ;
      DBC_VALUES.resize( idx+1, x ) ;
   }
   else
   {
      size_t const idx = DBC_IDX( n, ic ) ;
      DBC_FLAG( idx ) = true ;
      DBC_VALUES( idx ) = x ;
   }

   PEL_CHECK_POST( DOF_has_imposed_value( n, ic ) ) ;
   PEL_CHECK_POST( DOF_imposed_value( n, ic ) == x ) ;
}

//----------------------------------------------------------------------
double
PDE_DiscreteField:: DOF_imposed_value( size_t n, size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: DOF_imposed_value" ) ;
   PEL_CHECK_PRE( n < nb_nodes() ) ;
   PEL_CHECK_PRE( ic < nb_components() ) ;
   PEL_CHECK_PRE( DOF_has_imposed_value( n, ic ) ) ;

   return( DBC_VALUES( DBC_IDX( n, ic ) ) ) ;
}

//----------------------------------------------------------------------
double
PDE_DiscreteField:: DOF_value( size_t level, size_t n, size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: DOF_value" ) ;
   PEL_CHECK_PRE( n < nb_nodes() ) ;
   PEL_CHECK_PRE( ic < nb_components() ) ;
   PEL_CHECK_PRE( level < storage_depth() ) ;
   
   return( VALUES( n, ic, level ) ) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: extract_DOFs_value( size_t level,
                                        LA_SeqVector* vec,
                                        size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: extract_DOFs_value" ) ;
   PEL_CHECK_PRE( level < storage_depth() ) ;
   PEL_CHECK_PRE( ic < nb_components() ) ;
   PEL_CHECK_PRE( vec->nb_rows() == nb_nodes() ) ;
   
   for( size_t nn=0 ; nn<NB_NODES ; nn++ )
   {
      vec->set_item( nn, VALUES( nn, ic, level ) ) ;
   }

   vec->synchronize() ;

   PEL_CHECK_POST( vec->is_synchronized() ) ;
   PEL_CHECK_POST( 
      FORALL ( ( size_t i=0 ; i<vec->nb_rows() ; ++i ),
               vec->item(i) == DOF_value( level, i, ic ) ) ) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: extract_unknown_DOFs_value(
                                       size_t level,
                                       LA_SeqVector* vec,
                                       PDE_LinkDOF2Unknown const* link ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: extract_unknown_DOFs_value" ) ;
   PEL_CHECK_PRE( level < storage_depth() ) ;
   PEL_CHECK_PRE( link != 0 ) ;
   PEL_CHECK_PRE( link->field() == this ) ;
   PEL_CHECK_PRE( vec != 0 ) ;
   PEL_CHECK_PRE( vec->nb_rows() >= link->unknown_vector_size() ) ;
   
   size_t_vector const& fComp = link->components_table() ;

   vec->nullify() ;
   for( size_t nn=0 ; nn<nb_nodes() ; nn++ )
   {
      for( size_t j=0 ; j<fComp.size() ; j++ )
      {
         size_t ic = fComp( j ) ;
         if( link->DOF_is_unknown( nn, ic ) )
         {
            size_t ivec = link->unknown_linked_to_DOF( nn, ic ) ;
            vec->set_item( ivec, VALUES( nn, ic, level ) ) ;
         }
      }
   }
   vec->synchronize() ;

   PEL_CHECK_POST( vec->is_synchronized() ) ;
   PEL_CHECK_POST( 
      FORALL(
         ( size_t n=0 ; n<nb_nodes() ; ++n ), 
         FORALL(
            ( size_t j=0 ; j<link->components_table().size() ; ++j ),
            IMPLIES(
               link->DOF_is_unknown(n,link->components_table()(j)),
               vec->item( link->unknown_linked_to_DOF(
                                    n,link->components_table()(j)))== 
                DOF_value(level,n,link->components_table()(j)))))) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: set_DOF_value( size_t level, size_t n, double x, size_t ic )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: set_DOF_value" ) ;
   PEL_CHECK_PRE( level < storage_depth() ) ;
   PEL_CHECK_PRE( n < nb_nodes()  ) ;
   PEL_CHECK_PRE( ic < nb_components() ) ;
   
   VALUES( n, ic, level ) = x ;

   PEL_CHECK_POST( DOF_value(level,n,ic) == x ) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: set_DOF_imposed_value( size_t n, double x, size_t ic )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: set_DOF_imposed_value" ) ;
   PEL_CHECK_PRE( n < nb_nodes()  ) ;
   PEL_CHECK_PRE( ic < nb_components() ) ;
   PEL_CHECK_PRE( DOF_has_imposed_value( n, ic ) ) ;

   DBC_VALUES( DBC_IDX( n, ic ) ) = x ;

   PEL_CHECK_POST( DOF_imposed_value( n, ic ) == x ) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: update_DOFs_value( size_t level,
                                       LA_SeqVector const* vec,
                                       PDE_LinkDOF2Unknown const* link )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: update_DOFs_value" ) ;
   PEL_CHECK_PRE( level<storage_depth() ) ;
   PEL_CHECK_PRE( link != 0 ) ;
   PEL_CHECK_PRE( link->field() == this ) ;
   PEL_CHECK_PRE( vec != 0 ) ;
   PEL_CHECK_PRE( vec->nb_rows() >= link->unknown_vector_size() ) ;

   size_t_vector const& fComp = link->components_table() ;

   for( size_t nn=0 ; nn<nb_nodes() ; nn++ )
   {
      for( size_t j=0 ; j<fComp.size() ; j++ )
      {
         size_t ic = fComp( j ) ;
         if( link->DOF_is_unknown( nn, ic ) )
         {
            size_t ivec = link->unknown_linked_to_DOF( nn, ic ) ;
            double val = vec->item( ivec ) ;
            set_DOF_value( level, nn, val, ic ) ;
         }
      }
   }
   PEL_CHECK_POST(
      FORALL(
         ( size_t n=0 ; n<nb_nodes() ; ++n ),
         FORALL(
            ( size_t j=0 ; j<link->components_table().size() ; ++j ),
            IMPLIES(
               link->DOF_is_unknown(n,link->components_table()(j)),
               DOF_value(level,n,link->components_table()(j))==
                   vec->item(link->unknown_linked_to_DOF(
                                    n,link->components_table()(j))))))) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: update_free_DOFs_value( size_t level,
                                            LA_SeqVector const* vec,
                                            PDE_LinkDOF2Unknown const* link )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: update_free_DOFs_value" ) ;
   PEL_CHECK_PRE( level < storage_depth() ) ;
   PEL_CHECK_PRE( link != 0 ) ;
   PEL_CHECK_PRE( link->field() == this ) ;
   PEL_CHECK_PRE( vec != 0 ) ;
   PEL_CHECK_PRE( vec->nb_rows() >= link->unknown_vector_size() ) ;
   
   size_t_vector const& fComp = link->components_table() ;

   for( size_t nn=0 ; nn<nb_nodes() ; nn++ )
   {
      for( size_t j=0 ; j<fComp.size() ; j++ )
      {
         size_t ic = fComp( j ) ;
         if( link->DOF_is_unknown( nn, ic ) )
         {
            if( !DOF_is_fixed( nn, ic ) )
            {
               size_t ivec = link->unknown_linked_to_DOF( nn, ic ) ;
               double val = vec->item( ivec ) ;
               set_DOF_value( level, nn, val, ic ) ;
            }
         }
      }
   }
   PEL_CHECK_POST(
      FORALL(
         ( size_t n=0 ; n<nb_nodes() ; ++n ),
         FORALL(
            ( size_t j=0 ; j<link->components_table().size() ; ++j ),
            IMPLIES(
               link->DOF_is_unknown(n,link->components_table()(j)) &&
                         DOF_is_free(n,link->components_table()(j)),
               DOF_value(level,n,link->components_table()(j)) ==
                     vec->item(link->unknown_linked_to_DOF(
                                     n,link->components_table()(j))))))) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: add_to_free_DOFs_value( size_t level,
                                            LA_SeqVector const* vec,
                                            PDE_LinkDOF2Unknown const* link,
                                            double alpha )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: add_to_free_DOFs_value" ) ;
   PEL_CHECK_PRE( level < storage_depth() ) ;
   PEL_CHECK_PRE( link != 0 ) ;
   PEL_CHECK_PRE( link->field() == this ) ;
   PEL_CHECK_PRE( vec != 0 ) ;
   PEL_CHECK_PRE( vec->nb_rows() >= link->unknown_vector_size() ) ;
   
   size_t_vector const& fComp = link->components_table() ;

   for( size_t nn=0 ; nn<nb_nodes() ; nn++ )
   {
      for( size_t j=0 ; j<fComp.size() ; j++ )
      {
         size_t ic = fComp( j ) ;
         if( link->DOF_is_unknown( nn, ic ) )
         {
            if( !DOF_is_fixed( nn, ic ) )
            {
               size_t ivec = link->unknown_linked_to_DOF( nn, ic ) ;
               double val =  DOF_value( level, nn, ic ) + 
                             alpha * vec->item( ivec );
               set_DOF_value( level, nn, val, ic ) ;
            }
         }
      }
   }
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: enforce_imposed_values_to_DOFs( size_t level )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: enforce_imposed_values_to_DOFs" ) ;
   PEL_CHECK_PRE( level<storage_depth() ) ;

   if( DBC_IDX.index_bound( 0 ) != 0 )
   {
      for( size_t n=0 ; n<NB_NODES ; ++n )
      {
         for( size_t ic=0 ; ic<NB_COMPS ; ++ic )
         {
            size_t idx = DBC_IDX( n, ic ) ;
            if( idx!=PEL::bad_index() && DBC_FLAG( idx ) )
            {
               VALUES( n, ic, level ) = DBC_VALUES( idx ) ;
            }
         }
      }
   }

   PEL_CHECK_POST(      
      FORALL( ( size_t n=0 ; n<nb_nodes() ; ++n ),
         FORALL( ( size_t ic=0 ; ic<nb_components() ; ++ic ),
	    IMPLIES( DOF_has_imposed_value(n,ic),
		     DOF_value(level,n,ic) == DOF_imposed_value(n,ic) )))) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: copy_DOFs_value( size_t source_level, size_t target_level )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: copy_DOFs_value" ) ;
   PEL_CHECK_PRE( source_level<storage_depth() ) ;
   PEL_CHECK_PRE( target_level<storage_depth() ) ;

   for( size_t n=0 ; n<NB_NODES ; n++ )
   {
      for( size_t ic=0 ; ic<NB_COMPS ; ic++ )
      {
         VALUES( n, ic, target_level ) = VALUES( n, ic, source_level ) ;
      }
   }

   PEL_CHECK_POST(
      FORALL( ( size_t n=0 ; n<nb_nodes() ; ++n ),
         FORALL( ( size_t ic=0 ; ic<nb_components() ; ++ic ),
	    DOF_value(target_level,n,ic) == DOF_value(source_level,n,ic) ))) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: copy_unknown_DOFs_value( size_t source_level,
                                             size_t target_level,
                                             PDE_LinkDOF2Unknown const* link )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: copy_unknown_DOFs_value" ) ;
   PEL_CHECK_PRE( source_level < storage_depth() ) ;
   PEL_CHECK_PRE( target_level < storage_depth() ) ;
   PEL_CHECK_PRE( link->field() == this ) ;

   size_t_vector const& fComp = link->components_table() ;

   for( size_t nn=0 ; nn<nb_nodes() ; nn++ )
   {
      for( size_t j=0 ; j<fComp.size() ; j++ )
      {
         size_t ic = fComp( j ) ;
         if( link->DOF_is_unknown( nn, ic ) )
         {
            VALUES( nn, ic, target_level ) = VALUES( nn, ic, source_level ) ;
         }
      }
   }
   PEL_CHECK_POST(
      FORALL(
         ( size_t n=0 ; n<nb_nodes() ; ++n ),
         FORALL(
            ( size_t j=0 ; j<link->components_table().size() ; ++j ),
            IMPLIES(
               link->DOF_is_unknown(n,link->components_table()(j)),
               DOF_value(target_level,n,link->components_table()(j))==
                  DOF_value(source_level,n,link->components_table()(j)))))) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: copy_unknown_free_DOFs_value( size_t source_level,
                                                  size_t target_level,
                                                  PDE_LinkDOF2Unknown const* link )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: copy_unknown_free_DOFs_value" ) ;
   PEL_CHECK_PRE( source_level < storage_depth() ) ;
   PEL_CHECK_PRE( target_level < storage_depth() ) ;
   PEL_CHECK_PRE( link->field() == this ) ;

   size_t_vector const& fComp = link->components_table() ;

   for( size_t nn=0 ; nn<nb_nodes() ; nn++ )
   {
      for( size_t j=0 ; j<fComp.size() ; j++ )
      {
         size_t ic = fComp( j ) ;
         if( link->DOF_is_unknown( nn, ic ) && !DOF_is_fixed( nn, ic ) )
         {
            VALUES( nn, ic, target_level ) = VALUES( nn, ic, source_level ) ;
         }
      }
   }
   PEL_CHECK_POST(
      FORALL(
         ( size_t n=0 ; n<nb_nodes() ; ++n ),
         FORALL(
            ( size_t j=0 ; j<link->components_table().size() ; ++j ),
            IMPLIES(
               link->DOF_is_unknown(n,link->components_table()(j)) &&
                                  DOF_is_free(n,link->components_table()(j)),
               DOF_value(target_level,n,link->components_table()(j))==
                    DOF_value(source_level,n,link->components_table()(j)))))) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: remove_constraint_for_DOF( size_t slave_n, size_t slave_ic )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: remove_constraint_for_DOF" ) ;
   PEL_CHECK_PRE( slave_n < nb_nodes() ) ;
   PEL_CHECK_PRE( slave_ic < nb_components() ) ;
   PEL_CHECK_PRE( DOF_is_constrained( slave_n, slave_ic ) ) ;
   
   CSTR->remove( slave_n, slave_ic ) ;
   
   PEL_CHECK_POST( !DOF_is_constrained( slave_n, slave_ic ) ) ; 
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: add_constraint_for_DOF( size_t slave_n, size_t slave_ic,
                                            size_t master_n, size_t master_ic, 
                                            double coef )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: add_constraint_for_DOF" ) ;
   PEL_CHECK_PRE( slave_n   < nb_nodes() ) ;
   PEL_CHECK_PRE( slave_ic  < nb_components() ) ;
   PEL_CHECK_PRE( master_n  < nb_nodes() ) ;
   PEL_CHECK_PRE( master_ic < nb_components() ) ;
   PEL_CHECK_PRE( ! DOF_is_constrained( master_n, master_ic ) ) ;
   PEL_CHECK_PRE( ! DOF_has_imposed_value( slave_n, slave_ic ) ) ;
   
   if( CSTR == 0 )
   {
      CSTR = new PDE_DOFconstraints( this, NB_NODES, NB_COMPS ) ;
   }
   
   CSTR->add( slave_n, slave_ic, master_n, master_ic, coef ) ;
   
   PEL_CHECK_POST( DOF_is_constrained( slave_n, slave_ic ) ) ; 
}

//----------------------------------------------------------------------
bool
PDE_DiscreteField:: DOF_is_constrained( size_t n, size_t ic ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: DOF_is_constrained" ) ;
   PEL_CHECK_PRE( n < nb_nodes() ) ;
   PEL_CHECK_PRE( ic < nb_components() ) ;
   
   bool result = ( ( CSTR != 0 ) && ( CSTR->has( n, ic ) ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_DOFconstraintsIterator*
PDE_DiscreteField:: create_constraints_iterator( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: create_constraints_iterator" ) ;

   PDE_DOFconstraintsIterator* result = 
                     new PDE_DOFconstraintsIterator( a_owner, this, CSTR ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->field() == this ) ;
   PEL_CHECK_POST( !result->is_valid() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: enforce_constraints_for_DOFs( size_t level )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: enforce_constraints_for_DOFs" ) ;
   PEL_CHECK_PRE( level < storage_depth() ) ;
   
   if( CSTR != 0 )
   {
      if( CSTR_IT == 0 ) CSTR_IT = create_constraints_iterator( this ) ;
      for( size_t n=0 ; n<NB_NODES ; ++n )
      {
         for( size_t ic=0 ; ic<NB_COMPS ; ++ic )
         {
            if( DOF_is_constrained( n, ic ) )
            {
               CSTR_IT->start( n, ic ) ;
               double val = 0.0 ;
               for( ; CSTR_IT->is_valid() ; CSTR_IT->go_next() )
               {
                  size_t cstr_nn = CSTR_IT->node_of_constraining_DOF() ;
                  size_t cstr_ic = CSTR_IT->component_of_constraining_DOF() ;
                  double cstr_xx = CSTR_IT->constraint_coefficient() ;
                  val += cstr_xx * VALUES( cstr_nn, cstr_ic, level ) ;
               }
               VALUES( n, ic, level ) = val ;
            }
         }
      }
   }
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: save_state( PEL_ObjectWriter* writer ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: save_state" ) ;
   PEL_CHECK_PRE( save_state_PRE( writer ) ) ;

   //??? here constraints are ignored and should not be
   PEL_ASSERT( CSTR == 0 ) ;

   writer->start_new_object( "PDE_DiscreteField" ) ;
   
   writer->add_entry( "name", PEL_String::create( 0, FNAME ) ) ;
   writer->add_entry( "nb_nodes", PEL_Int::create( 0, NB_NODES ) ) ;
   writer->add_entry( "nb_components", PEL_Int::create( 0, NB_COMPS ) ) ;
   writer->add_entry( "nb_levels", PEL_Int::create( 0, STO_DEPTH ) ) ;
   
   // Saving values
   {
      doubleVector vec( NB_NODES*NB_COMPS*STO_DEPTH ) ;
      size_t idx = 0 ;
      for( size_t n = 0 ; n<NB_NODES ; n++ )
      {
         for( size_t ic=0 ; ic<NB_COMPS ; ic++ )
         {
            for( size_t level=0 ; level<STO_DEPTH ; level++ )
            {
               vec( idx++ ) = VALUES( n, ic, level ) ;
            }
         }
      }
      writer->add_entry( "values", PEL_DoubleArray3D::create( 0, VALUES ) ) ;
   }
   
   // Active nodes:
   bool const has_inactive_nodes = ( ACTIVE_NODES.size() != 0 ) ;
   writer->add_entry( "has_inactive_nodes",
                      PEL_Bool::create( 0, has_inactive_nodes ) ) ;
   if( has_inactive_nodes )
   {
      writer->add_entry( "active_nodes",
                         PEL_BoolVector::create( 0, ACTIVE_NODES ) ) ;
   }
   
   // Saving imposed DOF
   bool const has_dbc = ( DBC_IDX.index_bound(0) != 0 ) ;
   writer->add_entry( "has_dbc",
                      PEL_Bool::create( 0, has_dbc ) ) ;
   if( has_dbc )
   {
      intArray2D dbc_indx( NB_NODES, NB_COMPS ) ;
      for( size_t i=0 ; i<NB_NODES ; ++i )
      {
         for( size_t j=0 ; j<NB_COMPS ; ++j )
         {
            size_t const idx = DBC_IDX( i,j ) ;
            if( idx == PEL::bad_index() )
            {
               dbc_indx(i,j) = -1 ;
            }
            else
            {
               PEL_ASSERT( idx< (size_t) PEL::max_int() ) ;
               dbc_indx(i,j) = (int) idx ;
            }
         }
      }
      writer->add_entry( "dbc_indx", PEL_IntArray2D::create( 0, dbc_indx ) ) ;
      writer->add_entry( "dbc_flag", PEL_BoolVector::create( 0, DBC_FLAG ) ) ;
      writer->add_entry( "dbc_values", PEL_DoubleVector::create( 0, DBC_VALUES ) ) ;
   }
   
   writer->finalize_object() ;

   PEL_CHECK_POST( save_state_POST( writer ) ) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: restore_state( PEL_ObjectReader* reader )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: restore_state" ) ;
   PEL_CHECK_PRE( restore_state_PRE( reader ) ) ;

   //??? here constraints are ignored and should not be
   PEL_ASSERT( CSTR == 0 ) ;

   reader->start_object_retrieval( "PDE_DiscreteField" ) ;

   // Retrieving stored datas :
   string const& read_name = reader->data_of_entry( "name" )->to_string() ;
   int const read_nb_nodes = reader->data_of_entry( "nb_nodes" )->to_int() ;
   int const read_nb_comps = reader->data_of_entry( "nb_components" )->to_int() ;
   int const read_nb_levels = reader->data_of_entry( "nb_levels" )->to_int() ;

   // Does some checks
   PEL_ASSERT( read_name==FNAME ) ;
   PEL_ASSERT( read_nb_nodes==(int) NB_NODES ) ;
   PEL_ASSERT( read_nb_comps==(int) NB_COMPS ) ;
   PEL_ASSERT( read_nb_levels==(int) STO_DEPTH ) ;

   // Retrieving values
   {
      VALUES = reader->data_of_entry( "values" )->to_double_array3D() ;
      PEL_ASSERT( VALUES.index_bound(0) == NB_NODES ) ;
      PEL_ASSERT( VALUES.index_bound(1) == NB_COMPS ) ;
      PEL_ASSERT( VALUES.index_bound(2) == STO_DEPTH ) ;
   }

   // Discretization representation:
   bool const has_inactive_nodes =
                  reader->data_of_entry( "has_inactive_nodes" )->to_bool() ;
   if( has_inactive_nodes )
   {
      ACTIVE_NODES = reader->data_of_entry( "active_nodes" )->to_bool_vector() ;
      PEL_ASSERT( ACTIVE_NODES.size() == NB_NODES ) ;
   }
   else
   {
      ACTIVE_NODES.re_initialize( 0 ) ;
   }   

   // Retrieving imposed DOF
   bool const has_dbc = reader->data_of_entry( "has_dbc" )->to_bool() ;
   if( has_dbc )
   {
      intArray2D const& dbc_indx =
                 reader->data_of_entry( "dbc_indx" )->to_int_array2D() ;
      PEL_ASSERT( dbc_indx.index_bound(0) == NB_NODES ) ;
      PEL_ASSERT( dbc_indx.index_bound(1) == NB_COMPS ) ;    
      DBC_FLAG = reader->data_of_entry( "dbc_flag" )->to_bool_vector() ;
      DBC_VALUES = reader->data_of_entry( "dbc_values" )->to_double_vector() ;
      size_t const nb_dbc = DBC_FLAG.size() ;
      PEL_ASSERT( DBC_VALUES.size() == nb_dbc ) ;
      DBC_IDX.re_initialize( NB_NODES, NB_COMPS ) ;
      for( size_t i=0 ; i<NB_NODES ; ++i )
      {
         for( size_t j=0 ; j<NB_COMPS ; ++j )
         {
            int const idx = dbc_indx( i,j ) ;
            if( idx == -1 )
            {
               DBC_IDX(i,j) = PEL::bad_index() ;
            }
            else
            {
               DBC_IDX(i,j) = (size_t) idx ;
               PEL_ASSERT( DBC_IDX(i,j)<nb_dbc ) ;
            }
         }
      }
   }
   else
   {
      DBC_IDX.re_initialize( 0, 0 ) ;
      DBC_FLAG.re_initialize( 0 ) ;
      DBC_VALUES.re_initialize( 0 ) ;
   }
   
   reader->end_object_retrieval() ;

   PEL_CHECK_POST( restore_state_POST( reader ) ) ;
}

//----------------------------------------------------------------------
void
PDE_DiscreteField:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_DiscreteField:: print" ) ;
   std::string space( indent_width, ' ' ) ;
   os << space << "PDE_DiscreteField n°" << ID
      << "   ->   " << FNAME << std::endl << space ;

   for( size_t n=0 ; n<NB_NODES ; n++ )
   {
      os << "   Node n°" << n ;
      if( is_distributed() )
      {
         os << "  (global n°" << DF->global_node_index(n) << ")" ;
      }
      os << std::endl << space ;
      std::ios_base::fmtflags original_flags = os.flags() ;
      os.setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
      os << std::setprecision(7) ;
      for( size_t ic=0 ; ic<NB_COMPS ; ic++ )
      {
         for( size_t level=0 ; level<storage_depth() ; level++ )
         {
           
            os << "         " ;
            os << "X(" << level << ")= " << std::setw(14) 
                << DOF_value( level, n, ic ) ;
         }
        
         if( DOF_has_imposed_value( n, ic ) )
         {
            os << "     " << "Xdbc=" << std::setw(14) 
                << DOF_imposed_value( n, ic ) ;
         }
         if( !node_is_active( n ) )
         {
            os << "     NOT ACTIVE" ;
         }        
         os << std::endl << space ;
      }
      os.flags( original_flags ) ;
   }
}

