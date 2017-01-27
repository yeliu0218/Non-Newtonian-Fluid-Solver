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

#include <PDE_SetOfBasisFunctions.hh>

#include <PDE_BasisFunction.hh>
#include <PDE_BasisFunctionCell.hh>
#include <PDE_CrossProcessBFNumbering.hh>
#include <PDE_ReferenceElement.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_IndexSet.hh>
#include <PEL_ObjectReader.hh>
#include <PEL_ObjectWriter.hh>
#include <PEL_Vector.hh>

#include <stringVector.hh>

#include <iostream>
#include <sstream>

using std::ostringstream ; using std::endl ;

//----------------------------------------------------------------------
PDE_SetOfBasisFunctions*
PDE_SetOfBasisFunctions:: create( PEL_Object* a_owner )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBasisFunctions:: create" ) ;

   PDE_SetOfBasisFunctions* result =
                            new PDE_SetOfBasisFunctions( a_owner ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->nb_basis_functions() == 0 ) ;
   PEL_CHECK_POST( FORALL( ( size_t e=0 ;
       e < result->nb_ref_elts_grps() ; ++e ),
       !result->is_valid( e ) ) ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
PDE_SetOfBasisFunctions:: PDE_SetOfBasisFunctions(
                                         PEL_Object* a_owner )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , REF_ELTS( 0 )
   , NUMBERING( 0 )
   , BFS_I( 0, PEL::bad_index() )
{
}

//----------------------------------------------------------------------
PDE_SetOfBasisFunctions:: ~PDE_SetOfBasisFunctions( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PDE_SetOfBasisFunctions:: add( PDE_BasisFunctionCell* bf, size_t e )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBasisFunctions:: add" ) ;
   PEL_CHECK_PRE( bf != 0 ) ;
   PEL_CHECK_PRE( bf->owner() == 0 ) ;
   PEL_CHECK_PRE( bf->id_number() == PEL::bad_index() ) ;
   PEL_CHECK_PRE( bf->ref_elts_grp_index() == PEL::bad_index() ) ;
   PEL_SAVEOLD( size_t, bfs, nb_basis_functions() ) ;
   PEL_SAVEOLD( size_t, bfse, nb_basis_functions( e ) ) ;

   bf->set_ref_elts_grp_index( e ) ;
   for( size_t i = BFS.size() ; i < e+1 ; ++i )
   {
      PEL_Vector* vec = PEL_Vector::create( this, 0 ) ;
      BFS.push_back( vec ) ;
   }
   BFS_I.re_initialize( BFS.size(), PEL::bad_index() ) ;

   bf->set_id_number( BFS[e]->index_limit() ) ;
   BFS[e]->append( bf ) ;
   bf->set_owner( this ) ;

   if( is_distributed() ) NUMBERING->set_unsynchronized_state() ;
   
   PEL_CHECK_POST( bf->owner() == this ) ;
   PEL_CHECK_POST( nb_basis_functions() == OLD( bfs ) + 1 ) ;
   PEL_CHECK_POST( nb_basis_functions( e ) == OLD( bfse ) + 1 ) ;
   PEL_CHECK_POST( bf->ref_elts_grp_index() == e ) ;
   PEL_CHECK_POST( bf->id_number()
                   == nb_basis_functions( bf->ref_elts_grp_index() ) - 1 ) ;
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i < nb_ref_elts_grps() ; ++i ),
                      !is_valid( i ) ) ) ;
}

//----------------------------------------------------------------------
void
PDE_SetOfBasisFunctions:: extend_ref_elts_grp( PEL_Vector const* elms,
                                               size_t& e )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBasisFunctions:: extend_ref_elts_grp" ) ;
   PEL_CHECK_PRE( elms != 0 ) ;

   size_t_vector elt_ids( elms->index_limit() ) ;
   for( size_t ii = 0 ; ii< elms->index_limit() ; ++ii )
   {
      PDE_ReferenceElement const* elt =
            static_cast< PDE_ReferenceElement* >( elms->at( ii ) ) ;
      elt_ids( ii ) = elt->id_number() ;
   }

   PEL_IndexSet* elts_ids_set =
                 PEL_IndexSet::create( 0, elt_ids, PEL::bad_index() ) ;
   
   e = index_of_ref_elts_grp( elts_ids_set ) ;
   if( e == PEL::bad_index() )
   {
      e = REF_ELTS.size() ;
      elts_ids_set->set_owner( this ) ;
      REF_ELTS.push_back( elts_ids_set ) ;
   }
   else
   {
      elts_ids_set->destroy() ; 
   }
   
   PEL_CHECK_POST( e < nb_ref_elts_grps() ) ;
}

//----------------------------------------------------------------------
size_t
PDE_SetOfBasisFunctions:: index_of_ref_elts_grp(
                                    PEL_IndexSet const* ref_elts ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBasisFunctions:: index_of_ref_elts_grp" ) ;

   size_t result = PEL::bad_index() ;

   for( size_t i=0 ; i<nb_ref_elts_grps() ; i++ )
   {
      PEL_IndexSet const* elts = REF_ELTS[ i ] ;
      if( ref_elts->is_equal( elts ) )
      {
         result = i ;
         break ;
      }
   }

   return( result ) ;
}

//----------------------------------------------------------------------
PDE_BasisFunction*
PDE_SetOfBasisFunctions:: item( size_t id_number,size_t e ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBasisFunctions:: item" ) ;
   PEL_CHECK_PRE( has( id_number, e ) ) ;

   PDE_BasisFunction* result =
                static_cast<PDE_BasisFunction*>( BFS[e]->at( id_number ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PDE_CrossProcessBFNumbering*
PDE_SetOfBasisFunctions:: cross_process_numbering( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBasisFunctions:: cross_process_numbering" ) ;
   PEL_CHECK_PRE( is_distributed() ) ;

   PDE_CrossProcessBFNumbering* result = NUMBERING ;

   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==this ) ;
   PEL_CHECK_POST( result->set_of_basis_functions()==this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_SetOfBasisFunctions:: build_cross_process_globalization(
                                PDE_CrossProcessBFNumbering* numbering )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBasisFunctions:: build_cross_process_globalization" ) ;
   PEL_CHECK_PRE( !is_distributed() ) ;
   PEL_CHECK_PRE( numbering != 0 ) ;
   PEL_CHECK_PRE( numbering->owner() == this ) ;

   NUMBERING = numbering ;
   NUMBERING->globalize() ;

   PEL_CHECK_POST( is_distributed() ) ;
   PEL_CHECK_POST( cross_process_numbering() == numbering ) ;
}

//----------------------------------------------------------------------
bool
PDE_SetOfBasisFunctions:: is_distributed( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBasisFunctions:: is_distributed" ) ;
   return( NUMBERING != 0 ) ;
}

//----------------------------------------------------------------------
bool
PDE_SetOfBasisFunctions:: has( size_t id_number, size_t e ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBasisFunctions:: has" ) ;
   PEL_CHECK_PRE( e < nb_ref_elts_grps() ) ;

   return( id_number < nb_basis_functions( e ) ) ;
}

//----------------------------------------------------------------------
size_t
PDE_SetOfBasisFunctions:: nb_basis_functions( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBasisFunctions:: nb_basis_functions" ) ;

   size_t result = 0 ;
   for( size_t e = 0 ; e < nb_ref_elts_grps() ; ++e )
   {
      result += nb_basis_functions( e ) ;
   }
   return( result ) ;
}


//----------------------------------------------------------------------
size_t
PDE_SetOfBasisFunctions:: nb_ref_elts_grps( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBasisFunctions:: nb_ref_elts_grps" ) ;

   return( REF_ELTS.size() ) ;
}

//----------------------------------------------------------------------
size_t
PDE_SetOfBasisFunctions:: nb_basis_functions( size_t e ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBasisFunctions:: nb_basis_functions" ) ;

   size_t result = 0 ;
   if( e < BFS.size() ) result = BFS[e]->index_limit() ;

   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_SetOfBasisFunctions:: start( size_t e ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBasisFunctions:: start" ) ;
   PEL_CHECK_PRE( e < nb_ref_elts_grps() ) ;

   if( e < BFS_I.size() ) BFS_I( e ) = 0 ;
}

//----------------------------------------------------------------------
void
PDE_SetOfBasisFunctions:: go_next( size_t e ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBasisFunctions:: go_next" ) ;
   PEL_CHECK_PRE( e < nb_ref_elts_grps() ) ;

   if( e < BFS_I.size() ) BFS_I( e )++ ;
}

//----------------------------------------------------------------------
bool
PDE_SetOfBasisFunctions:: is_valid( size_t e ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBasisFunctions:: is_valid" ) ;
   PEL_CHECK_PRE( e < nb_ref_elts_grps() ) ;

   bool result = ( e < BFS_I.size() )
              && ( BFS_I( e ) < nb_basis_functions( e ) ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
PDE_BasisFunction*
PDE_SetOfBasisFunctions:: item( size_t e ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfBasisFunctions:: item()" ) ;
   PEL_CHECK_PRE( is_valid( e ) ) ;

   PDE_BasisFunction* result =
             static_cast<PDE_BasisFunction*>( BFS[e]->at( BFS_I(e) ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_SetOfBasisFunctions:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string space( indent_width, ' ' ) ;
   os << space << nb_basis_functions() << " basis functions : "
      << std::endl ;
   for( size_t e=0 ; e < nb_ref_elts_grps() ; e++ )
   {
      os << space << "   **elements de reference : " ;
      PEL_IndexSet const* inds = REF_ELTS[ e ] ;
      inds->print( os, 1 ) ;
      PEL::out() << std::endl ;
      for( start( e ) ; is_valid( e ) ; go_next( e ) )
      {
         os << space << "     - id_number \"" << item( e )->id_number() << "\" : " ;
         os << std::endl ;
         item( e )->print( os, indent_width+7 ) ;
      }
   }
}
