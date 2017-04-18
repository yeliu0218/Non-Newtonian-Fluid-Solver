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

#include <FE_OneStepIterationOpen.hh>

#include <PEL_DistributedPartition.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_assertions.hh>

#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_LinkDOF2Unknown.hh>
#include <PDE_SystemNumbering.hh>

#include <FE_TimeIterator.hh>

#include <iostream>
#include <sstream>

using std::cout ; using std::endl ;
using std::string ; using std::ostringstream ;

std::map< std::string, FE_OneStepIterationOpen* > 
FE_OneStepIterationOpen:: OBJS ;

struct FE_OneStepIterationOpen_ERROR {
   static void n0( std::string const& a_name ) ;
   static void n1( std::string const& a_name ) ;
   static void n2( void ) ;
} ;

//---------------------------------------------------------------------------
FE_OneStepIterationOpen:: FE_OneStepIterationOpen( std::string const& a_name )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_name )
{
}

//---------------------------------------------------------------------------
FE_OneStepIterationOpen:: FE_OneStepIterationOpen(
                                             PEL_Object* a_owner,
                                             PDE_DomainAndFields const* dom,
                                             PEL_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom,  exp )
   , NN( "" )
   , DOM( dom )
{
   if( exp->has_entry( "name" ) )
   {
      NN = exp->string_data( "name" ) ;
      if( NN.substr(0,7) == "unnamed" ) FE_OneStepIterationOpen_ERROR::n2() ;
      if( OBJS.count( NN ) != 0 ) FE_OneStepIterationOpen_ERROR::n0( NN ) ;
      OBJS[ NN ] = this ;
   }
   else
   {
      static size_t i_unnamed = 0 ;
      ostringstream nnn ;
      nnn << "unnamed_" << i_unnamed ;
      NN = nnn.str() ;
      PEL_ASSERT( OBJS.count( NN ) == 0 ) ;
      OBJS[ NN ] = this ;
      ++i_unnamed ;
   }
}

//---------------------------------------------------------------------------
FE_OneStepIterationOpen:: ~FE_OneStepIterationOpen( void )
//---------------------------------------------------------------------------
{
   OBJS.erase( NN ) ;
}

//--------------------------------------------------------------------------
FE_OneStepIterationOpen*
FE_OneStepIterationOpen:: object( std::string const& a_name )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIterationOpen:: object" ) ;
   std::map<std::string,FE_OneStepIterationOpen*>::const_iterator it = 
                                                         OBJS.find( a_name ) ;
   if( it == OBJS.end() ) FE_OneStepIterationOpen_ERROR::n1( a_name ) ;
   
   FE_OneStepIterationOpen* result = (*it).second ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->name() == a_name ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
FE_OneStepIterationOpen:: nb_objects( void )
//---------------------------------------------------------------------------
{
   size_t result = OBJS.size() ;
   return( result ) ;
}

//---------------------------------------------------------------------------
FE_OneStepIterationOpen*
FE_OneStepIterationOpen:: object( size_t i )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIterationOpen:: object" ) ;
   PEL_CHECK_PRE( i < nb_objects() ) ;
   
   std::map<std::string,FE_OneStepIterationOpen*>::const_iterator it = 
                                                              OBJS.begin() ;
   size_t i_cur=0 ;
   for( ; it!=OBJS.end() ; ++it )
   {
      if( i_cur == i ) break ;
      ++i_cur ;
   }
   PEL_CHECK( it != OBJS.end() ) ;
   FE_OneStepIterationOpen* result = (*it).second ;
   
   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
std::string const&
FE_OneStepIterationOpen:: name( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIterationOpen:: name" ) ;

   return( NN ) ;
}

//---------------------------------------------------------------------------
PDE_DomainAndFields const*
FE_OneStepIterationOpen:: domain( void ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIterationOpen:: domain" ) ;
   
   PDE_DomainAndFields const* result = DOM ;
   
   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
size_t
FE_OneStepIterationOpen:: index_of_field( std::string const& fname ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIterationOpen:: index_of_field" ) ; 
   
   size_t result = PEL::bad_index() ;
   for( size_t i=0 ; i<nb_unknowns() ; ++i )
   {
      if( field( i )->name() == fname )
      {
         result = i ;
         break ;
      }
   }
   
   PEL_CHECK_POST( result==PEL::bad_index() ||
                   field( result )->name() == fname ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
FE_OneStepIterationOpen:: build_function_and_jacobian( 
                                              FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIterationOpen:: build_function_and_jacobian" ) ;
   
   PEL_Error::object()->raise_not_implemented( this, 
                                               "build_function_and_jacobian" ) ;
}

//---------------------------------------------------------------------------
LA_SeqVector const*
FE_OneStepIterationOpen:: create_function( PEL_Object* a_owner,
                                           size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIterationOpen:: create_function" ) ;
   PEL_CHECK_PRE( create_function_PRE( a_owner, i_unk ) ) ;
   
   PEL_Error::object()->raise_not_implemented( this, "create_function" ) ;
   
   LA_SeqVector const* result = 0 ;
   PEL_CHECK_POST( create_function_POST( result, a_owner, i_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
LA_SeqMatrix const*
FE_OneStepIterationOpen:: create_jacobian( PEL_Object* a_owner,
                                           size_t i_eq,
                                           size_t j_unk ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIterationOpen:: create_jacobian" ) ;
   PEL_CHECK_PRE( create_jacobian_PRE( a_owner, i_eq, j_unk ) ) ;
   
   PEL_Error::object()->raise_not_implemented( this, "create_jacobian" ) ;
   
   LA_SeqMatrix const* result = 0 ;
   PEL_CHECK_POST( create_jacobian_POST( result, a_owner, i_eq, j_unk ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
void
FE_OneStepIterationOpen:: assemble_contribution( FE_TimeIterator const* t_it,
                                                 LA_Matrix* matrix,
                                                 LA_Vector* vector,
                                                 PDE_SystemNumbering const* nmb,
                                                 size_t i_link ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIterationOpen:: assemble_contribution" ) ;
   PEL_CHECK_PRE( assemble_contribution_PRE( t_it, matrix, vector, 
                                             nmb, i_link ) ) ;
   
   PEL_Error::object()->raise_not_implemented( this, "assemble_contribution" ) ;
}

//---------------------------------------------------------------------------
void
FE_OneStepIterationOpen:: update_DOFs( LA_Vector* vector,
                                       PDE_SystemNumbering const* nmb,
                                       size_t i_link ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_OneStepIterationOpen:: update_DOFs" ) ;
   PEL_CHECK_PRE( update_DOFs_PRE( vector, nmb, i_link ) ) ;
   
   PEL_Error::object()->raise_not_implemented( this, "update_DOFs" ) ;
}

//---------------------------------------------------------------------------
void
FE_OneStepIterationOpen:: check_assembled_system( 
                                PDE_SystemNumbering const* nmb,
                                std::string const& function_name ) const
//---------------------------------------------------------------------------
{
   if( nmb == 0 )
   {
      ostringstream mesg ;
      mesg << endl << "*** " << type_name() << endl << endl ;
      mesg << "   an instance of" << endl ;
      mesg << "      PDE_SystemNumbering" << endl ;
      mesg << "   should be available when calling" << endl ;
      mesg << "      " << function_name  ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }
}

//---------------------------------------------------------------------------
bool
FE_OneStepIterationOpen:: field_PRE( size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( i_unk < nb_unknowns() ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
FE_OneStepIterationOpen:: field_POST( PDE_DiscreteField const* result,
                                      size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result == link_DOF_2_unknown( i_unk )->field() ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
FE_OneStepIterationOpen:: level_of_field_PRE( size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( i_unk < nb_unknowns() ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
FE_OneStepIterationOpen:: level_of_field_POST( size_t result, 
                                               size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( result < field( i_unk )->storage_depth() ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
FE_OneStepIterationOpen:: link_DOF_2_unknown_PRE( size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( i_unk < nb_unknowns() ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
FE_OneStepIterationOpen:: link_DOF_2_unknown_POST( 
                                 PDE_LinkDOF2Unknown const* result,
                                 size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->field( ) == field( i_unk ) ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
FE_OneStepIterationOpen:: create_function_PRE( PEL_Object* a_owner,
                                                size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( i_unk < nb_unknowns() ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
FE_OneStepIterationOpen:: create_function_POST( LA_SeqVector const* result,
                                                PEL_Object* a_owner,
                                                size_t i_unk ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( result->nb_rows() == 
               link_DOF_2_unknown( i_unk )->unknown_vector_size() ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
FE_OneStepIterationOpen:: create_jacobian_PRE( PEL_Object* a_owner,
                                                size_t i_eq,
                                                size_t j_unk) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( i_eq < nb_unknowns() ) ;
   PEL_ASSERT( j_unk < nb_unknowns() ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
FE_OneStepIterationOpen:: create_jacobian_POST( LA_SeqMatrix const* result,
                                                 PEL_Object* a_owner,
                                                 size_t i_eq,
                                                 size_t j_unk) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( result != 0 ) ;
   PEL_ASSERT( result->owner() == a_owner ) ;
   PEL_ASSERT( result->nb_rows() == 
               link_DOF_2_unknown( i_eq )->unknown_vector_size() ) ;
   PEL_ASSERT( result->nb_cols() == 
               link_DOF_2_unknown( j_unk )->unknown_vector_size() ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
FE_OneStepIterationOpen:: assemble_contribution_PRE( 
                                                FE_TimeIterator const* t_it,
                                                LA_Matrix* matrix,
                                                LA_Vector* vector,
                                                PDE_SystemNumbering const* nmb,
                                                size_t i_link  ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( t_it != 0 ) ;
   PEL_ASSERT( t_it->is_started() ) ;
   PEL_ASSERT( !t_it->is_finished() ) ;
   PEL_ASSERT( nmb != 0 ) ;
   PEL_ASSERT( matrix != 0 ) ;
   PEL_ASSERT( vector != 0 ) ;
   PEL_ASSERT( matrix->nb_rows() == nmb->nb_global_unknowns() ) ;
   PEL_ASSERT( matrix->nb_cols() == nmb->nb_global_unknowns() ) ;
   PEL_ASSERT( vector->nb_rows() == nmb->nb_global_unknowns() ) ;
   PEL_ASSERT( matrix->state() != LA::NotSync_set ) ;
   PEL_ASSERT( vector->state() != LA::NotSync_set ) ;
   PEL_ASSERT( i_link < nmb->nb_links() ) ;
   return( true ) ;
}

//---------------------------------------------------------------------------
bool
FE_OneStepIterationOpen:: update_DOFs_PRE( LA_Vector* vector,
                                           PDE_SystemNumbering const* nmb,
                                           size_t i_link ) const
//---------------------------------------------------------------------------
{
   PEL_ASSERT( vector != 0 ) ;
   PEL_ASSERT( nmb != 0 ) ;
   PEL_ASSERT( nmb->scatters_are_defined() ) ;
   PEL_ASSERT( i_link < nmb->nb_links() ) ;
   PEL_ASSERT( vector != 0 ) ;
   PEL_ASSERT( vector->is_synchronized() ) ;
   PEL_ASSERT( vector->implementation() == 
               nmb->scatter( i_link )->implementation() ) ;
   PEL_ASSERT(
      IMPLIES( vector->row_distribution() != 0,
               nmb->scatter( i_link )->distribution() != 0 ) ) ;
   PEL_ASSERT(
      IMPLIES( vector->row_distribution() != 0,
               vector->row_distribution()->is_compatible( 
                       nmb->scatter( i_link )->distribution() ) ) ) ;
   return( true ) ;
}

//internal-------------------------------------------------------------------
void
FE_OneStepIterationOpen_ERROR:: n0( std::string const& a_name )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "FE_OneStepIterationOpen:" << endl ;
   mesg << "   attempt to create two instances with the" << endl ;
   mesg << "   same name :\"" << a_name << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
FE_OneStepIterationOpen_ERROR:: n1( std::string const& a_name )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "FE_OneStepIterationOpen:" << endl ;
   mesg << "   there is no recorded instance of name :\"" << a_name << "\"" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal-------------------------------------------------------------------
void
FE_OneStepIterationOpen_ERROR:: n2( void )
//internal-------------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "FE_OneStepIterationOpen:" << endl ;
   mesg << "   the name of an instance cannot start with \"unnamed\"" << endl ;
   mesg << "   (\"unnamed\" is reserved for internal purposes)" ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

