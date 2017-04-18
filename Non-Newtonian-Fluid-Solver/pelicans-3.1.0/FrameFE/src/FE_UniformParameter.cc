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

#include <FE_UniformParameter.hh>

#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_ObjectReader.hh>
#include <PEL_ObjectWriter.hh>

FE_UniformParameter const* 
FE_UniformParameter:: PROTOTYPE = new FE_UniformParameter() ;

//-------------------------------------------------------------------------
FE_UniformParameter:: FE_UniformParameter( void )
//-------------------------------------------------------------------------
   : FE_Parameter( "FE_UniformParameter" )
   , VAL( 0 )
{
}

//-------------------------------------------------------------------------
FE_UniformParameter*
FE_UniformParameter:: create( PEL_Object* a_owner,
                              std::string const& a_name,
                              doubleVector const& values )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: create" ) ;

   FE_UniformParameter* result = new FE_UniformParameter( a_owner, 
                                                          a_name, 
                                                          values ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_UniformParameter:: FE_UniformParameter( PEL_Object* a_owner,
                                           std::string const& a_name,
                                           doubleVector const& values )
//-------------------------------------------------------------------------
   : FE_Parameter( a_owner, a_name )
   , VAL( values )
{
}

//-------------------------------------------------------------------------
FE_UniformParameter*
FE_UniformParameter:: create_replica( PEL_Object* a_owner,
                                      PDE_DomainAndFields const* dom,
                                      PEL_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   FE_UniformParameter* result = new FE_UniformParameter( a_owner, dom, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
FE_UniformParameter:: FE_UniformParameter( PEL_Object* a_owner,
                                           PDE_DomainAndFields const* dom,
                                           PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : FE_Parameter( a_owner, exp->string_data( "name" ) )
   , VAL( exp->doubleVector_data( "value" ) )
{
}

//-------------------------------------------------------------------------
FE_UniformParameter:: ~FE_UniformParameter( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
FE_UniformParameter:: reset( PEL_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: reset" ) ;
   if( exp != 0 )
   {
      if( !exp->has_entry( "value" ) )
      {
         std::string msg ;
         msg += "*** FE_UniformParameter<" + name() + "> error\n" ;
         msg += "    reset should be called with :\n" ;
         msg += "       - \"exp\" != 0\n" ;
         msg += "       - \"exp\" contains a doubleVector data\n" ;
         msg += "                 of name \"value\"" ;
         PEL_Error::object()->raise_internal( msg ) ;
      }
      VAL = exp->doubleVector_data( "value" ) ;
   }
}

//-------------------------------------------------------------------------
size_t
FE_UniformParameter:: nb_components( void ) const
//-------------------------------------------------------------------------
{
   size_t result = VAL.size() ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: cell_value( FE_TimeIterator const* t_it,
                                  PDE_LocalFEcell const* fe,
                                  size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: cell_value" ) ;
   PEL_CHECK_PRE( cell_value_PRE( t_it, fe, ic ) ) ;

   double result = VAL(ic) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: cell_gradient( FE_TimeIterator const* t_it,
                                     PDE_LocalFEcell const* fe,
                                     size_t a,
                                     size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: cell_gradient" ) ;
   PEL_CHECK_PRE( cell_gradient_PRE( t_it, fe, a, ic ) ) ;

   return( 0.0 ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: cell_value_at_IP( FE_TimeIterator const* t_it,
                                        PDE_LocalFEcell const* fe,
                                        size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: cell_value_at_IP" ) ;
   PEL_CHECK_PRE( cell_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double result = VAL(ic) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: cell_gradient_at_IP( FE_TimeIterator const* t_it,
                                           PDE_LocalFEcell const* fe,
                                           size_t a,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: cell_gradient_at_IP" ) ;
   PEL_CHECK_PRE( cell_gradient_at_IP_PRE( t_it, fe, a, ic ) ) ;

   return( 0.0 ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: cell_value_at_pt( FE_TimeIterator const* t_it,
                                        PDE_LocalFEcell const* fe,
                                        size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: cell_value_at_pt" ) ;
   PEL_CHECK_PRE( cell_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double result = VAL(ic) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: cell_gradient_at_pt( FE_TimeIterator const* t_it,
                                           PDE_LocalFEcell const* fe,
                                           size_t a,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: cell_gradient_at_pt" ) ;
   PEL_CHECK_PRE( cell_gradient_at_pt_PRE( t_it, fe, a, ic ) ) ;

   return( 0.0 ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: bound_value( FE_TimeIterator const* t_it,
                                   PDE_LocalFEbound const* fe,
                                   size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: bound_value" ) ;
   PEL_CHECK_PRE( bound_value_PRE( t_it, fe, ic ) ) ;

   double result = VAL(ic) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: bound_gradient( FE_TimeIterator const* t_it,
                                      PDE_LocalFEbound const* fe,
                                      size_t a,
                                      size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: bound_gradient" ) ;
   PEL_CHECK_PRE( bound_gradient_PRE( t_it, fe, a, ic ) ) ;

   return( 0.0 ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: bound_value_at_IP( FE_TimeIterator const* t_it,
                                         PDE_LocalFEbound const* fe,
                                         size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: bound_value_at_IP" ) ;
   PEL_CHECK_PRE( bound_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double result = VAL(ic) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: bound_gradient_at_IP( FE_TimeIterator const* t_it,
                                            PDE_LocalFEbound const* fe,
                                            size_t a,
                                            size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: bound_gradient_at_IP" ) ;
   PEL_CHECK_PRE( bound_gradient_at_IP_PRE( t_it, fe, a, ic ) ) ;

   return( 0.0 ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: bound_value_at_pt( FE_TimeIterator const* t_it,
                                         PDE_LocalFEbound const* fe,
                                         size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: bound_value_at_pt" ) ;
   PEL_CHECK_PRE( bound_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double result = VAL(ic) ;

   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: bound_gradient_at_pt( FE_TimeIterator const* t_it,
                                            PDE_LocalFEbound const* fe,
                                            size_t a,
                                            size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: bound_gradient_at_pt" ) ;
   PEL_CHECK_PRE( bound_gradient_at_pt_PRE( t_it, fe, a, ic ) ) ;

   return( 0.0 ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: side_value( FE_TimeIterator const* t_it,
                                  PDE_CursorFEside const* fe,
                                  size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: side_value" ) ;
   PEL_CHECK_PRE( side_value_PRE( t_it, fe, ic ) ) ;

   double result = VAL(ic) ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: side_gradient( FE_TimeIterator const* t_it,
                                     PDE_CursorFEside const* fe,
                                     size_t a,
                                     size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: side_gradient" ) ;
   PEL_CHECK_PRE( side_gradient_PRE( t_it, fe, a, ic ) ) ;

   return( 0.0 ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: side_value_at_pt( FE_TimeIterator const* t_it,
                                          PDE_CursorFEside const* fe,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: side_value_at_pt" ) ;
   PEL_CHECK_PRE( side_value_at_pt_PRE( t_it, fe, ic ) ) ;

   double result = VAL(ic) ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: side_gradient_at_pt( FE_TimeIterator const* t_it,
                                           PDE_CursorFEside const* fe,
                                           size_t a,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: side_gradient_at_pt" ) ;
   PEL_CHECK_PRE( side_gradient_at_pt_PRE( t_it, fe, a, ic ) ) ;

   return( 0.0 ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: side_value_at_IP( FE_TimeIterator const* t_it,
                                          PDE_CursorFEside const* fe,
                                          size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: side_value_at_IP" ) ;
   PEL_CHECK_PRE( side_value_at_IP_PRE( t_it, fe, ic ) ) ;

   double result = VAL(ic) ;
   
   return( result ) ;
}

//-------------------------------------------------------------------------
double
FE_UniformParameter:: side_gradient_at_IP( FE_TimeIterator const* t_it,
                                           PDE_CursorFEside const* fe,
                                           size_t a,
                                           size_t ic ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: side_gradient_at_IP" ) ;
   PEL_CHECK_PRE( side_gradient_at_IP_PRE( t_it, fe, a, ic ) ) ;

   return( 0.0 ) ;
}

//----------------------------------------------------------------------
void
FE_UniformParameter:: save_state( PEL_ObjectWriter* writer ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: save_state" ) ;
   PEL_CHECK_PRE( save_state_PRE( writer ) ) ;
   
   writer->start_new_object( "FE_UniformParameter" ) ;

   // Saving current time
   writer->add_entry( "VAL", PEL_DoubleVector::create( 0, VAL ) ) ;
   
   writer->finalize_object() ;

   PEL_CHECK_POST( save_state_POST( writer ) ) ;
}

//----------------------------------------------------------------------
void
FE_UniformParameter:: restore_state( PEL_ObjectReader* reader )
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_UniformParameter:: restore_state" ) ;
   PEL_CHECK_PRE( restore_state_PRE( reader ) ) ;

   reader->start_object_retrieval( "FE_UniformParameter" ) ;
   
   VAL = reader->data_of_entry( "VAL" )->to_double_vector() ;
   reader->end_object_retrieval() ;

   PEL_CHECK_POST( restore_state_POST( reader ) ) ;
}
