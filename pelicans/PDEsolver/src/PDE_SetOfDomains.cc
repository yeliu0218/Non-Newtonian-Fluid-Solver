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

#include <PDE_SetOfDomains.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <GE_Color.hh>
#include <GE_Mpolyhedron.hh>
#include <GE_Point.hh>

#include <PDE_BoundFE.hh>
#include <PDE_DomainBuilder.hh>
#include <PDE_DomainAndFields.hh>
#include <PDE_GridFE.hh>
#include <PDE_InterfaceBuilder.hh>
#include <PDE_LocalFEmortarSide.hh>

#include <iostream>
#include <sstream>

using std::cout ;
using std::endl ;
using std::ostringstream ;

//----------------------------------------------------------------------------
PDE_SetOfDomains*
PDE_SetOfDomains:: create( PEL_Object* a_owner,
                           PEL_ModuleExplorer* exp )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDomains:: create" ) ;
   PEL_CHECK_PRE( exp != 0 ) ;

   PDE_SetOfDomains* result = new PDE_SetOfDomains( a_owner, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PDE_SetOfDomains:: PDE_SetOfDomains( PEL_Object* a_owner,
                                     PEL_ModuleExplorer* exp )
//----------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , VERB( exp->int_data( "verbose_level" ) ) 
   , DOMAINS( PEL_Vector::create( this, 0 ) )
   , INTERFACES( PEL_Vector::create( this, 0 ) )
   , I_DOM_0( 0 )
   , I_DOM_1( 0 )
{
   PEL_ModuleExplorer* ee = 
                exp->create_subexplorer( 0, "list_of_PDE_DomainAndFields" ) ;

   ee->start_module_iterator() ;
   for( ; ee->is_valid_module() ; ee->go_next_module() )
   {
      PEL_ModuleExplorer* sexp = ee->create_subexplorer( this ) ;
      std::string nn = sexp->name() ;
      if( nn.substr( 0, 19 ) != "PDE_DomainAndFields" )
      {
         PEL_Error::object()->raise_plain( "invalid module name : " + nn ) ;
      }

      PDE_DomainAndFields* daf = PDE_DomainAndFields::create( this, sexp ) ;
      DOMAINS->append( daf ) ;
   }
   ee->destroy() ;

   if( exp->has_module( "list_of_PDE_InterfaceAndFields" ) )
   {
      ee = exp->create_subexplorer( 0, "list_of_PDE_InterfaceAndFields" ) ;

      ee->start_module_iterator() ;
      for( ; ee->is_valid_module() ; ee->go_next_module() )
      {
         PEL_ModuleExplorer* sexp = ee->create_subexplorer( this ) ;
         std::string nn = sexp->name() ;
         if( nn.substr( 0, 22 )  != "PDE_InterfaceAndFields" )
         {   
            PEL_Error::object()->raise_plain( "invalid module name : "+nn ) ;
         }

         PDE_InterfaceAndFields* iaf =
                                PDE_InterfaceAndFields::create( this, sexp ) ;

         INTERFACES->append( iaf ) ;
         for( size_t i=0 ; i<nb_domains() ; ++i)
         {
            if( domain(i) == iaf->adjacent_domain( 0 ) ) 
               I_DOM_0.append( i ) ;
            else if( domain(i) == iaf->adjacent_domain( 1 ) ) 
               I_DOM_1.append( i ) ;
            else
               PEL_Error::object()->raise_plain( "invalid interface" ) ;
         }
         PEL_ASSERT( INTERFACES->index_limit() == I_DOM_0.size() ) ;
         PEL_ASSERT( INTERFACES->index_limit() == I_DOM_1.size() ) ;
      }

      ee->destroy() ;
   }

   if( exp->has_module( "list_of_conformal_adjacencies" ) )
   {
      ee = exp->create_subexplorer( 0, "list_of_conformal_adjacencies" ) ;

      if( VERB!=0 ) PEL::out()<< "*** Conformal adjacencies" << endl ;

      ee->start_module_iterator() ;
      for( ; ee->is_valid_module() ; ee->go_next_module() )
      {
         PEL_ModuleExplorer* sexp = ee->create_subexplorer( 0 ) ;
         resolve_conformal_adjacencies( sexp ) ;
         sexp->destroy() ;
      }
      ee->destroy() ;
   }
}


//----------------------------------------------------------------------------
PDE_SetOfDomains:: ~PDE_SetOfDomains( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
size_t
PDE_SetOfDomains:: nb_domains( void ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDomains:: nb_domains" ) ;

   return( DOMAINS->index_limit() ) ;
}

//----------------------------------------------------------------------------
PDE_DomainAndFields*
PDE_SetOfDomains:: domain( size_t i ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDomains:: domain" ) ;
   PEL_CHECK_PRE( i < nb_domains() ) ;

   PDE_DomainAndFields* result = 
                         static_cast<PDE_DomainAndFields*>( DOMAINS->at(i) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PDE_DomainAndFields*
PDE_SetOfDomains:: domain( std::string const& name ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDomains:: domain" ) ;

   PDE_DomainAndFields* result = 0 ;
   PEL_VectorIterator* it = DOMAINS->create_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PDE_DomainAndFields* dom = 
                           static_cast<PDE_DomainAndFields*>( it->item() ) ;
      if( dom->name() == name )
      {
         result = dom ;
         break ;
      }
   }
   it->destroy() ;

   if( result == 0 )
   {
      ostringstream mesg ;
      mesg << "PDE_SetOfDomains : " << endl ;
      mesg << "   no \"PDE_DomainAndFields\" object of name" << endl ;
      mesg << "   \"" << name << "\" is recorded" ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
size_t
PDE_SetOfDomains:: nb_interfaces( void ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDomains:: nb_interfaces" ) ;

   return( INTERFACES->index_limit() ) ;
}

//----------------------------------------------------------------------------
PDE_InterfaceAndFields*
PDE_SetOfDomains:: interface( size_t i ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDomains:: interface" ) ;
   PEL_CHECK_PRE( i < nb_interfaces() ) ;

   PDE_InterfaceAndFields* result =
                static_cast<PDE_InterfaceAndFields*>( INTERFACES->at( i ) ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
PDE_InterfaceAndFields*
PDE_SetOfDomains:: interface( std::string const& name ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDomains:: interface" ) ;

   PDE_InterfaceAndFields* result = 0 ;
   PEL_VectorIterator* it = INTERFACES->create_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PDE_InterfaceAndFields* interf = 
                           static_cast<PDE_InterfaceAndFields*>( it->item() ) ;
      if( interf->name() == name )
      {
         result = interf ;
         break ;
      }
   }
   it->destroy() ;

   if( result == 0 )
   {
      ostringstream mesg ;
      mesg << "PDE_SetOfDomains : " << endl ;
      mesg << "   no \"PDE_InterfaceAndFields\" object of name" << endl ;
      mesg << "   \"" << name << "\" is recorded" ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == this ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
size_t
PDE_SetOfDomains:: index_of_interface_adjacent_domain( size_t i,
                                                       size_t i_adj ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDomains:: index_of_interface_adjacent_domain" ) ;
   PEL_CHECK_PRE( i < nb_interfaces() ) ;
   PEL_CHECK_PRE( i_adj==0 || i_adj==1 ) ;

   size_t result = ( i_adj==0 ? I_DOM_0( i ) : I_DOM_1( i ) ) ;

   PEL_CHECK_POST( result < nb_domains() ) ;
   PEL_CHECK_POST( domain( result ) == 
                   interface( i )->adjacent_domain( i_adj ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
PDE_SetOfDomains:: resolve_conformal_adjacencies( 
                                               PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_SetOfDomains:: resolve_conformal_adjacencies" ) ;

   std::string const& nn0 = exp->string_data( "adjacent_domain_0" ) ;
   std::string const& nn1 = exp->string_data( "adjacent_domain_1" ) ;

   PDE_DomainBuilder* builder_0 = PDE_DomainBuilder::object( nn0 ) ;
   PDE_DomainBuilder* builder_1 = PDE_DomainBuilder::object( nn1 ) ;

   PEL_VectorIterator* it0 = PEL_VectorIterator::create( 0, 
                                builder_0->finite_element_grid()->bounds() ) ;
   PEL_VectorIterator* it1 = PEL_VectorIterator::create( 0,
                                builder_1->finite_element_grid()->bounds() ) ;

   size_t nb_adjacent_bounds = 0 ;
   for( it0->start() ; it0->is_valid() ; it0->go_next() )
   {
      PDE_BoundFE* bd0 = static_cast<PDE_BoundFE*>( it0->item() ) ;
      GE_Mpolyhedron* poly0 = bd0->polyhedron() ;

      GE_Point const* pt0 = poly0->center() ;
      for( it1->start() ; it1->is_valid() ; it1->go_next() ) 
      {
         PDE_BoundFE* bd1 = static_cast<PDE_BoundFE*>( it1->item() ) ;
         GE_Mpolyhedron const* poly1 = bd1->polyhedron() ;
         GE_Point const* pt1 = poly1->center() ;
         if( pt1->distance( pt0 ) < 1.e-8 ) //??????????????
         {
            poly0->reorder_vertices_according_to( poly1 ) ;
            bd0->insert_adjacent_bound( bd1 ) ;
            bd1->insert_adjacent_bound( bd0 ) ;
            nb_adjacent_bounds++ ;
            break ;
         }
      }
   }
   it0->destroy() ;
   it1->destroy() ;

   if( nb_adjacent_bounds==0 )
   {
      PEL_Error::object()->raise_plain( "no adjacent bound found between \"" +
                                        nn0 + "\" and \"" + nn1 + "\"" ) ;
   }

   builder_0->insert_conformal_adjacency( builder_1 ) ;
   builder_1->insert_conformal_adjacency( builder_0 ) ;

   if( VERB!=0 )
      PEL::out() << "    " << nb_adjacent_bounds 
                 << " adjacent bounds between \"" 
                 << nn0 << "\" and \"" << nn1 << "\"" << endl ;
}
