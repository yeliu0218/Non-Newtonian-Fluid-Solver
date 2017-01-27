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

#include <DOC_Package.hh>
#include <DOC_Class.hh>
#include <DOC_Tools.hh>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <PEL_Root.hh>

using std::cerr ;
using std::string ;

std::string DOC_Package::application_name ="Unknown application" ;


//--------------------------------------------------------------------
DOC_Package*
DOC_Package:: create( std::string const& a_name,
                      DOC_Package* a_father,
                      std::string const& comment ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Package:: create" ) ;
   PEL_CHECK( a_father!=0 ) ;
   
   DOC_Package * result ;
   
   result = a_father->find_sub_package( a_name ) ;
   if( result==0 )
   {
      result = new DOC_Package( packages(), a_name, a_father, comment ) ;
      packages()->append( result ) ;
   }
   PEL_CHECK( result!=0 ) ;
   PEL_CHECK( result->name()==a_name ) ;
   PEL_CHECK( result->father()==a_father ) ;
   
   return result ;
}



//--------------------------------------------------------------------
DOC_Package const*
DOC_Package:: attached_to( DOC_Class const* classe,
                           std::string& comment_court  ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Package:: attached_to" ) ;
   PEL_CHECK_PRE( classe!=0 ) ;
   
   // Retourne le package auquel est associee la classe
   DOC_Package * result = 0 ;
   PEL_Iterator* it = packages()->create_iterator(0) ;
   
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      DOC_Package*p = static_cast<DOC_Package*>( it->item() ) ;
      size_t j=p->elements.index_of(classe->name()) ;
      if( j<p->elements.size() )
      {
         result=p ;
         comment_court = p->element_comments(j) ;
         break ;
      }

   }
// *** commented out so that "No package items
// *** do not appear in the documentation
// **    if( result==0 )
// ***   {
// *** //       std::cout << "miscellanous requested for " << classe->name() 
// *** //                 << std::endl;
// ***      result = miscellanous() ;
// ***   }
// ***
   it->destroy() ;
   return result ;
}



//--------------------------------------------------------------------
std::string const&
DOC_Package:: name( void ) const
//--------------------------------------------------------------------
{
   return my_name ;
}


//--------------------------------------------------------------------
DOC_Package const*
DOC_Package:: top( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Package:: top" ) ;
   PEL_CHECK_PRE( father()!=0 ) ;
   
   DOC_Package const* result = this ;
   while( result->my_father != master() )
      result = result->my_father ;
   
   PEL_CHECK_POST( result->father()==master() ) ;
   return result ;
}


//--------------------------------------------------------------------
DOC_Package*
DOC_Package:: father( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Package:: father" ) ;
   PEL_CHECK_PRE( this != master() ) ;
   
   DOC_Package* result = my_father ;
   
   PEL_CHECK_POST( result != 0 ) ;
   return result ;
}


//--------------------------------------------------------------------
DOC_Package* 
DOC_Package:: find( std::string const& a_name )
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Package:: find" ) ;

   DOC_Package * result = 0 ;
   PEL_Iterator *it = packages()->create_iterator(0) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      DOC_Package*p = static_cast<DOC_Package*>( it->item() ) ;
      if( p->name() == a_name )
      {
         result = p ;
         break ;
      }
   }
   it->destroy() ;
   
   PEL_CHECK_POST( IMPLIES( result!=0, result->name() == a_name ) ) ;
   
   return result ;
}


//--------------------------------------------------------------------
DOC_Package* 
DOC_Package:: find_sub_package( std::string const& a_name ) const
//--------------------------------------------------------------------
{
   DOC_Package * result = 0 ;
   PEL_Iterator *it = subpackages->create_iterator(0) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      DOC_Package*p = static_cast<DOC_Package*>( it->item() ) ;
      if( p->name() == a_name )
      {
         result = p ;
         break ;
      }
   }
   it->destroy() ;
   return result ;
}


//--------------------------------------------------------------------
string const&
DOC_Package:: comment( void ) const
//--------------------------------------------------------------------
{
   return my_comment ;
}


//--------------------------------------------------------------------
DOC_Package:: DOC_Package(
   PEL_Object* a_owner,
   std::string const& a_name,
   DOC_Package* a_father,
   std::string const& com ) 
//--------------------------------------------------------------------
      : PEL_Object( a_owner ),
        my_name( a_name ),
        my_comment( com ),
        my_father( a_father ),
        subpackages(PEL_List::create(this)),
        elements(0), element_comments(0)
{
   PEL_LABEL( "DOC_Package:: DOC_Package" ) ;
   
   PEL_CHECK( a_father!=0 ) ;
   
   my_father->subpackages->append( this ) ;
}


//--------------------------------------------------------------------
DOC_Package:: DOC_Package( PEL_Object* a_owner ) 
//--------------------------------------------------------------------
      : PEL_Object( a_owner ),
        my_name( "MASTER" ),
        my_comment( "FrameWork" ),
        my_father( 0 ),
        subpackages(PEL_List::create(this)),
        elements(0), element_comments(0)
{
   // Singleton enforcement
   static bool prem = true ;
   PEL_ASSERT( prem ) ;
   prem = false ;
}


//--------------------------------------------------------------------
DOC_Package:: ~DOC_Package( void ) 
//--------------------------------------------------------------------
{
}


//--------------------------------------------------------------------
PEL_List*
DOC_Package:: packages( void ) 
//--------------------------------------------------------------------
{
   static PEL_List* result = PEL_List::create( PEL_Root::object() );
   return result ;
}


//--------------------------------------------------------------------
DOC_Package*
DOC_Package:: master( void ) 
//--------------------------------------------------------------------
{
   static DOC_Package* result = new DOC_Package( packages() ) ;
   return result ;
}


//--------------------------------------------------------------------
DOC_Package*
DOC_Package:: miscellanous( void ) 
//--------------------------------------------------------------------
{
   static DOC_Package* result = 0 ;
   if( result==0 )
   {
      result = create( "Miscellaneous", master(), "No packaged items" ) ;
   }
   return result ;
}


//--------------------------------------------------------------------
PEL_List const*
DOC_Package:: sub_packages( void ) const
//--------------------------------------------------------------------
{
   return subpackages ;
}


//--------------------------------------------------------------------
void
DOC_Package:: read( std::string const& file ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Package:: read" ) ;
   
   std::ifstream f( file.c_str(), std::ios::in ) ;
   if( !f )
   {
      PEL_Error::object()->raise_file_handling( file.c_str(), "open" ) ;
   }
   getline( f, application_name ) ;
   size_t nbl =0 ;
   while( !f.eof() )
   {
      string rubrique ;
      
      getline( f, rubrique ) ; nbl++ ;
      bool package = ( rubrique.find( "Package inventory" )<rubrique.length() ) ;
      bool classe = ( rubrique.find( "Class inventory" )<rubrique.length() ) ;
      bool white = rubrique[0]==' ' || rubrique[0]=='$' || rubrique[0]=='#' ;
      
      if( !( package || classe || white ) )
      {
         XWarningE( file, nbl, "Bad header in packaging file " + rubrique ) ;
      }
      if( !white )
      {
         bool prem = true ;
         size_t id_name=(size_t)-1, id_comment=(size_t)-1, id_pack=(size_t)-1 ;
      
         while( !f.eof() )
         {
            string nname, ncomment, npack ;
            string buff ;
            getline( f, buff ) ; nbl++ ;
            size_t nb = 0 ;
         
            while( !buff.empty() ) 
            {
               size_t idx = buff.find( ";" ) ;
            
               string item ;
               if( !( idx<buff.length() ) )
               {
                  item = buff ;
                  buff="" ;
               }
               else
               {
                  item = buff.substr( 0, idx ) ;
                  buff = buff.substr( idx+1, buff.length()-idx-1 ) ;
               }
               if( prem )
               {
                  if( item=="Name" ) id_name=nb ;
                  if( item=="Comment" ) id_comment=nb ;
                  if( item=="Package" ) id_pack=nb ;
               }
               else
               {
                  if( nb==id_name ) nname=item ;
                  if( nb==id_comment ) ncomment=item ;
                  if( nb==id_pack ) npack=item ;
               }
               nb++ ;
            }
            
            if( nb<2 ) break ;
            if( !prem )
            {
               
               if( package )
               {
                  if( !npack.empty() && find( npack )==0 )
                  {
                     XWarningE( file, nbl, "Bad package " + npack ) ;
                  }
                  else
                  {
                     DOC_Package* a_father =
                        ( npack.empty() ? master(): find( npack ) ) ;
                     create( nname, a_father, ncomment ) ;
                  }
                  
               }
               
               if( classe )
               {
                  DOC_Package * pack = find( npack ) ;
                  if( pack==0 )
                  {
                     XWarningE( file, nbl, "Bad package " + npack ) ;
                  }
                  else
                  {
                     pack->elements.append( nname ) ;
                     pack->element_comments.append( ncomment ) ;
                  }
               }  
            }
            prem = false ;
         }
      }
   }
}

//---------------------------------------------------------------------------
bool
DOC_Package:: own_class( DOC_Class const* cl ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Package:: own_class" ) ;
   PEL_CHECK_PRE( cl!=0 ) ;
   
   bool result = ( cl->package()==this ) ;

   PEL_Iterator *it = subpackages->create_iterator(0) ;
   for( it->start() ; it->is_valid() && !result ; it->go_next() )
   {
      DOC_Package*p = static_cast<DOC_Package*>( it->item() ) ;
      result = p->own_class( cl ) ;
   }
   it->destroy() ;
   return( result ) ;
}
   
//---------------------------------------------------------------------------
bool
DOC_Package:: own_non_external_classes( void ) const
//---------------------------------------------------------------------------
{
   bool result = false ;

   PEL_Iterator* cl_it = DOC_Class::list()->create_iterator( 0 ) ;
   for( ; cl_it->is_valid() ; cl_it->go_next() )
   {
      DOC_Class* cl = static_cast<DOC_Class* >( cl_it->item() ) ;
      if( own_class( cl  ) && !cl->is_external() )
      {
         result = true ;
         break ;
      }
   }
   cl_it->destroy() ;

   return( result ) ;      
}

//---------------------------------------------------------------------------
std::string const&
DOC_Package:: application( void )
//---------------------------------------------------------------------------
{
   return( application_name ) ;
}



//--------------------------------------------------------------------
void
DOC_Package:: print( std::ostream& os, size_t ident_witdth ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Package:: print" ) ;
   std::string bl( ident_witdth, ' ' ) ;
   os<<bl<<" Package "<<name()<<std::endl;
   PEL_Iterator *it = subpackages->create_iterator(0) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      DOC_Package*p = static_cast<DOC_Package*>( it->item() ) ;
      p->print( os, ident_witdth+1 ) ;
   }
   it->destroy() ;
}
