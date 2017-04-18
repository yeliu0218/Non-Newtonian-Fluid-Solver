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

#include <DOC_Category.hh>
#include <DOC_Tools.hh>

#include <PEL_List.hh>
#include <PEL_Root.hh>

#include <sstream>

size_t DOC_Category::global_index = 0 ;
double const default_rank = 1000.0 ;

//--------------------------------------------------------------------
PEL_List*
DOC_Category::categories( void ) 
//--------------------------------------------------------------------
{
   static PEL_List* result = PEL_List::create( PEL_Root::object() ) ;
   return result ;
}


//--------------------------------------------------------------------
DOC_Category const*
DOC_Category:: create( std::string const& a_name ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Category:: create" ) ;
   
   std::string the_name = DOC_Tools::text_comment( a_name ) ;
   size_t ord = the_name.find_last_of("(") ;
   double rank = default_rank ;
   bool has_rank = ord<the_name.length() && the_name[the_name.length()-1]==')';
   if( has_rank )
   {
      std::string sub = the_name.substr( ord+1, the_name.length()-ord-2 ) ;
      std::istringstream is( sub ) ;
      is >> rank ;
      if( is.fail() )
         rank=default_rank ;
      else
         the_name = the_name.substr( 0, ord ) ;
   }
   DOC_Category * result = 0 ;
   PEL_Iterator* it = categories()->create_iterator(0) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      DOC_Category* cat = static_cast<DOC_Category*>( it->item() ) ;
      if( cat->name() == the_name )
      {
         result = cat ;
      }
   }
   if( result==0 )
   {
      result = new DOC_Category( categories(), the_name, rank ) ;
      categories()->append( result ) ;
   }
   else
   {
      double old_rank = result->rank ;
      if( rank!=default_rank && old_rank!=rank )
      {
         if( old_rank==default_rank )
         {
            result->rank = rank ;
         }
         else
         {
            std::ostringstream mesg ;
            mesg << "the rank of category \"" << the_name << "\" " ;
            mesg << "was previoulsy set at " << old_rank ;
            mesg << " instead of " << rank ;
            XWarningE( DOC_Tools::file(), 0 , mesg.str() ) ;
         }
      }
   }
   
   it->destroy() ;
   
   PEL_CHECK_POST( FORALL( ( size_t i=0 ; i<global_index ; i++ ),
                           dynamic_cast<DOC_Category * >( categories()->at( i ) )->index()==i ) ) ;
   
   return result ;
}



//--------------------------------------------------------------------
std::string const&
DOC_Category:: name( void ) const
//--------------------------------------------------------------------
{
   return my_name ;
}


//--------------------------------------------------------------------
int
DOC_Category:: index( void ) const
//--------------------------------------------------------------------
{
   return my_index ;
}


//--------------------------------------------------------------------
DOC_Category:: DOC_Category( PEL_Object* a_owner,
                             std::string const& a_name,
                             double a_rank ) 
//--------------------------------------------------------------------
      : PEL_Object( a_owner ),
        my_index( global_index++ ),
        my_name( a_name ),
        rank( a_rank )
{
}


//--------------------------------------------------------------------
DOC_Category:: ~DOC_Category( void ) 
//--------------------------------------------------------------------
{
}


//--------------------------------------------------------------------
void
DOC_Category:: display_list( std::ostream& os, size_t indent_width ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Category:: display_list" ) ;
   categories()->sort() ;
   std::string bl( indent_width, ' ' ) ;
   os <<bl<< "List of categories : "<<std::endl;
   PEL_Iterator* it = categories()->create_iterator(0) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      DOC_Category* cat = static_cast<DOC_Category*>( it->item() ) ;
      os <<bl<< "-"<<cat->name();
      if( cat->rank!=default_rank ) os<<" ("<<cat->rank<<")";
      os<<std::endl;
   }
   it->destroy() ;
}


//--------------------------------------------------------------------
int
DOC_Category:: three_way_comparison( PEL_Object const* other ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Category:: three_way_comparison" ) ;
   PEL_CHECK_PRE( three_way_comparison_PRE( other ) ) ;
   DOC_Category const* cat = static_cast<DOC_Category const*>( other ) ;
   
   int result = 0 ;
   if( rank < cat->rank )
   {
      result = -1 ;
   }
   else if( rank > cat->rank )
   {
      result = 1 ;
   }
   else
   {
      if( index()<cat->index() )
      {
         result = -1 ;
      }
      else if( index()>cat->index() )
      {
         result = 1 ;
      }
   }
   PEL_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return result ;
}


