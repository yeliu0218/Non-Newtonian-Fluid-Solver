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

#include <DOC_Class.hh>

#include <algorithm>
#include <functional>

#include <DOC_Argument.hh>
#include <DOC_Text.hh>
#include <DOC_Category.hh>
#include <DOC_Text.hh>
#include <DOC_ClassItem.hh>
#include <DOC_Sequence.hh>
#include <DOC_Method.hh>
#include <DOC_Tools.hh>
#include <DOC_Package.hh>
#include <DOC_Type.hh>

#include <PEL_Root.hh>

using std::endl ;


//--------------------------------------------------------------------
PEL_List*
DOC_Class::list( void )
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::liste" ) ;
      
   static PEL_List* result = PEL_List::create( PEL_Root::object() ) ;
   
   PEL_CHECK_POST( result!=0 ) ;
   return result ;
}


//--------------------------------------------------------------------
stringVector DOC_Class::upper_DOC_Classs = stringVector(0) ;
//--------------------------------------------------------------------



//--------------------------------------------------------------------
DOC_Class*
DOC_Class::to_class( void ) 
//--------------------------------------------------------------------
{
   return this ;
}



//--------------------------------------------------------------------
DOC_Class::DOC_Class( PEL_Object* a_owner,
                      std::string const& a_name, 
                      DOC_Text* a_comment,
                      std::string const& a_file,
                      int a_line_nb,
                      DOC_Sequence* arg,
                      DOC_Class * inherit_from ) 
//--------------------------------------------------------------------
      : DOC_Symbol( a_owner ),
        myName( a_name ), 
        mother_class( inherit_from ),
        virtuel( false ),
        abstract( false ),
        file(a_file),
        line_number(a_line_nb),
        modele( false ),
        upper_DOC_Class( false ),
        elements(  ),
        inherited_list( PEL_List::create(this) ),
        pub_src( false ),
        COMPONENT( 0 )

{
   if( arg!=0 )
   {
      elements = arg->list() ;
      arg->set_owner( this ) ;
   }
   else
   {
      elements=PEL_List::create(this) ;
   }
   
   PEL_Iterator* it = elements->create_iterator(this) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      DOC_ClassItem* el = dynamic_cast<DOC_ClassItem* >( it->item() ) ;
      PEL_ASSERT( el!=0 ) ;
      add_element( el ) ;
   }
   
   if( inherit_from!=0 )
   {
      inherit_from->add_descendant( this ) ;
   }
   if( a_comment!=0 )
   {
      strcomment = DOC_Tools::text_comment( a_comment->text() ) ;
      modele = strcomment.find( "FRAMEWORK INSTANTIATION" )<strcomment.length() ;
      pub_src = strcomment.find( "PUBLISHED" )<strcomment.length() ;
   }
   upper_DOC_Class = upper_DOC_Classs.has( myName ) ;
   
}



//--------------------------------------------------------------------
DOC_Class::~DOC_Class( void ) 
//--------------------------------------------------------------------
{
}



//--------------------------------------------------------------------
void
DOC_Class::declare_upper( stringVector const& class_name )
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::declare_upper" ) ;
   upper_DOC_Classs = class_name ;
}



//--------------------------------------------------------------------
DOC_Class *
DOC_Class::search( std::string const& a_name ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::search" ) ;
   
   PEL_Iterator* it = list()->create_iterator( 0 ) ;
   DOC_Class * result = 0 ;
   
   for( it->start() ; it->is_valid();  it->go_next() )
   {
      DOC_Class* cl = static_cast<DOC_Class*>( it->item() ) ;
      if( cl->myName == a_name )
      {
	 result = cl ;
         break ;
      }
   }
   it->destroy() ;
   PEL_CHECK_POST( IMPLIES( result!=0, result->name()==a_name ) ) ;
   return result ;
}



//--------------------------------------------------------------------
DOC_ClassItem const*
DOC_Class::search_and_find( std::string const& expr ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::search_and_find" ) ;
   DOC_ClassItem const* result = 0 ;
   
   size_t deb = expr.find( "::" ) ;
   if( deb<expr.length() )
   {
      string DOC_Class_n = expr.substr( 0, deb ) ;
      DOC_Class const* cl = DOC_Class::search( DOC_Class_n ) ;
      if( cl!=0 )
      {
         string type_n = expr.substr( deb+2, expr.length()-deb-2 ) ;
         result = cl->find( type_n ) ;
         if( result==0 )
         { 
            cl->warning( "No method " + type_n + " found" ) ;
	 }
      }
   }
   
   return result ;
}



//--------------------------------------------------------------------
DOC_Class*
DOC_Class::create( std::string const& a_name, 
                   DOC_Text* a_comment,
                   std::string const& a_file,
                   int a_line_nb,
                   DOC_Sequence* arg,
                   DOC_Class * inherit_from  )
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::create" ) ;
   PEL_CHECK_PRE( search(a_name)==0 ) ;
   
   DOC_Class* result =
      new DOC_Class( list(), a_name, a_comment, a_file, a_line_nb, arg, inherit_from ) ;
   list()->append( result ) ;
   return result ;
}



//--------------------------------------------------------------------
bool
DOC_Class::is_upper_class( void ) const
//--------------------------------------------------------------------
{
   return upper_DOC_Class ;
}



//--------------------------------------------------------------------
bool
DOC_Class::is_external( void ) const
//--------------------------------------------------------------------
{
   return elements->count()==0 ;
}



//--------------------------------------------------------------------
std::string 
DOC_Class::short_description( void ) const 
//--------------------------------------------------------------------
{
   std::string result ;
   DOC_Package::attached_to( this, result ) ;
 
   return result ;
}



//--------------------------------------------------------------------
bool
DOC_Class::is_component( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::is_component" ) ;
   
   if( COMPONENT==0 )
   {
      DOC_Class const* c0 = this ;
      while( c0->mother_class !=0 ) c0 = c0->mother_class ;
      COMPONENT = ( c0->name()=="PEL_Object" ? 1 : -1 ) ;
   }
   return COMPONENT > 0  ;
}



//--------------------------------------------------------------------
void
DOC_Class::add_element( DOC_ClassItem * element  )
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::add_element" ) ;
   PEL_CHECK( element!=0 ) ;
   
   element->attach( this ) ;
   DOC_Method* meth ;
   
   if( ( meth=element->method() ) !=0 )
   {
      if( meth->is_virtual() )
      {
         virtuel = true ;
      } 
      if( meth->is_abstract() )
      {
         abstract = true ;
      } 
   }
}



//--------------------------------------------------------------------
void
DOC_Class::add_descendant( DOC_Class const* element )
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::add_descendant" ) ;
   PEL_CHECK( element!=0 ) ;
   
   inherited_list->append( const_cast<DOC_Class *>(element) ) ;
   inherited_list->sort() ;
   virtuel = true ;
}



//--------------------------------------------------------------------
string 
const&
DOC_Class::name( void  ) const
//--------------------------------------------------------------------
{
   return myName ;
}



//--------------------------------------------------------------------
DOC_Package const*
DOC_Class::package( void  ) const
//--------------------------------------------------------------------
{
   std::string comm_court ;
   DOC_Package const* result =
      DOC_Package::attached_to( this, comm_court ) ;
   return result ;
}



//--------------------------------------------------------------------
void
DOC_Class::display_all( void )
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::display_all" ) ;
   
   list()->sort() ;

   PEL_Iterator* it = list()->create_iterator( 0 ) ;
   
   for( it->start() ; it->is_valid();  it->go_next() )
   {
      DOC_Class* cl = static_cast<DOC_Class*>( it->item() ) ;
      PEL_CHECK( dynamic_cast<DOC_Class*>(cl)!=0 ) ;
      // std::cout << "display_all : " << cl->name() << std::endl ;
      // Check for void class
      if( !cl->is_external() )
      {
         DOC_Writer* red = DOC_Writer::create( 0, cl->name() ) ;
         cl->process_documentation( *red ) ;
         red->destroy() ;
      }
   }
   DOC_Writer* red = DOC_Writer::create( 0, "index" ) ;
   red->finalize() ;
   it->destroy() ;
   
   red->destroy() ;
}



//--------------------------------------------------------------------
DOC_Method*
DOC_Class::find_method( DOC_Method const* meth  ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::find_method" ) ;
   PEL_CHECK_PRE( meth!=0 ) ;
   
   DOC_Method * result = 0 ;
   PEL_Iterator* it = elements->create_iterator(0) ;
   
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      DOC_ClassItem* el = static_cast<DOC_ClassItem* >( it->item() ) ;
      DOC_Method * m = el->method() ;
      
      if( m && meth->is_equal( m ) )
      {
	 result = m ;
	 break ;
      }
   }
   it->destroy() ;
   PEL_CHECK_POST( IMPLIES( result!=0, result->is_equal(meth) ) ) ;
   
   return result ;
}



//--------------------------------------------------------------------
DOC_ClassItem*
DOC_Class::find( std::string const& expr) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::find" ) ;
   DOC_ClassItem* result = 0 ;
   PEL_Iterator* it = elements->create_iterator(0) ;
   
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      DOC_ClassItem* el = static_cast<DOC_ClassItem* >( it->item() ) ;
      if( el->declare(expr) )
      {
	 result = el ;
      }
   }
   it->destroy() ;
   if( result==0 && mother_class!=0 )
   {
      result = mother_class->find( expr ) ;
   }
   PEL_CHECK_POST( IMPLIES( result!=0, result->declare(expr) ) ) ;
   return result ;
}



//--------------------------------------------------------------------
void
DOC_Class::inherit( PEL_List* items  ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::inherit" ) ;
   PEL_CHECK_PRE( items!=0 ) ;
   
   items->copy( elements ) ;
   
   for( DOC_Class const* self = this->mother_class ;
	self!=0  ;
	self = self->mother_class )
   {
      PEL_Iterator *it = self->elements->create_iterator( 0 ) ;
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
         DOC_ClassItem* el = static_cast<DOC_ClassItem* >( it->item() ) ;
         PEL_CHECK( dynamic_cast<DOC_ClassItem* >( el )!=0 ) ;
         
	 if( el->is_inheritable() )
	 {
            DOC_Method * m1 = el->method() ;
            if( m1!=0 )
            {
               DOC_Method* m = 0 ;
               PEL_Iterator *it2 = items->create_iterator( 0 ) ;
               for( it2->start() ; it2->is_valid() ; it2->go_next() )
               {
                  PEL_CHECK( dynamic_cast<DOC_ClassItem* >( it2->item() )!=0 ) ;
                  DOC_Method * m2 = static_cast<DOC_ClassItem* >(it2->item())->method() ;
                  if( m2!=0 && m2->is_equal( m1 ) )
                  {
                     m = m2 ;
                     break ;
                  }
               }
               it2->destroy() ;
               if( m == 0 )
               {                  
                  if( !self->is_upper_class() )
                  {
                     items->append( m1 ) ;
                  }
               }
               else if( m!=m1 )
               {
                  m->overide( m1 ) ;
               }
            }
            
	 }
         else
         {
         }
         
      }
      it->destroy() ;
   }
} 


	 
//--------------------------------------------------------------------
void
DOC_Class::verify( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::verify" ) ;
   
   PEL_Iterator* it = elements->create_iterator(0) ;
   
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      DOC_ClassItem* el = static_cast<DOC_ClassItem* >( it->item() ) ;
      DOC_Method * m = el->method() ;
      
      if( m !=0 && m->def_class()==this )
      {
         m->verify() ;
      }
   }
   it->destroy() ;
   if( is_component() ) 
   { 
      verify_default_implementation() ;
      verify_protected_area() ;
   }
}
      

	 
//--------------------------------------------------------------------
void
DOC_Class::verify_default_implementation( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::verify_default_implementation" ) ;
   
   static DOC_Type * void_t = DOC_Type::create( PEL_Root::object(), "" ) ;
   static DOC_Sequence * liste_vide = DOC_Sequence::create( PEL_Root::object() ) ;
   static DOC_Text * vide = DOC_Text::create( PEL_Root::object(), "" ) ;
   
   string mess = "class " ;
   mess += name() ;
   
   DOC_Method const* constructeur_void =
      DOC_Method::create( name(),             // Name 
                     void_t,           // DOC_Type
                     DOC_Method::Private,  // Protection
                     0,                 // DOC_Category
                     vide,             // Comment
                     liste_vide,       // DOC_Arguments
                     liste_vide ) ;    // Modifiers

   DOC_Method const* meth = find_method( constructeur_void ) ;
   if( meth==0 )
   {
      warning( mess + " doesn't declare its non public default constructor" ) ;
   }
   else
   {
      if( !meth->comment().empty() )
         meth->warning( DOC_Method::Definition, 
                        "constructor " + meth->signature() + 
                        " does not needs any comment" ) ;
   }
   constructeur_void->destroy() ;
   
   DOC_Type* const_et = DOC_Type::create( 0, name() ) ;
   const_et->set_constant() ;
   const_et->set_reference() ;
   DOC_Argument* arg_const_et = DOC_Argument::create( 0, const_et, "" ) ;
   DOC_Sequence* liste_arg_const_et = DOC_Sequence::create( 0 ) ;
   liste_arg_const_et->list()->append( arg_const_et ) ;
   
   DOC_Method const* copy_constructeur =
      DOC_Method
      ::create( name(),             // Name 
                       void_t,           // DOC_Type
                       DOC_Method::Private,  // Protection
                       0,                 // DOC_Category
                       vide,             // Comment
                       liste_arg_const_et,     // DOC_Arguments
                       liste_vide ) ;    // Modifiers
   
   meth = find_method( copy_constructeur ) ;
   if( meth==0 )
   {
      warning( mess + " doesn't declare private copy constructor" ) ;
   }
   else
   {
      if( !meth->comment().empty() )
         meth->warning( DOC_Method::Definition, 
                        meth->signature() + " does not needs any comment" ) ;
   }
   copy_constructeur->destroy() ;
   
   DOC_Method const* operateur_affectation =
      DOC_Method::create( "operator=",       // Name 
                       const_et,         // DOC_Type
                       DOC_Method::Private,  // Protection
                       0,                 // DOC_Category
                       vide,             // Comment
                       liste_arg_const_et,     // DOC_Arguments
                       liste_vide ) ;    // Modifiers
   meth = find_method( operateur_affectation ) ;
   if( meth==0 )
   {
      warning( mess + " doesn't declare private operator=" ) ;
   }
   else
   {
      if( !meth->comment().empty() )
         meth->warning( DOC_Method::Definition, 
                        meth->signature() + " does not needs any comment" ) ;
   }
   operateur_affectation->destroy() ;

   bool a_destructeur = false ;
   PEL_Iterator* it = elements->create_iterator(0) ;
   
   for( it->start() ; !a_destructeur && it->is_valid() ; it->go_next() )
   {
      DOC_ClassItem* el = static_cast<DOC_ClassItem* >( it->item() ) ;
       a_destructeur = a_destructeur ||
         ( el->method() && el->method()->is_destructor() ) ;
      if( a_destructeur && !el->comment().empty() )
      {
         el->method()->warning( DOC_Method::Definition, 
                                      el->method()->signature() + 
                                      " does not needs any comment" ) ;
      }

   }
   it->destroy() ;
   if( !a_destructeur )
   {
      warning( mess + " doesn't declare destructor" ) ;
   }              
   const_et->destroy() ;
   arg_const_et->destroy() ;
   liste_arg_const_et->destroy() ;
}



//--------------------------------------------------------------------
void
DOC_Class:: verify_protected_area( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::verify_protected_area" ) ;
   
   bool has_protected_destuctor = false ;
   bool has_protected_constructor = false ;
   bool has_invariant = false ;
   bool invariant_is_protected = false ;
   bool has_protected_method = false ;

   PEL_Iterator* it = elements->create_iterator( 0 ) ;
   
   for( it->start() ; it->is_valid();  it->go_next() )
   {
      DOC_ClassItem* el = static_cast<DOC_ClassItem* >( it->item() ) ;
      DOC_Method* meth = el->method() ;
      if( meth!=0 && meth->def_class()==this )
      {
         if( meth->is_destructor() )
	 {
            if( meth->protection() == DOC_ClassItem::Public )
               meth->warning( DOC_Method::Definition, 
                              meth->signature() + " should not be public" ) ;
            if( meth->protection() == DOC_ClassItem::Protected )
            {
               has_protected_destuctor = true ;
               if( !meth->is_virtual() )
                  meth->warning( DOC_Method::Definition, 
                              "a protected destructor should be virtual" ) ;
	    }
            if( meth->protection() == DOC_ClassItem::Private && 
                meth->is_virtual() )
	    {
               meth->warning( DOC_Method::Definition, 
                              "a private destructor should not be virtual" ) ;
	    }
                                                 
	 }
         if( meth->is_constructor() )
         {
            if( meth->protection() == DOC_ClassItem::Public )
               meth->warning( DOC_Method::Definition, 
                              meth->signature() + " should not be public" ) ;
            if( meth->protection() == DOC_ClassItem::Protected )
               has_protected_constructor = true ;
         }
         if( meth->is_invariant() )
         {
            has_invariant = true ;
            if( meth->protection()==DOC_ClassItem::Protected ) 
               invariant_is_protected = true ;
	 }
         if( meth->protection()==DOC_ClassItem::Protected )
            has_protected_method = true ;
      }
   }
   it->destroy() ;
   if( has_protected_method )
   {
      bool consistent = false ;
      if( has_protected_destuctor && has_protected_constructor )
      {
         if( has_invariant && invariant_is_protected )
            consistent = true ;
         else if( !has_invariant )
            consistent = true ;
      }
      if( !consistent )
      {
         warning( name() + " has an inconsistent protected area" ) ;
      }
   }
}



//--------------------------------------------------------------------
void
DOC_Class:: warning( std::string const& mess ) const
//--------------------------------------------------------------------
{
   XWarningE( file, line_number , mess ) ;
}



//--------------------------------------------------------------------
bool
DOC_Class::has_to_document( DOC_ClassItem const* item ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::has_to_document" ) ;
   PEL_CHECK_PRE( item!=0 ) ;
   
   return( item->has_to_document(this) ) ;
}



//--------------------------------------------------------------------
bool
DOC_Class::is_virtual( void ) const
//--------------------------------------------------------------------
{
   return( virtuel ) ;
}



//--------------------------------------------------------------------
DOC_Class const*
DOC_Class::mother( void ) const
//--------------------------------------------------------------------
{
   return mother_class ;
}



//--------------------------------------------------------------------
std::string const& 
DOC_Class::comment( void ) const
//--------------------------------------------------------------------
{
   return strcomment ;
}



//--------------------------------------------------------------------
PEL_List const*
DOC_Class::inherited_classes( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::inherited_classes" ) ;
   PEL_List const* result = inherited_list ;
   PEL_CHECK_POST( result!=0 ) ;
   return result ;
}



//--------------------------------------------------------------------
PEL_List const*
DOC_Class::owned_components( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::owned_components" ) ;
   PEL_List const* result = elements ;
   PEL_CHECK_POST( result!=0 ) ;
   return result ;
}


//--------------------------------------------------------------------
void
DOC_Class::process_documentation( DOC_Writer& sullizer )
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class::process_documentation" ) ;

   /* Recuperation des methodes heritees */
   PEL_List* list_inherite = PEL_List::create( 0 ) ;
   inherit( list_inherite ) ;
   list_inherite->sort() ;
   elements->copy( list_inherite ) ;
   list_inherite->destroy() ;
   
   /* Recuperation des pre-conditions et post-conditions a travers les
      appels de fonction. */
   PEL_Iterator* it = elements->create_iterator( 0 ) ;
   stringVector given(0) ;
   
   for( it->start() ; it->is_valid();  it->go_next() )
   {
      DOC_ClassItem* el = static_cast<DOC_ClassItem* >( it->item() ) ;
      
      DOC_Method* meth ;
      if( (meth=el->method()) !=0 )
      {
         /* Recuperation des pre-conditions */
         meth->recover_conditions( this, DOC_Method::pre ) ;
         /* Recuperation des post-conditions */
         meth->recover_conditions( this, DOC_Method::post ) ;
      }

      if( el->bookmark().empty() )
      {
         std::string a_name = el->name() ;
         while( given.has( a_name ) )
         {
            a_name += "_" ;
         }
         given.append( a_name ) ;
         el->set_bookmark( a_name ) ;
      }      
   }
   it->destroy() ;

   // Processus de verification
   verify() ;

   // Ecriture
   sullizer.process_documentation( this ) ;
  
}




//--------------------------------------------------------------------
bool
DOC_Class:: publish_source( void ) const 
//--------------------------------------------------------------------
{
   return pub_src ;
}



//--------------------------------------------------------------------
bool
DOC_Class:: is_plug_point( void ) const 
//--------------------------------------------------------------------
{
   return modele ;
}



//--------------------------------------------------------------------
int
DOC_Class:: three_way_comparison( PEL_Object const* other ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Class:: three_way_comparison" ) ;
   PEL_CHECK_PRE( dynamic_cast<DOC_Class const*>( other ) != 0 ) ;
   DOC_Class const* class2 = static_cast<DOC_Class const*>( other ) ;
   
   int result = 0 ;
   
   std::string const& n1 = name() ;
   std::string const& n2 = class2->name() ;
   if( n1<n2 )
   {
      result=-1 ;
   }
   else if( n1>n2 )
   {
      result=1 ;
   }
   return result ;
}

//--------------------------------------------------------------------
PEL_List const*
DOC_Class::all_elements( void  ) const
//--------------------------------------------------------------------
{
   return elements ;
}

