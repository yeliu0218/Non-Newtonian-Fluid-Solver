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

#include <DOC_Method.hh>

#include <PEL.hh>
#include <DOC_Argument.hh>
#include <DOC_Category.hh>
#include <DOC_Text.hh>
#include <DOC_Class.hh>
#include <DOC_Function.hh>
#include <DOC_Sequence.hh>
#include <DOC_Tools.hh>
#include <DOC_Type.hh>
#include <DOC_Typedef.hh>
#include <string>
#include <algorithm>
#include <sstream>

using std::cerr ;
using std::endl ;
using std::istringstream ;

//--------------------------------------------------------------------
DOC_Method*
DOC_Method::create( std::string const& a_name,
                DOC_Type * a_type,
                Protection protection,
                DOC_Category const* category,
                DOC_Text * comment,
                DOC_Sequence* arguments,
                DOC_Sequence* modifiers,
                int line,
                std::string const& file ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::create" ) ;
   
   PEL_CHECK_PRE( a_type!=0 ) ;
   PEL_CHECK_PRE( arguments!=0 ) ;
   PEL_CHECK_PRE( modifiers!=0 ) ;
   
   DOC_Method * result = 0 ;
   size_t del = a_name.find("::") ;
   if( del>= a_name.length() )
   {
      /* Definition */
      result = new DOC_Method( a_name,
                               a_type,
                               protection,
                               category,
                               comment,
                               arguments,
                               modifiers ) ;
      result->def_line = line ;
      result->def_fich = file ;
   }
   else
   {
      /* Implementation */
      string class_name = a_name.substr( 0, del ) ;
      string name_DOC_Method = a_name.substr( del+2, a_name.length()-del-2 ) ;
      DOC_Class* aclass = DOC_Class::search( class_name ) ;
      if( aclass!=0 )
      {
	 DOC_Method * tmp = new DOC_Method( name_DOC_Method,
                                            a_type,
                                            protection,
                                            category,
                                            comment,
                                            arguments,
                                            modifiers ) ;
         tmp->attach( aclass ) ;
	 result = aclass->find_method( tmp ) ;
         // Attention : il faut differencier une declaration friend de l'implementation
         if( result != 0 && DOC_Tools::is_implementation_file() )
         {
            if( result->implementation!=0 )
            {
               XWarningE( file, line, result->full_name() + " is already implemented " ) ;
               XWarningE( result->imp_fich, result->imp_line, " here ! " ) ;
            }
            result->implementation = tmp ;
            result->imp_line = line ;
            result->imp_fich = file ;
            tmp->set_owner( aclass ) ;
         }
         else
         {
            delete tmp ;
         }
      }
      if( result==0 )
      {
         string mess = "Method " ;
         mess += a_name ;
         mess += " hasn't been yet declared" ;
         XWarningE( file, line, mess ) ;
         XWarningE( file, line, "or has been declared with a different prototype." ) ;
      }
   }
   
   return result ;
}



//--------------------------------------------------------------------
DOC_Method::DOC_Method( std::string const& a_name,
                DOC_Type * a_type,
                Protection protection,
                DOC_Category const* category,
                DOC_Text * comment,
                DOC_Sequence* arguments,
                DOC_Sequence* modifiers ) 
//--------------------------------------------------------------------
   : DOC_Attribute( a_name, a_type, protection, category, comment ),
     args(PEL_List::create(this)),
     estVirtuel( false ), estVirtuel_pur( false ), is_constant( false ),
     surdefinition( 0 ),
     fait_PRE( false ),
     fait_POST( false ),
     constructeur( false ),
     destructeur( false ),
     implementation( 0 ),
     def_line( -1 ),
     imp_line( -1 ),
     preConditions(PEL_List::create(this)),
     postConditions(PEL_List::create(this)),
     proceded_preConditions(PEL_List::create(this)),
     proceded_postConditions(PEL_List::create(this)),
     referring_item( 0 )
{
   PEL_LABEL( "DOC_Method::DOC_Method" ) ;
   
   bool constructeur_destructeur = a_type->full_type_name().empty() ;
   if( constructeur_destructeur )
   {
      if( a_name[0] == '~' )
      {
         destructeur = true ;
      }
      else
      {
         constructeur = true ;
      }
   }
   PEL_Iterator* it = arguments->list()->create_iterator(0) ;
   for( it->start() ;  it->is_valid() ;  it->go_next()  )
   {
      args->append( it->item() ) ;
   }
   it->destroy() ;
   
   it = modifiers->list()->create_iterator(0) ;
   for( it->start() ;  it->is_valid() ;  it->go_next()  )
   {
      DOC_Symbol*s = static_cast<DOC_Symbol*>( it->item() ) ;
      string const& str = s->to_text()->text() ;
      if( str=="const" )
      {
	 is_constant = true ;
      }
      else if( str=="abstract" )
      {
	 estVirtuel_pur = true ;
      }
      else
      {
         Erreur( "Modificateur inconnu : " ) ;
      }
   }
   it->destroy() ;
}



//--------------------------------------------------------------------
DOC_Method::~DOC_Method( void ) 
//--------------------------------------------------------------------
{
}



//--------------------------------------------------------------------
DOC_Method*
DOC_Method::method( void )
//--------------------------------------------------------------------
{
   return this ;
}



//--------------------------------------------------------------------
void
DOC_Method::set_constant( void ) 
//--------------------------------------------------------------------
{
   is_constant = true ;
}



//--------------------------------------------------------------------
void
DOC_Method::set_virtual( void ) 
//--------------------------------------------------------------------
{
   estVirtuel = true ;
}



//--------------------------------------------------------------------
void
DOC_Method::overide( DOC_Method const* precedent ) 
//--------------------------------------------------------------------
{
   surdefinition = precedent ;
//   if( category()->name().empty() ) 
//   {
//      inherit_category( precedent ) ;
//   }
}



//--------------------------------------------------------------------
bool
DOC_Method::is_virtual( void ) const
//--------------------------------------------------------------------
{
   return estVirtuel ;
}



//--------------------------------------------------------------------
bool
DOC_Method::is_abstract( void ) const
//--------------------------------------------------------------------
{
   return estVirtuel_pur ;
}



//--------------------------------------------------------------------
void
DOC_Method::set_label( std::string const& a_label )
//--------------------------------------------------------------------
{
   my_label = a_label ;
}



//--------------------------------------------------------------------
std::string const&
DOC_Method::label( void ) const
//--------------------------------------------------------------------
{
   return my_label ;
}



//--------------------------------------------------------------------
string 
DOC_Method::signature( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::signature" ) ;
   string ret ;
   if( estVirtuel )
   {
      ret += "virtual " ;
   }
   ret += DOC_Attribute::signature( ) +"( ";
   bool prem = true ;
   if( args->count()>0 )
   {
      PEL_Iterator*it = args->create_iterator(0) ;
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
	 if( !prem )
	 {
	    ret += ", " ;
	 }
	 prem = false ;
         DOC_Argument*a=static_cast<DOC_Argument*>( it->item() ) ;
         PEL_CHECK( dynamic_cast<DOC_Argument*>(a)!=0) ;
	 ret += a->type()->full_type_name() ;
         ret += " " ;
	 ret += a->text() ;
      }
      it->destroy() ;
   }
   else
   {
      ret+="void";
   }
   ret += " )" ;
   if( is_constant )
   {
      ret += " const" ;
   }
   if( estVirtuel_pur )
   {
      ret += " = 0" ;
   }
   return ret ;
}



//--------------------------------------------------------------------
string 
DOC_Method::prototype( DOC_Writer& sullizer ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::prototype" ) ;
   size_t nb_bl = DOC_Attribute::signature().length() + 2 ;
   string ret ;
   if( estVirtuel )
   {
      ret += "virtual " ;
      nb_bl += 8 ;
   }
   ret += DOC_Attribute::prototype( sullizer ) ;
   ret += "( " ;
   
   bool prem = true ;
   if( args->count()>0 )
   {
      PEL_Iterator*it = args->create_iterator(0) ;
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
         DOC_Argument*a=static_cast<DOC_Argument*>( it->item() ) ;
	 if( !prem )
	 {
	    ret += ",\n" + string( nb_bl, ' ' ) ;
	 }
	 prem = false ;
         DOC_ClassItem* e = 0 ;
         if( ( e = def_class()->find(
                a->type()->name() ) )!=0 )
         {
            ret += sullizer.reference( e->def_class(),
                                       e->name(),
                                       a->type()->full_type_name() ) ;
         }
         else
         {
            ret += a->type()->type_reference( sullizer ) ;
         }
         ret += " " ;
         ret += sullizer.underlined_text( 1, a->text() ) ;
      }
      it->destroy() ;
   }
   else
   {
      ret += "void" ;
   }
   ret += " ) " ;
   if( is_constant )
   {
      ret += "const" ;
   }
   if( estVirtuel_pur )
   {
      ret += " = 0" ;
   }
   return ret ;
}


//--------------------------------------------------------------------
bool
DOC_Method::has_argument( std::string const& argstr ) const 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::has_argument" ) ;
   bool ret = false ;
   
   PEL_Iterator*it = args->create_iterator(0) ;
   for( it->start() ; !ret && it->is_valid() ; it->go_next() )
   {
      DOC_Argument* arg=static_cast<DOC_Argument* >(it->item()) ;
      ret = arg->silent_var()==argstr ;
   }
   it->destroy() ;
   return ret ;
}



//--------------------------------------------------------------------
void
DOC_Method::display( std::ostream& out ) const 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::display" ) ;
   if( def_class() )
   {
      out << "[" << def_class()->name() << "]" ;
   }
   
   if( !category()->name().empty() )
   {
      out << " Method category : " << category()->name() << endl ;
   }
   DOC_Attribute::display( out ) ;
   if( estVirtuel_pur )
   {
      out << " ABSTRACT " ;
   }
   else if( estVirtuel )
   {
      out << " VIRTUAL " ;
   }
   out << signature() ;
   if( preConditions->count()>0 )
   {
      PEL_Iterator*it = preConditions->create_iterator(0) ;
      out  << "Pre-conditions : " << endl ;
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
         DOC_Function*f = static_cast<DOC_Function*>( it->item() ) ;
	 out << f->text() << endl ;
      }
      it->destroy() ;
   }
   if( postConditions->count()>0 )
   {
      PEL_Iterator*it = postConditions->create_iterator(0) ;
      out << "Post-conditions : " << endl ;
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
         DOC_Function*f = static_cast<DOC_Function*>( it->item() ) ;
	 out << f->text() << endl ;
      }
      it->destroy() ;
   }

}



//--------------------------------------------------------------------
std::string
DOC_Method:: format_comment( DOC_Writer const& red,
                         size_t marqueur,
                         std::string const& com,
                         char debut,
                         char fin,
                         DOC_Class const* def ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::format_comment" ) ;
   string ret = com ;
   size_t f, e ;
   while( ( f = ret.find( debut ) ) < ret.length() )
   {
      string tmp = ret.substr( f+1, ret.length()-f ) ;
      e = tmp.find( fin ) ;
      if( e < tmp.length() )
      {
         string mot = ret.substr( f+1, e ) ;
         string mot_doc = mot ;
         size_t debpp = mot.find( '(' ) ;
         if( debpp<mot.length() ) mot = mot.substr( 0, debpp ) ;
         debpp = mot.find( "::" ) ;
         string cible = mot ;
         DOC_Class const* cl = 0 ;
         DOC_Class const* my_def = def ;
         if( debpp<mot.length() )
         {
            cible = mot.substr( debpp+2, mot.length()-debpp-2 ) ;
            cl = DOC_Class::search( mot.substr( 0, debpp ) ) ;

            if( cible.empty() ) mot_doc=mot.substr( 0, debpp ) ;
            
            if( cl!=0 ) my_def = cl ;
            else mot_doc=cible ;
         }
            
         if( cl!=0 || my_def->find( cible ) !=0 )
         {
            ret.replace( f, e+2, red.reference( my_def,
                                                cible,
                                                mot_doc ) ) ;
         }
         else
         {
            ret.replace( f, e+2, red.underlined_text( 1, mot_doc ) ) ;
         }
      }
      else
      {
         Erreur( "Bad comment : " << com << " : no termination symbol " ) ;
         break ;
      }
   }
   return ret ;
}






//--------------------------------------------------------------------
bool 
DOC_Method::is_inheritable( void ) const 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::is_inheritable" ) ;
   bool ret = DOC_Attribute::is_inheritable() && !is_static() && !constructeur && !destructeur ;
   return ret ;
}


//--------------------------------------------------------------------
DOC_Method const* 
DOC_Method::overiden_method( void ) const 
//--------------------------------------------------------------------
{
   return surdefinition ;
}


//--------------------------------------------------------------------
PEL_List const*
DOC_Method::pre_conditions( void ) const 
//--------------------------------------------------------------------
{
   return preConditions ;
}



//--------------------------------------------------------------------
PEL_List const*
DOC_Method::post_conditions( void ) const 
//--------------------------------------------------------------------
{
   return postConditions ;
}




//--------------------------------------------------------------------
bool 
DOC_Method::is_constructor( void ) const 
//--------------------------------------------------------------------
{
   return constructeur ;
}



//--------------------------------------------------------------------
bool 
DOC_Method::is_destructor( void ) const 
//--------------------------------------------------------------------
{
   return destructeur ;
}



//--------------------------------------------------------------------
bool 
DOC_Method::is_precondition( void ) const  
//--------------------------------------------------------------------
{
   return name().find( "_PRE" ) < name().length() ;
}



//--------------------------------------------------------------------
bool 
DOC_Method::is_invariant( void ) const  
//--------------------------------------------------------------------
{
   return name().find( "invariant" ) < name().length() ;
}



//--------------------------------------------------------------------
bool 
DOC_Method::is_postcondition( void ) const  
//--------------------------------------------------------------------
{
   return name().find( "_POST" ) < name().length() ;
}



//--------------------------------------------------------------------
string 
DOC_Method::short_name( void ) const  
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::short_name" ) ;
   string ret = name() ;
   size_t idx = ret.find_last_of( '_' ) ;
   if( idx<ret.length() ) ret = ret.substr( 0, idx ) ;
   
   return ret ;
}



//--------------------------------------------------------------------
void 
DOC_Method::add_condition( DOC_Function * test, Condition cond )  
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::add_condition" ) ;
   
   bool precond = cond==pre ;
   bool postcond = cond==post ;
   bool assertcond = cond==assert ;
   PEL_CHECK( precond || postcond || assertcond ) ;
   
   if( precond || ( is_precondition() && assertcond ) )
   {
      preConditions->append( test ) ;
   }
    if( postcond || ( is_postcondition() && assertcond ) )
   {
      postConditions->append( test ) ;
   }
  
    if( is_precondition() && !assertcond )
    {
       XWarningE( DOC_Tools::file(),
                  DOC_Tools::current_line_number(),
                  "Precondition not allowed in precondition implementation method" ) ;
    }
    if( is_postcondition() && !assertcond )
    {
       XWarningE( DOC_Tools::file(),
                  DOC_Tools::current_line_number(),
                  "Postcondition not allowed in postcondition implementation method" ) ;
    }

}



//--------------------------------------------------------------------
string 
DOC_Method::comment( void ) const 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::comment" ) ;
   string ret = DOC_Attribute::comment() ;
   if( surdefinition!=0 )
   {
      string ret1 = surdefinition->comment() ;
      size_t idx = ret1.find( "IMPLEMENTATION" ) ;
      if( idx<ret1.length() ) ret1 = ret1.substr( 0, idx ) ;
      if( !ret1.empty() )
      {
	 ret = ret1 + "\n" +  ret ;
      }
   }
   return ret ;
   
}



//--------------------------------------------------------------------
bool
DOC_Method::has_compatible_arguments( DOC_Method const* called ) const 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::has_compatible_arguments" ) ;
   bool result = true ;
   
   PEL_Iterator*it = args->create_iterator(0) ;
   PEL_Iterator*ita = called->args->create_iterator(0) ;
   for( it->start() ; result && it->is_valid() ; it->go_next() )
   {
      DOC_Argument* arg=static_cast<DOC_Argument* >(it->item()) ;
      bool found = false ;
      for( ita->start() ; !found && ita->is_valid() ; ita->go_next() )
      {
         DOC_Argument* aarg=static_cast<DOC_Argument* >(ita->item()) ;
         if( arg->silent_var() == aarg->silent_var() &&
             arg->type()->name() == aarg->type()->name() )
         {
            found = true ;
         }
      }
      
      result= found ;
   }
   it->destroy() ;
   ita->destroy() ;
   return result ;
}



//--------------------------------------------------------------------
DOC_Method*
DOC_Method::has_method( DOC_Class const* a_classe,
                        std::string const& method_name,
                        std::string const& ext ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::has_method" ) ;
   DOC_Method* ret = 0 ;
   std::string m_name = method_name ;

   size_t idx = m_name.find( "::" ) ;
   if( ext.empty() || m_name.find( ext ) == m_name.length()-ext.length() )
   {
      if( idx < m_name.length() )
      {
         string nameDOC_Class = m_name.substr( 0, idx ) ;
         m_name = m_name.substr( idx+2, m_name.length()-idx-2 ) ;
         a_classe = DOC_Class::search( nameDOC_Class ) ;
      }
      
      while( ret==0 && a_classe!=0 )
      {
         PEL_List* choices = PEL_List::create( 0 ) ;
         PEL_Iterator * it = a_classe->all_elements()->create_iterator( 0 ) ;
         for( it->start() ; it->is_valid() ; it->go_next() ) 
         {
            DOC_ClassItem* e = static_cast<DOC_ClassItem*>( it->item() ) ;
            DOC_Method* m = e->method() ;
            if( m!=0 && (m->def_class()==a_classe) && m->name()==m_name )
            {
               choices->append( m ) ;
               ret = m ;
            }
            
         }
               
         it->destroy() ;
         if( choices->count()>1 )
         {  // au moins 2 choix : il faut departager
            ret = 0 ;
            it = choices->create_iterator( 0 ) ;
            for( it->start() ; it->is_valid() ; it->go_next() ) 
            {
               DOC_Method* m = static_cast<DOC_Method*>( it->item() ) ;
               bool good_m = has_compatible_arguments( m ) ;
               size_t nb_args = args->count() ;
               
               if( m->is_postcondition() && type()->name()!="void" )
               {
                   nb_args++ ;
               }
               good_m = good_m && m->args->count() == nb_args ;
                   
               if( good_m ) 
               {
                  if( ret!=0 ) 
                  {
                     warning( ( implementation!=0 ? Implementation : Definition ), 
                              "The method " + full_name() +
                              " seems to refer to several methods : " ) ;
                     ret->warning( Definition, "First one is here" ) ;
                     m->warning( Definition, "Second one is there" ) ;
                  }
                  else
                  {
                     ret = m ;
                  }
               }
            }
            it->destroy() ;
          }
         choices->destroy() ;
            
         a_classe = a_classe->mother() ;
      }
   }
   return ret ;
}



//--------------------------------------------------------------------
void
DOC_Method::recover_conditions( DOC_Class * aclass,
                                Condition preCondition ) 
//--------------------------------------------------------------------
{
//    cout << "Recuperation des conditions de " << def_class()->name()
//         << "::" << name()  << endl ;
   PEL_LABEL( "DOC_Method::recover_conditions" ) ;
   
   PEL_List* lst ;
   bool * fait ;
   string ext ;
   PEL_List* newList ;
   
   if( preCondition==pre )
   {
      lst = preConditions ;
      newList = proceded_preConditions ;
      fait = &fait_PRE ;
      ext = "_PRE" ;
   }
   else
   {
      lst = postConditions ;
      newList = proceded_postConditions ;
      fait = &fait_POST ;
      ext = "_POST" ;
   }
   if( *fait ) newList->clear() ;
   
//   if( ! (*fait) )
   {
      if( is_abstract() && implementation== 0 )
      {
         DOC_Method* meth = has_method( aclass, name()+ext, "" ) ;
//         PEL::out() << "Search for "<<name()+ext<<" in "<<aclass->name()<<" = "<<meth<<std::endl;
         
         if( meth!=0 )
         {
            if( meth->def_class()==def_class() )
               meth->referring_item = this ;
            meth->recover_conditions( aclass, preCondition ) ;
            newList->copy( meth->conditions_list( preCondition ) ) ;
         }
      }
      else
      {
         PEL_Iterator* it=lst->create_iterator(0) ;
         for( it->start() ; it->is_valid() ; it->go_next() )
         {
            DOC_Function* f = static_cast<DOC_Function* >( it->item() ) ;
            
            DOC_Method* meth = has_method( aclass, f->name(), ext ) ;

            if( meth==0 &&
                ( f->name().find( "_PRE" ) < f->name().length()
                  || f->name().find( "_POST" ) < f->name().length() ) &&
                DOC_Tools::warn_unreachable_condition() )
            {
               std::string cond = ( preCondition==pre ? "precondition " : "postcondition " ) ;
               
               warning( Implementation, 
                        "The " + cond
                        + f->name() + " can't be found" ) ;
            }
            
            
            if( meth==0 || ! ( meth->is_precondition() || meth->is_postcondition() ) )
            {
               newList->append( it->item() ) ;
            }
            else if( meth==this ) 
            {
               cerr << "Find recursivity when searching for condition in method : " << endl ;
               meth->display( cerr ) ;
               Erreur( "Bad file" ) ;
            }
            else
            {
               if( meth->def_class()==def_class() )
                  meth->referring_item = this ;
               meth->recover_conditions( aclass, preCondition ) ; 
               PEL_List const* functionList = meth->conditions_list( preCondition ) ; 
               PEL_Iterator* itf=functionList->create_iterator(0) ;
               for( itf->start() ; itf->is_valid() ; itf->go_next() )
               {
                  DOC_Function* ff = static_cast<DOC_Function* >( itf->item() ) ;
                  newList->append( ff ) ;
               }
               itf->destroy() ;
            }
         }
         it->destroy() ;
      }
      *fait = true ;
   }
}


//--------------------------------------------------------------------
PEL_List const*
DOC_Method::conditions_list( Condition preCondition ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::conditions_list" ) ;
   
   PEL_List* result ;
   
   if( preCondition==pre )
   {
      PEL_CHECK( fait_PRE ) ;
      result = proceded_preConditions ;
   }
   else
   {
      PEL_CHECK( fait_POST ) ;
      result = proceded_postConditions ;
   }

   return result ;
}




//--------------------------------------------------------------------
string
DOC_Method::full_name( void ) const
//--------------------------------------------------------------------
{
   return def_class()->name() + "::" +
      signature() ;
}



//--------------------------------------------------------------------
std::string const&
DOC_Method::filename( File fich ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::filename" ) ;
   
   return ( fich==Implementation ? imp_fich : def_fich ) ;
}


//--------------------------------------------------------------------
DOC_Method const*
DOC_Method::referring_method( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::referring_method" ) ;
   DOC_Method const* result = 0 ;
   if( referring_item!=0 )
      result = const_cast<DOC_ClassItem*>(referring_item)->method() ;
   
   return ( result ) ;
}


//--------------------------------------------------------------------
void
DOC_Method::warning( File fich, std::string mess ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::warning" ) ;
   int line = ( fich==Implementation ? imp_line : def_line ) ;
   string file = ( fich==Implementation ? imp_fich : def_fich ) ;
   
   if( line>0 )
   {
      XWarningE( file, line, mess ) ;
   }
   else
   {
     Erreur( mess ) ;
   }
}



//--------------------------------------------------------------------
void
DOC_Method::verify( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::verify" ) ;

   // Verify les declarations virtuelles
   if( surdefinition!=0 )
   {
      if( !is_virtual() )
      {
         warning( Definition, "Method " + full_name() +
                  " should be declared virtual but isn't" ) ;
      }
      if( !surdefinition->is_virtual() )
      {
         warning( Definition, "Method " + surdefinition->full_name() +
                  " should be declared virtual but isn't" ) ;
      }
      if( surdefinition->category()->name() != category()->name() )
      {
         warning( Definition, "Method " + full_name() ) ;
         surdefinition->warning( Definition, " redeclares " + surdefinition->full_name() +
                                 " with not the same category name" ) ;
      }
      
   }

   // Verify la syntaxe query-command
   if( def_class()->is_component() )
   {
      if( has_to_document(def_class()) )
      {
         check_header( comment() ) ;
      }
      if( !is_static() && 
          !is_destructor() && !is_constructor() &&
          has_to_document(def_class()) )
      {
	 //   check_header( comment() ) ;
         // c'est un createur
         if( name().substr(0,6)=="create" )
         {
            if( type()->full_type_name() == "void" )
            {
               warning( Definition,
                        "a method whose name starts with create should return what is created " ) ;
            }
            check_command_header( comment() ) ;
         }
         // c'est a_ query
         else if( type()->full_type_name() != "void" )
         {
            if( !is_constant &&
                ( category()->name() != "Assignment attempt" ) )
            {
               warning( Definition,
                        "Command "
                        + full_name()
                        + " should not return something ( except assignment attempt )" ) ;
            }
            check_query_header( comment() ) ;
         }
         // c'est a_ command
         else
         {
            check_command_header( comment() ) ;
         }
      }
   }

   // Coherence de la definition avec l'implementation
   if( implementation!= 0 )
   {
      PEL_Iterator *it = args->create_iterator(0) ;
      PEL_Iterator *itOther = implementation->args->create_iterator(0) ;
      itOther->start() ;
      for( it->start() ; it->is_valid() ; it->go_next() )
      {
         DOC_Argument*arg = static_cast<DOC_Argument*>( it->item() ) ;
         DOC_Argument*argOther = static_cast<DOC_Argument*>( itOther->item() ) ;
         if( arg->silent_var() != argOther->silent_var() )
         {
            string me = "The parameter list in the definition of the method " ;
            me += implementation->full_name() ;
            warning( Implementation, me ) ;
            me = "is different from that of its declaration :\n " ;
            me += full_name() ;
            warning( Definition, me ) ;
            break ;
         
         }
         itOther->go_next() ;
      }
      it->destroy() ;
      itOther->destroy() ;
   }

   // PRE et POST conditions reservées aux méthodes publiques
   if( !has_to_document(def_class()) )
   {
      if( category()->name()!="Hidden" )
      {
//          if( preConditions->count()>0 && !is_precondition() )
//          {
//             warning( Implementation, 
//                      "The preconditions of the non public method " + full_name() +
//                      " should not be implemented with the macro \"PEL_CHECK_PRE\", but "
//                      " with the macro \"PEL_CHECK\"" ) ;
//          }
      }
   }
   else
   {
      if( label().empty() && !is_abstract()
          && !is_precondition() && !is_postcondition() && !is_invariant() )
      {
//         warning( Implementation, "Implementation for method " + full_name() +
//                  " should declare label" ) ;
      }
      else if( !label().empty() && !is_abstract()
               && ( is_precondition() || is_postcondition() || is_invariant() ) )
      {
         warning( Implementation, "PEL_LABEL not allowed in " + full_name() ) ;
         
      }
   }
   
   if( !label().empty() )
   {
      if( !( label().find( def_class()->name() )<label().length() &&
             label().find( name() )<label().length() ) )
      {
         warning( Implementation, "Bad label for " + full_name() + 
                  " : " + label() ) ;
      }
   }
   // Problem in polymorphism of pre and post conditions
   if( is_precondition() || is_postcondition() )
   {
      if( referring_item==0 ) 
      {
         DOC_Method const* a_method = surdefinition ;
         
         while( a_method!=0 ) 
         {
            if( a_method->referring_item!=0 )
            {
               referring_item = a_method->referring_item ;
            }
            a_method = a_method->surdefinition ;
         }
         if( referring_item==0 )
         {
            std::string tx = ( is_precondition() ? "precondition" : "postcondition" ) ;
            warning( Implementation, 
                  "The "+tx+" method " + full_name() +
                     " is not referred by any method" ) ;
         }
      }
   }
}



//--------------------------------------------------------------------
bool
DOC_Method::is_equal( PEL_Object const* other ) const 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::is_equal" ) ;
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;
   DOC_Method const* m = static_cast<DOC_Method const*>( other ) ;

   bool result = name()== m->name();
   result = result && args->count() == m->args->count() &&
      is_constant==m->is_constant ;
   if( result )
   {
      PEL_Iterator* it = args->create_iterator(0) ;
      it->start() ;
      PEL_Iterator* itOther = m->args->create_iterator(0) ;
      itOther->start() ;
      for( ;it->is_valid() ; it->go_next() )
      {
         DOC_Argument* a = static_cast<DOC_Argument*>( it->item() ) ;
         DOC_Argument* oa = static_cast<DOC_Argument*>( itOther->item() ) ;
         
         result = result && a->type()->is_equal( oa->type() ) ;
         itOther->go_next() ;
      }
      it->destroy() ;
      itOther->destroy() ;
      
   }

   return result ;
}


//--------------------------------------------------------------------
void
DOC_Method:: check_header( std::string const& header ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::check_header" ) ;
//?????????? il y a beaucoup des redites avec formate_comment???????????
//?????????? et ca n'est pas forcement coherent
//?????????? et ca doit pouvoir etre utilise pour process_documentationr des classes
   size_t n = header.length() ;
   size_t start, stop ;

   start = header.find( '`' ) ;
   while( start<n )
   {
      stop = header.find( '\'', start+1 ) ;
      if( stop>=n) warning( Definition, " unbalanced quotes " ) ;
      string identif = header.substr( start+1, stop-start-1 ) ;
      size_t debpp = identif.find("::") ;
      size_t iend = identif.find( '(' ) ;
      if( iend >= identif.length() ) iend = identif.length() ;
      if( debpp < identif.length() )
      {
         DOC_Class const* cl = def_class() ;
	 if( debpp!=0 )
	 {
            string cc = identif.substr( 0, debpp ) ;
            cl = DOC_Class::search( cc ) ;
            if( cl  == 0 )
	    {
               warning( Definition, " unknown class : " + cc ) ;
	    }
	 }
         if( debpp != identif.length()-2 )
	 {
            string mm = identif.substr( debpp+2, iend-debpp-2 ) ;
            if( cl!=0 && cl->find( mm ) == 0 )
	    {
               warning( Definition, " unknown method : " + mm ) ;
	    }
         }
      }
      else if( identif != "self"   )
      {
         if( !has_argument( identif ) )
               warning( Definition, " unknown method argument : " + identif ) ;
      }
      start = header.find( '`', stop+1 ) ;
   }
}



//--------------------------------------------------------------------
void
DOC_Method:: check_query_header( std::string const& header ) const
//--------------------------------------------------------------------
{
//   if( count( header.begin(), header.end(), '.') != 0 )
//   {
//      warning( Definition, 
//               " header comments for queries should not end with a period " ) ;
//   }
   PEL_LABEL( "DOC_Method::check_query_header" ) ;
   istringstream is( header.c_str() ) ;
   string first_word ;
   if( is >> first_word )
   {
      string last_word = first_word ;
      string word ;
      while( is>>word && word!="IMPLEMENTATION" )
      {
         last_word = word ;
      }
      char first_char = first_word[0] ;
      char last_char  = last_word[ last_word.length()-1 ] ;
      if( last_char=='?' )
      {
         if( !isupper(first_char) )
            warning( Definition,
                     " invalid boolean query header comment " ) ;
      }
      else
      {
         if( !( !isupper(first_char) && last_char!='.' ) ) 
            warning( Definition,
                     " invalid non-boolean query header comment " ) ;
      }
   }
}



//--------------------------------------------------------------------
void
DOC_Method:: check_command_header( std::string const& header ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Method::check_command_header" ) ;

   istringstream is( header.c_str() ) ;
   string first_word ;
   if( is >> first_word )
   {
      string last_word = first_word ;
      string word ;
      while( is>>word && word!="IMPLEMENTATION" )
      {
         last_word = word ;
      }
      char first_char = first_word[0] ;
      char last_char  = last_word[ last_word.length()-1 ] ;
      if( ! ( last_char=='.' && isupper(first_char) ) )
      {
         warning( Definition,
                     " invalid command header comment " ) ;
      }
   }
}
