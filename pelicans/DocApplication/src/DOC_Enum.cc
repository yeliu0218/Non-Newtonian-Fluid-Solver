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

#include <DOC_Enum.hh>

#include <DOC_Argument.hh>
#include <DOC_Category.hh>
#include <DOC_Text.hh>
#include <DOC_Class.hh>
#include <DOC_Function.hh>
#include <DOC_Sequence.hh>
#include <DOC_Tools.hh>
#include <DOC_Type.hh>


//--------------------------------------------------------------------
DOC_Enum*
DOC_Enum::create( std::string const& a_name,
              Protection maprotection,
              DOC_Category const* macategory,
              DOC_Text * mycomment,
              PEL_List const* lst ) 
//--------------------------------------------------------------------
{
   return new DOC_Enum( a_name, maprotection, macategory, mycomment,  lst ) ;
}


//--------------------------------------------------------------------
DOC_Enum::DOC_Enum( std::string const& a_name,
            Protection maprotection,
            DOC_Category const* macategory,
            DOC_Text * mycomment,
            PEL_List const* lst ) 
//--------------------------------------------------------------------
   : DOC_ClassItem( maprotection, macategory, mycomment ),
     myName( a_name ),
     my_list( lst )
{
}



//--------------------------------------------------------------------
DOC_Enum::~DOC_Enum( void ) 
//--------------------------------------------------------------------
{
}



//--------------------------------------------------------------------
string const&
DOC_Enum::name( void ) const
//--------------------------------------------------------------------
{
   return myName ;
}

//--------------------------------------------------------------------
bool
DOC_Enum::declare( std::string const& a__name ) const
//--------------------------------------------------------------------
{
   bool result = DOC_ClassItem::declare( a__name ) ;
   PEL_Iterator* it=my_list->create_iterator(0) ;
   
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      DOC_Symbol* s=static_cast<DOC_Symbol* >( it->item() ) ;
      result = result || s->text()==a__name ;
      
   }
   it->destroy() ;
   return result ;
}


//--------------------------------------------------------------------
string 
DOC_Enum::signature( void ) const
//--------------------------------------------------------------------
{
   string ret = "enum " + myName + " { " ;
   bool prem = true ;
   PEL_Iterator* it=my_list->create_iterator(0) ;
   
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      DOC_Symbol* s=static_cast<DOC_Symbol* >( it->item() ) ;
      if(!prem) ret+= ", " ;
      prem = false ;
      ret += s->text() ;
   }
   it->destroy() ;
   ret += " }" ;
   
   return ret ;
}



//--------------------------------------------------------------------
string 
DOC_Enum::prototype( DOC_Writer& sullizer ) const
//--------------------------------------------------------------------
{
   return signature() ;
}



//--------------------------------------------------------------------
void
DOC_Enum::display( std::ostream& out ) const 
//--------------------------------------------------------------------
{
   out << signature() ;
}






