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

#include <DOC_Argument.hh>
#include <DOC_Class.hh>
#include <DOC_Type.hh>

//--------------------------------------------------------------------
DOC_Argument*
DOC_Argument:: to_argument( void ) 
//--------------------------------------------------------------------
{
   return this ;
}



//--------------------------------------------------------------------
DOC_Argument*
DOC_Argument:: create( PEL_Object* a_owner,
                       DOC_Type * typ,
                       std::string const& name ) 
//--------------------------------------------------------------------
{
   return new DOC_Argument( a_owner, typ, name ) ;
}


//--------------------------------------------------------------------
DOC_Argument:: DOC_Argument( PEL_Object* a_owner,
                             DOC_Type * typ,
                             std::string const& name ) 
//--------------------------------------------------------------------
   : DOC_Symbol(a_owner),
     nameMuet( name ),
     myDOC_Type( typ )
{
}



//--------------------------------------------------------------------
DOC_Argument:: ~DOC_Argument( void ) 
//--------------------------------------------------------------------
{
}



//--------------------------------------------------------------------
void 
DOC_Argument:: initialize( std::string const& init ) 
//--------------------------------------------------------------------
{
   initialisation = init ;
}



//--------------------------------------------------------------------
DOC_Type const*
DOC_Argument:: type( void ) const
//--------------------------------------------------------------------
{
   return myDOC_Type ;
}



//--------------------------------------------------------------------
std::string const& 
DOC_Argument:: silent_var( void ) const
//--------------------------------------------------------------------
{
   return nameMuet ;
}



//--------------------------------------------------------------------
std::string
DOC_Argument:: text( void ) const
//--------------------------------------------------------------------
{
   string ret = nameMuet ;
   if( !initialisation.empty() )
   {
      ret += " = " ;
      ret += initialisation ;
   }
   return ret ;
}



//--------------------------------------------------------------------
bool 
DOC_Argument:: is_equal( PEL_Object const* other ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Argument:: is_equal" ) ;
   PEL_CHECK_PRE( is_equal_PRE( other ) ) ;
   DOC_Argument const* arg = static_cast<DOC_Argument const*>( other ) ;
      
   bool result = myDOC_Type->is_equal( arg->myDOC_Type ) &&
      nameMuet==arg->nameMuet ;
   
   return result ;
}
