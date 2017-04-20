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

#include <DOC_Typedef.hh>

#include <DOC_Argument.hh>
#include <DOC_Category.hh>
#include <DOC_Text.hh>
#include <DOC_Class.hh>
#include <DOC_Function.hh>
#include <DOC_Sequence.hh>
#include <DOC_Tools.hh>
#include <DOC_Type.hh>


//--------------------------------------------------------------------
DOC_Typedef*
DOC_Typedef::create( std::string const& a_name,
                 DOC_Type * a_type,
                 Protection maprotection,
                 DOC_Category const* macategory,
                 DOC_Text * mycomment ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Typedef::create" ) ;
   return new DOC_Typedef( a_name,
                       a_type,
                       maprotection,
                       macategory,
                       mycomment ) ;
   
}

//--------------------------------------------------------------------
DOC_Typedef::DOC_Typedef( std::string const& a_name,
                  DOC_Type * a_type,
                  Protection a_protection,
                  DOC_Category const* a_category,
                  DOC_Text * a_comment ) 
//--------------------------------------------------------------------
   : DOC_ClassItem( a_protection, a_category, a_comment ),
     myName( a_name ),
     myDOC_Type( a_type )
{
}



//--------------------------------------------------------------------
DOC_Typedef::~DOC_Typedef( void ) 
//--------------------------------------------------------------------
{
}



//--------------------------------------------------------------------
string const&
DOC_Typedef::name( void ) const
//--------------------------------------------------------------------
{
   return myName ;
}



//--------------------------------------------------------------------
string 
DOC_Typedef::signature( void ) const
//--------------------------------------------------------------------
{
   string ret = "typedef " + myDOC_Type->full_type_name() + " " + myName ;
   return ret ;
}



//--------------------------------------------------------------------
string 
DOC_Typedef::prototype( DOC_Writer& sullizer ) const
//--------------------------------------------------------------------
{
   string ret = "typedef " + myDOC_Type->type_reference( sullizer ) + " " + myName ;
   return ret ;
}



//--------------------------------------------------------------------
void
DOC_Typedef::display( std::ostream& out ) const 
//--------------------------------------------------------------------
{
   out << signature() ;
}



