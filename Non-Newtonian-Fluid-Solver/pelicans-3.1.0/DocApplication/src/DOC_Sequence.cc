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

#include <DOC_Sequence.hh>
#include <PEL_assertions.hh>

//--------------------------------------------------------------------
DOC_Sequence*
DOC_Sequence:: to_sequence( void ) 
//--------------------------------------------------------------------
{
   return this ;
}


//--------------------------------------------------------------------
DOC_Sequence*
DOC_Sequence:: create( PEL_Object* a_owner ) 
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Sequence::create" ) ;
   
   return new DOC_Sequence( a_owner ) ;
}



//--------------------------------------------------------------------
DOC_Sequence:: DOC_Sequence( PEL_Object* a_owner ) 
//--------------------------------------------------------------------
      : DOC_Symbol( a_owner ),
        maDOC_Sequence(PEL_List::create(this))
{
}



//--------------------------------------------------------------------
DOC_Sequence:: ~DOC_Sequence( void ) 
//--------------------------------------------------------------------
{
}



//--------------------------------------------------------------------
PEL_List*
DOC_Sequence:: list( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Sequence::list" ) ;
   
   PEL_List* result = maDOC_Sequence ;
   PEL_CHECK_POST( result!=0 ) ;
   
   return maDOC_Sequence ;
}




//--------------------------------------------------------------------
std::string
DOC_Sequence:: text( void ) const
//--------------------------------------------------------------------
{
   PEL_LABEL( "DOC_Sequence::text" ) ;
   
   std::string ret ;
   PEL_Iterator* it=maDOC_Sequence->create_iterator(0) ;
   bool prem = true ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      DOC_Symbol*s = static_cast<DOC_Symbol*>( it->item() ) ;
      if( !prem )
      {
         ret += ", " ;
      }
      else
      {
         prem = false ;
      }
      ret += s->text() ;
   }
   it->destroy() ;
   return ret ;
}
