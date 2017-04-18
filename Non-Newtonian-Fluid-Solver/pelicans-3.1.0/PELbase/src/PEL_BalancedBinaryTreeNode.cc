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

#include <PEL_BalancedBinaryTreeNode.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>

//-------------------------------------------------------------------------
PEL_BalancedBinaryTreeNode:: PEL_BalancedBinaryTreeNode( PEL_Object * value )
//-------------------------------------------------------------------------
   :   PEL_Object(0), height(0), val(value), less(0), more(0)
{
   PEL_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
PEL_BalancedBinaryTreeNode:: ~PEL_BalancedBinaryTreeNode( void )
//-------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
   clear() ;
}



//-------------------------------------------------------------------------
PEL_Object *
PEL_BalancedBinaryTreeNode:: search( const PEL_Object * value ) const
//-------------------------------------------------------------------------
{
   PEL_CHECK_PRE( value!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_Object * ret = 0 ;
   int cmp = val->three_way_comparison(value) ;
   
   if( cmp < 0 )
   {
      if( more!=0 )
      {
         ret = more->search( value );
      }
   }
   else
   {
      if( cmp == 0 )
      {
         ret = val ;
      }
      else
      {
         if( less!=0 )
         {
            ret = less->search( value ) ;
         }
      }
   }
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( ret==0 || ret->is_equal( value ) ) ;
   return( ret ) ;
}



//-------------------------------------------------------------------------
PEL_Object *
PEL_BalancedBinaryTreeNode:: addChild(
   PEL_Object * value,
   PEL_BalancedBinaryTreeNode * & top ) 
//-------------------------------------------------------------------------
{
   PEL_CHECK_PRE( value!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_Object * ret = 0 ;
   if( val == 0 )
   {
      val = value ;
      ret = val ;
   }
   else
   {
      int cmp = val->three_way_comparison( value ) ;
   
      if( cmp<0 )
      {
         if( more!=0 )
         {
            ret = more->addChild( value, more ) ;
         }
         else
         {
            more = new PEL_BalancedBinaryTreeNode( value ) ;
            ret = value ;
         }
      }
      else
      {
         if( cmp==0 )
         {
            ret=val ;
         }
         else
         {
            if( less!=0 )
            {
               ret = less->addChild( value, less )  ;
            }
            else
            {
               less =  new PEL_BalancedBinaryTreeNode( value ) ;
               ret = value ;
            }
         }
      }
      balance( top ) ;
   }
   PEL_CHECK_INV( invariant() ) ;
   return( ret ) ;
}



//-------------------------------------------------------------------------
PEL_BalancedBinaryTreeNode *
PEL_BalancedBinaryTreeNode:: remove(
   const PEL_Object * value,
   PEL_BalancedBinaryTreeNode * & top ) 
//-------------------------------------------------------------------------
{
   PEL_CHECK_PRE( value!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;
   PEL_BalancedBinaryTreeNode * ret = 0 ;
   int cmp = val->three_way_comparison( value ) ;
   
   if( cmp<0 )
   {
      if( more!=0 )
      {
         ret = more->remove( value, more ) ;
         balance( top ) ;
      }
      else
      {
         PEL_Error::object()->raise_plain( "Remove : error : item doesn't exist" ) ;
      }
   }
   else if( cmp>0 )
   {
      if( less!=0 )
      {
         ret = less->remove( value, less )  ;
         balance( top ) ;
      }
      else
      {
         PEL_Error::object()->raise_plain( "Remove : error : item doesn't exist" ) ;
      }
   }
   else // cmp == 0
   {
      top = 0 ;
      if( less!=0 )
      {
         top = less ;
         if( more!=0 )
         {
            less->copySubTree( more, top ) ;
         }
      }
      else if( more!=0 )
      {
         top = more ;
      }
      ret = this ;
      less = 0 ;
      more = 0 ;
   }
   PEL_CHECK_INV( invariant() ) ;   
   return( ret ) ;
}



//-------------------------------------------------------------------------
void
PEL_BalancedBinaryTreeNode:: balance( PEL_BalancedBinaryTreeNode * & top )
//-------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
   // Balancing
   int lessH = -1 ;
   if( less!=0 )
   {
      lessH = less->height ;
   }
   int moreH = -1 ;
   if( more!=0 )
   {
      moreH = more->height ;
   }
   //calculHeight() ;
   if( lessH > moreH+1 )
   {
      //PEL::out() << "Before less balancing : " << *top ;
      PEL_BalancedBinaryTreeNode * lessMore = less->more ;
      less->more = this ;
      if( this==top )
      {
         top = less ;
      }
      PEL_BalancedBinaryTreeNode * oldLess = less ;
      less = lessMore ;
      calculHeight() ;
      oldLess->calculHeight() ;
      //PEL::out() << "After balancing : " << *top ;
   } 
   else if( moreH > lessH+1 )
   {
      //PEL::out() << "Before more balancing : " << *top ;
      PEL_BalancedBinaryTreeNode * moreLess = more->less ;
      more->less = this ;
      if( this==top )
      {
         top = more ;
      }
      PEL_BalancedBinaryTreeNode * oldMore = more ;
      more = moreLess ;
      calculHeight() ;
      oldMore->calculHeight() ;
      //PEL::out() << "After balancing : " << *top ;
   }
   else
   {
      calculHeight() ;
   }
   PEL_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
void
PEL_BalancedBinaryTreeNode:: copySubTree( PEL_BalancedBinaryTreeNode * source,
                                          PEL_BalancedBinaryTreeNode *& top)
//-------------------------------------------------------------------------
{
   if( source->less )
   {
      copySubTree( source->less, top ) ;
   }
   addChild( source->val, top ) ;
   if( source->more )
   {
      copySubTree( source->more, top ) ;
   }
}

      
//-------------------------------------------------------------------------
void
PEL_BalancedBinaryTreeNode:: calculHeight( void )
//-------------------------------------------------------------------------
{
   int lessH = -1 ;
   if( less!=0 )
   {
      lessH = less->height ;
   }
   int moreH = -1 ;
   if( more!=0 )
   {
      moreH = more->height ;
   }
   height=(size_t)PEL::max( lessH, moreH )+1 ;
}



//-------------------------------------------------------------------------
void
PEL_BalancedBinaryTreeNode:: clear( void )
//-------------------------------------------------------------------------
{
   val = 0 ;
   
   if( more!=0 )
   {
      more->destroy() ;
      more = 0 ;
   }
   if( less!=0 )
   {
      less->destroy() ;
      less = 0 ;
   }
}



//-------------------------------------------------------------------------
size_t
PEL_BalancedBinaryTreeNode:: getNbElem( void ) const
//-------------------------------------------------------------------------
{
   // These methodis used for invariant verification purpose
   size_t nbElem = 1 ;
   PEL_ASSERT( val!=0 ) ;
   
   if( less!=0 )
   {
      nbElem += less->getNbElem() ;
   }
   if( more!=0 )
   {
      nbElem += more->getNbElem() ;
   }
   return nbElem ;
}



//-------------------------------------------------------------------------
bool
PEL_BalancedBinaryTreeNode:: invariant( void ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   PEL_ASSERT( val!=0 ) ;

   return( true ) ;
}
