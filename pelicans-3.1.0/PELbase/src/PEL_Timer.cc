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

#include <PEL_Timer.hh>

#include <PEL.hh>
#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_System.hh>
#include <PEL_assertions.hh>

#include <iostream>

//-------------------------------------------------------------------------
PEL_Timer* PEL_Timer:: create( PEL_Object* a_owner )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Timer:: create" ) ;

   PEL_Timer* result = new PEL_Timer( a_owner ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
PEL_Timer:: PEL_Timer( PEL_Object* a_owner )
//-------------------------------------------------------------------------
   : PEL_Object( a_owner )
   , CUMUL_TIME( 0. )
   , CURRENT_TIME( PEL::bad_double() )
   , RUNNING( false )
   , CUMUL_ELAPSED_TIME( 0. )
   , CURRENT_ELAPSED_TIME( PEL::bad_double() )
{
   PEL_CHECK_INV( invariant() ) ;
}

//-------------------------------------------------------------------------
PEL_Timer:: ~PEL_Timer( void )
//-------------------------------------------------------------------------
{
   PEL_CHECK_INV( invariant() ) ;
   if( is_running() )
   {
      stop() ;
   }
}

//-------------------------------------------------------------------------
void
PEL_Timer:: start( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Timer:: start" ) ;
   PEL_CHECK_PRE( !is_running() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   CURRENT_TIME = PEL_System::user_time() ;
   CURRENT_ELAPSED_TIME = PEL_System::epoch_time() ;
   RUNNING = true ;

   PEL_CHECK_INV( invariant() ) ;
}

//-------------------------------------------------------------------------
void
PEL_Timer:: stop( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Timer:: stop" ) ;
   PEL_CHECK_PRE( is_running() ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   CUMUL_TIME += ( PEL_System::user_time() - CURRENT_TIME ) ;
   CURRENT_TIME = PEL::bad_double() ;
   CUMUL_ELAPSED_TIME += ( PEL_System::epoch_time() - CURRENT_ELAPSED_TIME ) ;
   CURRENT_ELAPSED_TIME = PEL::bad_double() ;
   RUNNING = false ;
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( !is_running() ) ;
}

//-------------------------------------------------------------------------
void
PEL_Timer:: reset( void )
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Timer:: reset" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   if( is_running() )
   {
      stop() ;
   }
   CUMUL_TIME = 0. ;
   CUMUL_ELAPSED_TIME = 0. ;
   
   PEL_CHECK_INV( invariant() ) ;
}

//-------------------------------------------------------------------------
double
PEL_Timer:: time( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Timer:: time" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   double result = CUMUL_TIME ;
   if( is_running() )
   {
      result += ( PEL_System::user_time() - CURRENT_TIME ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result>=0.0 ) ;
   return( result );
}
   
//-------------------------------------------------------------------------
double
PEL_Timer:: elapsed_time( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Timer:: elapsed_time" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   double result = CUMUL_ELAPSED_TIME ;
   if( is_running() )
   {
      result += ( PEL_System::epoch_time() - CURRENT_ELAPSED_TIME ) ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK_POST( result>=0.0 ) ;
   return( result );
}
   
//-------------------------------------------------------------------------
bool
PEL_Timer:: is_running( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Timer:: is_running" ) ;
   PEL_CHECK_INV( invariant() ) ;
   return( RUNNING ) ;
}

//-------------------------------------------------------------------------
void
PEL_Timer:: print( std::ostream& os, size_t indent_width ) const
//-------------------------------------------------------------------------
{
   static size_t size = PEL_Exec::communicator()->nb_ranks() ;
   
   print_time( time(), os, indent_width ) ;
   if( size>1 ) 
   {
      os << " (" ;
      print_time( elapsed_time(), os, 0 ) ;
      os << ")" ;
   }
}

//-------------------------------------------------------------------------
void
PEL_Timer:: print_time( double a_time, std::ostream& os, size_t indent_width ) 
//-------------------------------------------------------------------------
{
   std::string space( indent_width, ' ' ) ;
   os << space ;
   double t = a_time ;
   
   int const hh = (int) ( t/3600. ) ;
   int const mm = (int) ( (t-3600.*hh)/60. ) ;
   if( hh==0 && mm==0 )
   {
      os << t << " s" ;
   }
   else
   {
      if( hh > 0 )
      {
         os << hh << ":" ;
      }
      if( mm==0 )
      {
         os << "00:" ;
      }
      else if( mm<10 )
      {
         os << "0" << mm << ":" ;
      }
      else
      {
         os << mm << ":" ;
      }
      double const ss = ( (int) (10.*(t-3600.*hh-60.*mm)) )/10. ;
      if( ss==0 )
      {
         os << "00" ;
      }
      else if( ss<10 )
      {
         os << "0" << ss ;
      }
      else
      {
         os << ss ;
      }
      os << " h:m:s" ;
   }
}

//-------------------------------------------------------------------------
bool
PEL_Timer:: invariant( void ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( PEL_Object::invariant() ) ;
   PEL_ASSERT( !( RUNNING && CURRENT_TIME==PEL::bad_double() ) ) ;
   PEL_ASSERT( !( RUNNING && CURRENT_ELAPSED_TIME==PEL::bad_double() ) ) ;
   return( true ) ;
}
