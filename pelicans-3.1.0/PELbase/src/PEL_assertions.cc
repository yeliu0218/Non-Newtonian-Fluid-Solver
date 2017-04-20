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

#include <PEL_assertions.hh>

#include <PEL_Communicator.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL.hh>

#include <iostream>

#ifdef OUTLINE
#define inline
#include <PEL_assertions.icc>
#undef inline
#endif

//no doc---------------------------------------------------------------------
PEL_Assertion::CheckType PEL_Assertion::current_check = None ;
PEL_Assertion PEL_Assertion::unique_instance ;
bool PEL_Assertion::negation ;
bool PEL_Assertion::result ;
bool PEL_Assertion::eval ;
bool PEL_Assertion::short_cut ;
bool PEL_Assertion::bool_table[ MAXBOOL ] ;
int PEL_Assertion::checking_level = 0 ;
size_t PEL_Assertion::nb_bool = 0 ;
//no doc---------------------------------------------------------------------

//no doc---------------------------------------------------------------------
bool
PEL_Assertion:: test_implement_handling( char const* file,
                                         int line,
                                         std::string const& text )
//no doc---------------------------------------------------------------------
{
   PEL_Error::object()->raise_not_tested( file, line, text ) ;
   return false ;
}

//no doc---------------------------------------------------------------------
bool
PEL_Assertion:: action( const char* file, int line, const char* text )
//no doc---------------------------------------------------------------------
{
   CheckType check = current_check ;
   current_check = None ;
   switch( check )
   {
      case Check :
         PEL_Error::object()->raise_assertion_violation( file, line, text ) ;
         break ;
      case Precondition :
	 PEL_Error::object()->raise_precondition_violation( text ) ;
         break ;
      case Postcondition :
	 PEL_Error::object()->raise_postcondition_violation( file, line, text ) ;
         break ;
      case Invariant :
	 PEL_Error::object()->raise_invariant_violation( file, line, text ) ;
         break ;
      default :
	 PEL_Error::object()->raise_assertion_violation( file, line, text ) ;
         break ;
   }
   return false ;
}

//no doc---------------------------------------------------------------------
const char* PEL_Marker::ring[ 256 ] ;
size_t PEL_Marker::ring_pos = 0 ;
//no doc---------------------------------------------------------------------

//no doc---------------------------------------------------------------------
void
PEL_Assertion:: add_handled_check( CheckType a_check )
//no doc---------------------------------------------------------------------
{
   checking_level |= a_check ;
   PEL_ASSERT( true && ( checking_level & a_check ) ) ;
}

//no doc---------------------------------------------------------------------
bool
PEL_Marker:: is_collective( int line )
//no doc---------------------------------------------------------------------
{   
   static int order = 1 ;
   
   bool result = PEL_Exec::communicator()->same_value_everywhere((int)line*order) ;
   order++ ;

   return result ;
}
