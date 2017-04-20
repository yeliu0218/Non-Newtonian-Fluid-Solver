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

#ifndef PEL_ASSERTIONS_HH
#define PEL_ASSERTIONS_HH

#include <PEL_export.hh>

#include <string>

class PEL_EXPORT PEL_Assertion
{
   public: //------------------------------------------------------------

      enum CheckType
      {          None          = 0,
         Precondition  = 1,
         Postcondition = 2,
         Invariant     = 4,
         Check         = 8,
         Objects       = 16
      } ;

      static PEL_Assertion& object( void ) ;

      ~PEL_Assertion( void ) ;

      static void add_handled_check( CheckType a_chec ) ;
      static bool is_handling_check( CheckType some_check ) ;

      static bool is_checking( void ) ;

      static bool test_implement_handling( const char* file,
                                           int line,
                                           std::string const& reason ) ;
      static bool action( const char* file, int line, const char* text ) ;
      static bool do_eval( bool& shortCut ) ;

      friend PEL_Assertion const& operator!( PEL_Assertion const& a ) ;
      friend bool operator&&( bool left, PEL_Assertion const& a ) ;
      friend bool operator||( bool left, PEL_Assertion const& a ) ;
      operator bool( void ) const ;

      static bool& push_bool( void ) ;
      static bool  pop_bool( void ) ;

      static CheckType current_check ;
      static bool result ;
      static bool negation ;
      static bool short_cut ;

   protected: //---------------------------------------------------------

   private: //-----------------------------------------------------------

      PEL_Assertion( void ) ;
      PEL_Assertion( PEL_Assertion const& other ) ;
      PEL_Assertion& operator=( PEL_Assertion const& other ) ;

      static PEL_Assertion unique_instance ;
      static bool eval ;
      static size_t const MAXBOOL = 256 ;
      static bool bool_table[ MAXBOOL ] ;
      static size_t nb_bool ;
      static int checking_level ;
} ;


class PEL_EXPORT PEL_Marker
{
   public: //-----------------------------------------------------------

      PEL_Marker( char const* name ) ;
     ~PEL_Marker( void ) ;

      static bool is_collective( int line ) ;
      
      static size_t nb_labels( void ) ;
      static char const* label( size_t i ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      PEL_Marker( void ) ;
      PEL_Marker( PEL_Marker const& other ) ;
      PEL_Marker const& operator=( PEL_Marker const& other ) ;

      static const char* ring[] ;
      static size_t ring_pos ;
} ;


//----------------------------------------------------------------------

#define asstLoopTest(asstFor,asstAll,predicate)  \
        PEL_Assertion::object() ; \
        { \
           if( PEL_Assertion::do_eval( PEL_Assertion::push_bool() ) ) \
           { \
              PEL_Assertion::push_bool() = PEL_Assertion::negation ; \
              PEL_Assertion::result = asstAll ; \
              for asstFor \
              { \
                 PEL_Assertion::result = predicate ; \
                 if( PEL_Assertion::result != asstAll ) break ; \
              } \
              if( PEL_Assertion::pop_bool() ) \
                 PEL_Assertion::result = !PEL_Assertion::result ;  \
          } \
          PEL_Assertion::short_cut = PEL_Assertion::pop_bool() ; \
        } \
        PEL_Assertion::result = PEL_Assertion::short_cut ? \
        PEL_Assertion::result : PEL_Assertion::result

#define PEL_ASSERTCOND(predicate,text) \
   PEL_Assertion::result = predicate, \
   PEL_Assertion::result || (PEL_Assertion::action(__FILE__,__LINE__,text))

#define PEL_CHECKUNCOND(predicate,check_type,text) \
   if( PEL_Assertion::current_check == PEL_Assertion::None ) \
   { \
      PEL_Assertion::current_check = check_type ; \
      PEL_ASSERTCOND(predicate,text) ; \
      PEL_Assertion::current_check = PEL_Assertion::None ; \
   }

#define PEL_CHECKCOND(predicate,check_type,text) \
   if( PEL_Assertion::is_handling_check( check_type ) && \
       ( PEL_Assertion::current_check == PEL_Assertion::None ) ) \
   { \
      PEL_Assertion::current_check = check_type ; \
      PEL_ASSERTCOND(predicate,text) ; \
      PEL_Assertion::current_check = PEL_Assertion::None ; \
   }

#if !defined(LEVEL) || LEVEL<0 || LEVEL>2
error
Macro LEVEL must be set when compiling Pelicans (opt: -DLEVEL=<value>).
LEVEL value meaning is :
0 : no assertion checking
1 : only preconditions will be checked (recommended).
2 : preconditions, postconditions, invariant and simple checks will be tested
    depending on dynamic level.
#endif

//----------------------------------------------------------------------
// PUBLIC MACROS
//----------------------------------------------------------------------

#define PEL_TEST_IMPLEMENTATION(X) \
   if(X) PEL_Assertion::test_implement_handling( __FILE__, __LINE__, #X )

#define PEL_ASSERT(X) PEL_ASSERTCOND(X,#X)

//--------
#if LEVEL>=2
#define PEL_CHECK_POST(X)  PEL_CHECKCOND(X,PEL_Assertion::Postcondition,#X)
#define OLD(asstname) old_##asstname
#define PEL_SAVEOLD(assttype,asstname,args) assttype old_##asstname = args
#define PEL_CHECK_INV(X) PEL_CHECKCOND(X,PEL_Assertion::Invariant,#X)
#define PEL_CHECK(X) PEL_CHECKCOND(X,PEL_Assertion::Check,#X)
#else
#define PEL_CHECK_POST(X) {}
#define OLD(asstname) {}
#define PEL_SAVEOLD(assttype,asstname,args) {}
#define PEL_CHECK_INV(X) {}
#define PEL_CHECK(X) {}
#endif
//--------

//--------
#if LEVEL>=1
#define PEL_CHECK_PRE(X) PEL_CHECKUNCOND(X,PEL_Assertion::Precondition,#X)
#define PEL_LABEL(X) PEL_Marker aSpy(X)
#define PEL_CHECK_COLLECTIVE(X) if(X){PEL_ASSERTCOND(PEL_Marker::is_collective(__LINE__),"unsynchronized collective operation") ; }
#else
#define PEL_CHECK_PRE(X) {}
#define PEL_LABEL(X) {}
#define PEL_CHECK_COLLECTIVE(X) {}
#endif
//--------

#define FORMAL(X) true
#define FORALL(asstFor,predicate)  asstLoopTest(asstFor,true,predicate)
#define EXISTS(asstFor,predicate)  asstLoopTest(asstFor,false,predicate)
#define IMPLIES(pred,res)          ( !(pred) || (res) )
#define EQUIVALENT(pred,res)       ( (pred) && (res) ) || ( !(res) && !(pred) )

#undef LEVEL

#ifndef OUTLINE
   #include <PEL_assertions.icc>
#endif

#endif
