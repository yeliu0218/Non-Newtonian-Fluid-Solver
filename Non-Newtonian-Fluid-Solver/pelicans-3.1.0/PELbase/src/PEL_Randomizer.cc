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

#include <PEL_Randomizer.hh>

#include <PEL_assertions.hh>

#include <size_t_vector.hh>
#include <iostream>

//----------------------------------------------------------------------
static const int IM1=2147483563 ;
static const int IM2=2147483399 ;
static const double AM=1./IM1 ;
static const int IMM1=IM1-1 ;
static const int IA1=40014 ;
static const int IA2=40692 ;
static const int IQ1=53668 ;
static const int IQ2=52774 ;
static const int IR1=12211 ;
static const int IR2=3791 ;
static const int NTAB=32 ;
static const int NDIV=1+IMM1/NTAB ;
static const double EPS=1.2E-7 ;
static const double RNMX=1.-EPS ;
//----------------------------------------------------------------------


//----------------------------------------------------------------------
PEL_Randomizer*
PEL_Randomizer:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Randomizer:: create_clone" ) ;
   PEL_Randomizer* result =  new PEL_Randomizer( a_owner, my_series ) ;

   PEL_CHECK_POST( create_clone_POST( result,  a_owner ) ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
PEL_Randomizer*
PEL_Randomizer:: create( PEL_Object* a_owner, int series )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Randomizer:: create" ) ;
   PEL_Randomizer* result =  new PEL_Randomizer( a_owner, series ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
PEL_Randomizer:: PEL_Randomizer( PEL_Object* a_owner,
                               int series )
//----------------------------------------------------------------------
   : PEL_Object( a_owner ),
     my_series( series ),
     IV(NTAB)
{
   PEL_LABEL( "PEL_Randomizer:: PEL_Randomizer" ) ;
   start() ;
}



//----------------------------------------------------------------------
PEL_Randomizer:: ~PEL_Randomizer( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Randomizer:: ~PEL_Randomizer" ) ;
}



//----------------------------------------------------------------------
void
PEL_Randomizer:: start( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Randomizer:: start" ) ;
   IDUM = my_series ;
   IDUM2 = my_series ;
   IY = 0 ;
   
   for( size_t j=NTAB+8; j>=1; j-- )
   {
      int K=IDUM/IQ1 ;
      IDUM=IA1*(IDUM-K*IQ1)-K*IR1 ;
      if( IDUM<0 ) IDUM=IDUM+IM1 ;
      if( j<=NTAB ) IV(j-1)=IDUM ;
   }
   IY=IV(0) ;
   go_next() ;
}



//----------------------------------------------------------------------
void
PEL_Randomizer:: go_next( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Randomizer:: go_next" ) ;
   
   int K=IDUM/IQ1 ;
   IDUM=IA1*(IDUM-K*IQ1)-K*IR1 ;
   if (IDUM<0) IDUM=IDUM+IM1 ;
   K=IDUM2/IQ2 ;
   IDUM2=IA2*(IDUM2-K*IQ2)-K*IR2 ;
   if(IDUM2<0) IDUM2=IDUM2+IM2 ;
   int J=1+IY/NDIV ;
   IY=IV(J-1)-IDUM2 ;
   IV(J-1)=IDUM ;
   if (IY<1) IY = IY+IMM1 ;
   value=AM*IY ;
   if (value>RNMX) value = RNMX ;
   
}

//----------------------------------------------------------------------
double
PEL_Randomizer:: item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Randomizer:: item" ) ;

   double result = value ;
   
   PEL_CHECK_POST( 0.0 <= result ) ;
   PEL_CHECK_POST( result <= 1.0 ) ;
   
   return( result ) ;
}


//----------------------------------------------------------------------
void
PEL_Randomizer:: build_permutation( size_t_vector& vec ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_Randomizer:: build_permutation" ) ;
   size_t size = vec.size() ;
   for( size_t j=0 ; j<size ; j++ )
   {
      vec(j)=j ;
   }
   
   for( size_t j=0 ; j<size-1 ; j++ )
   {
      size_t nj = size-j ;
      double r = item()*nj ;
      go_next() ;
      int ij = (int) r ;
      if( ij>0 )
      {
         size_t ja = j + ij ;
         PEL_CHECK( ja<size ) ;
         size_t tmp = vec(j) ;
         vec(j) = vec(ja) ;
         vec(ja) = tmp ;
      }
   }
   PEL_CHECK_POST( FORALL( (size_t i=0;i<vec.size();i++),
                           vec(i)<vec.size() ) ) ;
   PEL_CHECK_POST( FORALL( (size_t i=0;i<vec.size();i++),
                           FORALL( (size_t j=0;j<vec.size();j++),
                                   IMPLIES(i!=j,vec(i)!=vec(j) ) ) ) ) ;
   
}


