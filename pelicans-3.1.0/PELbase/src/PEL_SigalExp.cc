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

#include <PEL_SigalExp.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>

#include <doubleVector.hh>

#include <iostream>

PEL_SigalExp const* 
PEL_SigalExp::PROTOTYPE_Regu =
                          new PEL_SigalExp( Regu, "regular_vector" ) ;

PEL_SigalExp const* 
PEL_SigalExp::PROTOTYPE_Ratio =
                       new PEL_SigalExp( Ratio, "stretched_vector" ) ;

PEL_SigalExp const* 
PEL_SigalExp::PROTOTYPE_Concat = new PEL_SigalExp( Concat, "<<" ) ;

PEL_SigalExp const* 
PEL_SigalExp::PROTOTYPE_GeomSeq = new PEL_SigalExp( Geom, 
                                                    "geometric_sequence" ) ;

//----------------------------------------------------------------------
PEL_SigalExp:: PEL_SigalExp( SigalOperator a_op, std::string const& a_name ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , OP( a_op )
   , RESULT_I( 0 )
   , RESULT_D( 0 )
   , RESULT_S( 0 )
{
   PEL_LABEL( "PEL_SigalExp:: PEL_SigalExp" ) ;
}

//----------------------------------------------------------------------
PEL_SigalExp*
PEL_SigalExp:: create_replica( PEL_Object* a_owner,
                               PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SigalExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_SigalExp* result = new PEL_SigalExp( a_owner, OP, name(), argument_list ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_SigalExp:: PEL_SigalExp( PEL_Object* a_owner,
                             SigalOperator a_op,
                             std::string const& a_name,
                             PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , OP( a_op )
   , RESULT_I( 0 )
   , RESULT_D( 0 )
   , RESULT_S( 0 )
{
   PEL_LABEL( "PEL_SigalExp:: PEL_SigalExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}


//----------------------------------------------------------------------
PEL_SigalExp:: ~PEL_SigalExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SigalExp:: ~PEL_SigalExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
std::string const&
PEL_SigalExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result = "undefined" ;
   switch( OP )
   {
      case Regu :
         result = name() + "(DS,IS,DS) or " + name() + "(IS,IS,IS)" ;
         break ;
      case Ratio :
         result = name()+"(DS,DS,DS,DS)" ;
         break ;
      case Geom :
         result = name()+"(DS,DS,IS)" ;
         break ;
      case Concat :
         result = "(DV|IV|SV) " + name() + " (same type as first arg)" ;
         break ;
      default :
         result = "undefined" ;
         break ;
   }
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_SigalExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SigalExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result ;
   switch( OP )
   {
      case Regu :
         result = some_arguments->count()==3 &&
            ( extract_arg( some_arguments, 0 )->data_type()==Double ||
              extract_arg( some_arguments, 0 )->data_type()==Int ) &&
            extract_arg( some_arguments, 1 )->data_type()==Int &&
            extract_arg( some_arguments, 2 )->data_type()==
                      extract_arg( some_arguments, 0 )->data_type() ;
         break ;
      case Ratio :
         result = some_arguments->count()==4 &&
            extract_arg( some_arguments, 0 )->data_type()==Double &&
            extract_arg( some_arguments, 1 )->data_type()==Double &&
            extract_arg( some_arguments, 2 )->data_type()==Double &&
            extract_arg( some_arguments, 3 )->data_type()==Double ;
         break ;
      case Geom :
         result = some_arguments->count()==3 &&
            extract_arg( some_arguments, 0 )->data_type()==Double &&
            extract_arg( some_arguments, 1 )->data_type()==Double &&
            extract_arg( some_arguments, 2 )->data_type()==Int ;
         break ;
      case Concat :
         result = some_arguments->count()==2 &&
            (
               ( extract_arg( some_arguments, 0 )->data_type()==DoubleVector ||
                 extract_arg( some_arguments, 0 )->data_type()==IntVector ||
                 extract_arg( some_arguments, 0 )->data_type()==StringVector )
               &&
               ( extract_arg( some_arguments, 1 )->data_type()==
                              extract_arg( some_arguments, 0 )->data_type() ) );
         break ;
      default : result = false ;}
         
   return result ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_SigalExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SigalExp:: data_type" ) ;

   PEL_Data::Type result = Undefined ;
   switch( OP )
   {
      case Regu :
         {
            PEL_Data::Type tt = arg( 0 )->data_type() ;
            if( tt == PEL_Data::Double )
            {
               result = PEL_Data::DoubleVector ;
            }
            else if( tt == PEL_Data::Int )
            {
               result = PEL_Data::IntVector ;
            }
         }
         break ;
      case Ratio :
         result = PEL_Data::DoubleVector ;
         break ;
      case Geom :
         result = PEL_Data::DoubleVector ;
         break ;
      case Concat :
         result = arg( 0 )->data_type() ;
         break ;
      default :
         PEL_Error::object()->raise_internal( "not implemented" ) ;
         break ;
   }

   return result ;
}

//----------------------------------------------------------------------
stringVector const&
PEL_SigalExp:: to_string_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SigalExp:: to_string_vector" ) ;
   PEL_CHECK_PRE( to_string_vector_PRE( ct ) ) ;

   stringVector& result = RESULT_S ;

   if( OP == Concat )
   {
      stringVector const& v0 = arg(0)->to_string_vector(ct) ;
      stringVector const& v1 = arg(1)->to_string_vector(ct) ;
      size_t n = v0.size() + v1.size() ;
      RESULT_S.re_initialize( n ) ;
      for( size_t i=0 ; i<v0.size() ; i++ )
      {
         RESULT_S(i) = v0(i) ;
      }
      for( size_t i=0 ; i<v1.size() ; i++ )
      {
         RESULT_S(n-i-1) = v1(v1.size()-i-1) ;
      }
   }
   return result ;
}

//----------------------------------------------------------------------
intVector const&
PEL_SigalExp:: to_int_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SigalExp:: to_int_vector" ) ;
   PEL_CHECK_PRE( to_int_vector_PRE( ct ) ) ;

   intVector& result = RESULT_I ;

   if( OP == Concat )
   {
      intVector const& v0 = arg(0)->to_int_vector(ct) ;
      intVector const& v1 = arg(1)->to_int_vector(ct) ;
      size_t n = v0.size() + v1.size() ;
      if( v0.size()>0 && v1.size() > 0 && v0(v0.size()-1)==v1(0) )
      {
         n-- ;
      }
      RESULT_I.re_initialize( n ) ;
      for( size_t i=0 ; i<v0.size() ; i++ )
      {
         RESULT_I(i) = v0(i) ;
      }
      for( size_t i=0 ; i<v1.size() ; i++ )
      {
         RESULT_I(n-i-1) = v1(v1.size()-i-1) ;
      }
   }
   else
   {
      int alpha = arg( 0 )->to_int( ct ) ;
      int beta  = arg( 2 )->to_int( ct ) ;

      if( alpha==beta )
      {
         raise_error( "initial and final value should be different" ) ;
      }

      int nn = arg( 1 )->to_int( ct ) ;
      if( nn<0 )
      {
         raise_error( "Number of steps should be a positive value" ) ;
      }
      result.re_initialize( nn + 1 ) ;
      int dx = ( beta - alpha )/nn ;
      for( size_t i=0 ; i<(size_t)nn+1 ; ++i )
      {
         result( i ) = alpha + i*dx ;
      }
      if( result( nn ) != beta )
      {
         raise_error( "Bad definition" ) ;
      }
   }

   return( result ) ;
}

//----------------------------------------------------------------------
doubleVector const&
PEL_SigalExp:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SigalExp:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ;
   realize( ct ) ;
   doubleVector& result = RESULT_D ;
   return result ;
}

//----------------------------------------------------------------------
void
PEL_SigalExp:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string const space = std::string( indent_width, ' ' ) ;
   
   switch( OP )
   {
      case Regu :
      case Ratio :
         PEL_Expression::print( os, indent_width ) ;
         break ;
      case Geom :
         PEL_Expression::print( os, indent_width ) ;
         break ;
      case Concat :
         arg(0)->print( os, indent_width ) ;
         os << " " << name() << std::endl ;
         arg(1)->print( os, indent_width ) ;
         break ;
   }   
}

//----------------------------------------------------------------------
void
PEL_SigalExp:: realize( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SigalExp:: realize" ) ;
   if( OP==Regu )
   {
      double alpha = arg(0)->to_double(ct) ;
      double beta =  arg(2)->to_double(ct) ;
      if( alpha==beta )
      {
         raise_error( "initial and final value should be different" ) ;
      }
      doubleVector& v = RESULT_D ;
      int in = arg( 1 )->to_int( ct ) ;
      if( in<0 )
      {
         raise_error( "the mumber of steps should be a positive value" ) ;
      }
      size_t n = (size_t) in ;
      v.re_initialize(n+1) ;
      double dx = (beta-alpha)/n ;
      for( size_t i=1 ; i<n ; i++ )
      {
         v(i) = alpha + i*dx ;
      }
      v(0) = alpha ;
      v(v.size()-1) = beta ; 
   }
   else if( OP==Ratio )
   {
      double const x0 = arg(0)->to_double(ct) ;
      double const dx0 = arg( 1 )->to_double(ct) ;
      double const dx1 = arg( 2 )->to_double(ct) ;
      double const x1 = arg( 3 )->to_double(ct) ;

      // Check:
      if( dx0<=0. )
         raise_error( "length of first interval: should be a positive value" ) ;
      if( dx1<=0. )
         raise_error( "length of last interval: should be a positive value" ) ;
      if( x0 == x1 )
         raise_error( "initial and final values should be different" ) ;

      if( dx0<dx1 )
      {
         build_stretched_vector( x0, dx0, x1, dx1, RESULT_D ) ;
      }
      else
      {
         doubleVector dummy(0) ;
         build_stretched_vector( -x1, dx1, -x0, dx0, dummy ) ;
         size_t const n = dummy.size() ;
         RESULT_D.re_initialize( n ) ;
         for( size_t i=0 ; i<n ; ++i )
         {
            RESULT_D(i) = -dummy( n-i-1 ) ;
         }
      }
   }
   else if( OP==Geom )
   {
      double const x0 = arg( 0 )->to_double( ct ) ;
      double const rr = arg( 1 )->to_double( ct ) ;
      if( rr == 0.0 )
      {
         raise_error( "common ratio : should be a non negative value" ) ;
      }
      int in = arg( 2 )->to_int( ct ) ;
      if( in < 0 )
      {
         raise_error( "the number of intervals should be a positive value" ) ;
      }
      size_t n = (size_t) in ;

      RESULT_D.re_initialize( n+1 ) ;
      RESULT_D( 0 ) = x0 ;
      for( size_t i=0 ; i<(size_t)n ; ++i )
      {
         RESULT_D( i+1 ) = RESULT_D( i ) * rr ;
      }
   }
   else
   {
      doubleVector const& v0 = arg(0)->to_double_vector(ct) ;
      doubleVector const& v1 = arg(1)->to_double_vector(ct) ;
      size_t n = v0.size() + v1.size() ;
      if( v0.size()>0 && v1.size() > 0 && v0(v0.size()-1)==v1(0) )
      {
         n-- ;
      }
      RESULT_D.re_initialize( n ) ;
      for( size_t i=0 ; i<v0.size() ; i++ )
      {
         RESULT_D(i) = v0(i) ;
      }
      for( size_t i=0 ; i<v1.size() ; i++ )
      {
         RESULT_D(n-i-1) = v1(v1.size()-i-1) ;
      }
   }
}

//----------------------------------------------------------------------
void
PEL_SigalExp:: build_stretched_vector( double x0, double dx0,
                                       double x1, double dx1,
                                       doubleVector& result ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SigalExp:: build_stretched_vector" ) ;
   PEL_CHECK( OP==Ratio ) ;
   PEL_CHECK( dx0>0. ) ;
   PEL_CHECK( dx1>0. ) ;
   PEL_CHECK( dx0<=dx1 ) ;

   size_t n = PEL::bad_index() ;
   double pas0 = PEL::bad_double() ;
   double r = PEL::bad_double() ;
      
   tageom( dx0, dx1, PEL::abs(x1-x0), pas0, n ,r ) ;
   if( x1<x0 ) pas0 = -pas0 ;
   result.re_initialize( n+1 ) ;
   double pas = pas0 ;
   double deb = x0 ;
   result(0) = deb ;
   for( size_t i=1 ; i<n+1 ; ++i )
   {
      deb = deb+pas ;
      if( i==n ) deb = x1 ;
      result(i) = deb ;
      pas *= r ;
   }

   PEL_CHECK_POST( result(0) == x0 ) ;
   PEL_CHECK_POST( result(result.size()-1) == x1 ) ;
}

//----------------------------------------------------------------------
void
PEL_SigalExp:: tageom(double dx1, double dxn, double zl,
                      double& dx, size_t& n, double& r ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_SigalExp:: tageom" ) ;
   PEL_CHECK( dx1>0. ) ;
   PEL_CHECK( dxn>0. ) ;
   PEL_CHECK( zl>0. ) ;
   
   double const EPS = 1.E-6 ;

   if( dx1>=zl || dxn>=zl)
   {
      n=1 ;
      dx=zl ;
      r=1. ;
   }
   else if( dx1 == dxn )
   {
      // Constant step
      n=size_t( zl/dx1+0.5 ) ;
      dx=zl/n ;
      r=1. ;
   }
   else
   {
      // Number of intervals
      r=(zl-dx1)/(zl-dxn) ;
      n=size_t( 1.+PEL::log(dxn/dx1)/PEL::log((zl-dx1)/(zl-dxn))+0.5 ) ;
      // Final ratio
      double const a1=dx1/zl ;
      double const an=dxn/zl ;
      bool cont = true ;
      while( cont )
      {
         double f = rtaf( a1, an, n, r ) ;
         double df = ( rtaf( a1, an, n, r+EPS )-f )/EPS ;
         double dr = -f/df ;
         r += dr ;
         cont = ( PEL::abs(dr)>EPS ) ;
      }
      double u = 0. ;
      for( size_t i=0 ; i<n ; ++i )
      {
         u = u*r+1. ;
      }
      dx = zl/u ;
   }
}

//----------------------------------------------------------------------
double
PEL_SigalExp:: rtaf( double a1, double an, size_t n, double r ) const
//----------------------------------------------------------------------
{
   double du = 0. ;
   double u = 0. ;
   for( size_t i=0 ; i<n ; ++i )
   {
      du = u+du*r ;
      u = u*r+1. ;
   }
   return( -2.*a1*(1.-a1*u)*du
           +2.*(PEL::pow( r, n-1 )-an*u)*((n-1)*PEL::pow( r,n-2 )-an*du) ) ;
}
   
