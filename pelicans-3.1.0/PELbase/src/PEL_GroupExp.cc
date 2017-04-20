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

#include <PEL_GroupExp.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>

#include <doubleVector.hh>

#include <iostream>

bool PEL_GroupExp::OPT_EVAL = false ;

PEL_GroupExp const* 
PEL_GroupExp:: PROTOTYPE_UNIT_SORT = new PEL_GroupExp( "unit_sort",
                                                       unit_sort ) ;
PEL_GroupExp const* 
PEL_GroupExp:: PROTOTYPE_SEGM_SORT = new PEL_GroupExp( "segm_sort",
                                                       segm_sort ) ;
PEL_GroupExp const* 
PEL_GroupExp:: PROTOTYPE_SEGM2D_SORT = new PEL_GroupExp( "segm2D_sort",
                                                         segm2D_sort ) ;
PEL_GroupExp const* 
PEL_GroupExp:: PROTOTYPE_SEGM3D_SORT = new PEL_GroupExp( "segm3D_sort",
                                                         segm3D_sort ) ;

struct PEL_GroupExp_ERROR
{
   static void n0( std::string const& f_name, std::string const& op_name ) ;
} ;

//----------------------------------------------------------------------
void
PEL_GroupExp:: set_optimized_evaluation( void )
//----------------------------------------------------------------------
{
   OPT_EVAL = true ;
}

//----------------------------------------------------------------------
void
PEL_GroupExp:: unset_optimized_evaluation( void )
//----------------------------------------------------------------------
{
   OPT_EVAL = false ;
}

//----------------------------------------------------------------------
PEL_GroupExp:: PEL_GroupExp( std::string const& a_name, GroupOp a_op ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , OP( a_op )
   , INITIALIZED( false )
   , X( 0 )
   , Y( 0 )
   , Z( 0 )
{
}

//----------------------------------------------------------------------
PEL_GroupExp:: PEL_GroupExp( PEL_Object* a_owner,
			     std::string const& a_name,
			     PEL_Sequence const* argument_list,
                             GroupOp a_op ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , OP( a_op )
   , INITIALIZED( false )
   , X( 0 )
   , Y( 0 )
   , Z( 0 )
{
   PEL_LABEL( "PEL_GroupExp:: PEL_GroupExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_GroupExp:: ~PEL_GroupExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_GroupExp:: ~PEL_GroupExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
   switch( OP )
   {
      case unit_sort :
         PROTOTYPE_UNIT_SORT = 0 ;
         break ;
      case segm_sort :
         PROTOTYPE_SEGM_SORT = 0 ;
         break ;
      case segm2D_sort :
         PROTOTYPE_SEGM2D_SORT = 0 ;
         break ;
      case segm3D_sort :
         PROTOTYPE_SEGM3D_SORT = 0 ;
         break ;
      default :
         PEL_GroupExp_ERROR::n0( "~PEL_GroupExp", name() ) ;
         break ;
   }
}

//----------------------------------------------------------------------
PEL_GroupExp*
PEL_GroupExp:: create_replica( PEL_Object* a_owner,
                               PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_GroupExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_GroupExp* result = new PEL_GroupExp( a_owner, 
                                            name(), 
                                            argument_list,
                                            OP ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
PEL_GroupExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result = "undefined" ;
   switch( OP )
   {
      case unit_sort :
         result = name() + "(DS,DS,DS,IS[,BS])" ;
         break ;
      case segm_sort :
         result = name() + "(DS,DV,IS[,BS]) " ;
         break ;
      case segm2D_sort :
         result = name() + "(<DV with 2 elems>,DV,IS,DV,IS[,BS])" ;
         break ;
      case segm3D_sort :
         result = name() + "(<DV with 3 elems>,DV,IS,DV,IS,DV,IS[,BS])" ;
         break ;
      default :
         PEL_GroupExp_ERROR::n0( "usage", name() ) ;
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_GroupExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_GroupExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   bool result = false ;
   switch( OP )
   {
      case unit_sort :
         result = ( some_arguments->count() == 4 ||
                    some_arguments->count() == 5 ) ;
         if( result )
         {
            result &=
               extract_arg(some_arguments,0)->data_type() == PEL_Data::Double ;
            result &=
               extract_arg(some_arguments,1)->data_type() == PEL_Data::Double ;
            result &=
               extract_arg(some_arguments,2)->data_type() == PEL_Data::Double ;
            result &=
               extract_arg(some_arguments,3)->data_type() == PEL_Data::Int ;
         }
         if( some_arguments->count() == 5 )
         {
            result &=
               extract_arg(some_arguments,4)->data_type() == PEL_Data::Bool ;
         }   
         break ;
      case segm_sort :
         result = ( some_arguments->count() == 3 ||
                    some_arguments->count() == 4 ) ;
         if( result )
         {
            result &=
               extract_arg(some_arguments,0)->data_type() == PEL_Data::Double ;
            result &=
               extract_arg(some_arguments,1)->data_type() == PEL_Data::DoubleVector ;
            result &=
               extract_arg(some_arguments,2)->data_type() == PEL_Data::Int ;
         }
         if( some_arguments->count() == 4 )
         {
            result &=
               extract_arg(some_arguments,3)->data_type() == PEL_Data::Bool ;
         } 
         break ;
      case segm2D_sort :
         result = ( some_arguments->count() == 5 ||
                    some_arguments->count() == 6 ) ;
         if( result )
         {
            result &=
               extract_arg(some_arguments,0)->data_type() == PEL_Data::DoubleVector ;
            result &=
               extract_arg(some_arguments,1)->data_type() == PEL_Data::DoubleVector ;
            result &=
               extract_arg(some_arguments,2)->data_type() == PEL_Data::Int ;
            result &=
               extract_arg(some_arguments,3)->data_type() == PEL_Data::DoubleVector ;
            result &=
               extract_arg(some_arguments,4)->data_type() == PEL_Data::Int ;
         }
         if( some_arguments->count() == 6 )
         {
            result &=
               extract_arg(some_arguments,5)->data_type() == PEL_Data::Bool ;
         } 
         break ;
      case segm3D_sort :
         result = ( some_arguments->count() == 7 ||
                    some_arguments->count() == 8 ) ;
         if( result )
         {
            result &=
               extract_arg(some_arguments,0)->data_type() == PEL_Data::DoubleVector ;
            result &=
               extract_arg(some_arguments,1)->data_type() == PEL_Data::DoubleVector ;
            result &=
               extract_arg(some_arguments,2)->data_type() == PEL_Data::Int ;
            result &=
               extract_arg(some_arguments,3)->data_type() == PEL_Data::DoubleVector ;
            result &=
               extract_arg(some_arguments,4)->data_type() == PEL_Data::Int ;
            result &=
               extract_arg(some_arguments,5)->data_type() == PEL_Data::DoubleVector ;
            result &=
               extract_arg(some_arguments,6)->data_type() == PEL_Data::Int ;
         }
         if( some_arguments->count() == 8 )
         {
            result &=
               extract_arg(some_arguments,7)->data_type() == PEL_Data::Bool ;
         }
         break ;
      default :
         PEL_GroupExp_ERROR::n0( "matches_args", name() ) ;
         break ;
   } 
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_GroupExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   return Int ;
}

//----------------------------------------------------------------------
int
PEL_GroupExp:: to_int( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_GroupExp:: to_int" ) ;
   PEL_CHECK_PRE( to_int_PRE( ct ) ) ;

   int result = PEL::bad_int() ;

   switch( OP )
   {
      case unit_sort :
         {
            double vx = arg(0)->to_double( ct ) ;
            double const v1 = arg(1)->to_double( ct ) ;
            double const v2 = arg(2)->to_double( ct ) ;
            int const n = arg(3)->to_int( ct ) ;
            if( v1>=v2 )
               raise_error( "arg(2)<arg(3) expected" ) ;
            if( vx<v1 || vx>v2 )
               raise_error( "arg(2)<=arg(1)<=arg(3) expected" ) ;
            if( n<=0 )
               raise_error( "arg(4)>0 expected" ) ;
            bool shift = false ;
            if( nb_arguments() == 5 ) shift = arg(4)->to_bool( ct ) ;
            if( shift )
            {
               vx += 0.499999/(v2-v1)/n ;
               if( vx>v2 ) vx -= v2-v1 ;
            }
            PEL_ASSERT( vx>=v1 && vx<=v2 ) ;
            result = (int)(n*(vx-v1)/(v2-v1)) ;
            if( result>=n ) result=n-1 ;
            if( result<0 ) result=0 ;
         }
         break ;
      case segm_sort :
         {
            if( !INITIALIZED )
            {
               initialize( "arg(2)", arg(1)->to_double_vector( ct ),
                           "arg(3)", arg(2)->to_int( ct ),
                           X ) ;
               if( OPT_EVAL ) INITIALIZED = true ;
            }
            bool shift = false ;
            if( nb_arguments() == 4 ) shift = arg(3)->to_bool( ct ) ;
            result = index( arg(0)->to_double( ct ), X, shift ) ;
         }
         break ;
      case segm2D_sort :
         {
            if( !INITIALIZED )
            {
               initialize( "arg(2)", arg(1)->to_double_vector( ct ),
                           "arg(3)", arg(2)->to_int( ct ),
                           X ) ;
               initialize( "arg(4)", arg(3)->to_double_vector( ct ),
                           "arg(5)", arg(4)->to_int( ct ),
                           Y ) ;
               if( OPT_EVAL ) INITIALIZED = true ;
            }
            doubleVector const& vv = arg(0)->to_double_vector( ct ) ;
            if( vv.size() != 2 )
               raise_error( "arg(1) of size 2 expected" ) ;
            bool shift = false ;
            if( nb_arguments() == 6 ) shift = arg(5)->to_bool( ct ) ;
            int const nx = arg(2)->to_int( ct ) ;
            int const ix = index( vv(0), X, shift ) ;
            // int const ny = arg(4)->to_int( ct ) ;
            int const iy = index( vv(1), Y, shift ) ;
            result = iy*nx+ix ;
         }
         break ;
      case segm3D_sort :
         {
            if( !INITIALIZED )
            {
               initialize( "arg(2)", arg(1)->to_double_vector( ct ),
                           "arg(3)", arg(2)->to_int( ct ),
                           X ) ;
               initialize( "arg(4)", arg(3)->to_double_vector( ct ),
                           "arg(5)", arg(4)->to_int( ct ),
                           Y ) ;
               initialize( "arg(6)", arg(5)->to_double_vector( ct ),
                           "arg(7)", arg(6)->to_int( ct ),
                           Z ) ;
               if( OPT_EVAL ) INITIALIZED = true ;
            }
            doubleVector const& vv = arg(0)->to_double_vector( ct ) ;
            if( vv.size() != 3 )
               raise_error( "arg(1) of size 3 expected" ) ;
            bool shift = false ;
            if( nb_arguments() == 8 ) shift = arg(7)->to_bool( ct ) ;
            int const nx = arg(2)->to_int( ct ) ;
            int const ix = index( vv(0), X, shift ) ;
            int const ny = arg(4)->to_int( ct ) ;
            int const iy = index( vv(1), Y, shift ) ;
            // int const nz = arg(6)->to_int( ct ) ;
            int const iz = index( vv(2), Z, shift ) ;
            result = iz*nx*ny+iy*nx+ix ;
         }
         break ;
      default :
         PEL_GroupExp_ERROR::n0( "to_int", name() ) ;
         break ;
   }
   PEL_CHECK_INV( invariant() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_GroupExp:: initialize( std::string const& v_arg, doubleVector const& v,
                           std::string const& n_arg, int n,
                           doubleVector& x_table ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_GroupExp:: initialize" ) ;
   PEL_CHECK( !INITIALIZED ) ;
   PEL_CHECK( !v_arg.empty() ) ;
   PEL_CHECK( !n_arg.empty() ) ;

   if( n<=0 )
      raise_error( n_arg+">0 expected" ) ;
   size_t const vsize = v.size() ;
   if( vsize<=1 )
      raise_error( v_arg+" should have at least 2 elements" ) ;
   for( size_t i=0 ; i<vsize-1 ; ++i )
   {
      if( v(i)>=v(i+1) )
         raise_error( v_arg+": increasing table of values expected" ) ;
   }
   
   x_table.re_initialize( n+1 ) ;
   x_table(0) = v(0) ;
   double const dx = (v.size()-1.)/( (double) n ) ;
   for( size_t i=1 ; i<(size_t) n ; ++i )
   {
      double const idx = i*dx ;
      size_t const j = (size_t) idx ;
      double const alpha = idx-j ;
      x_table(i) = (1.-alpha)*v(j)+alpha*v(j+1) ;
   }
   x_table(n) = v( v.size()-1 ) ;
}

   
//----------------------------------------------------------------------
int
PEL_GroupExp:: index( double x, doubleVector const& x_table,
                      bool shift ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_GroupExp:: index" ) ;

   size_t const n = x_table.size()-1 ;
   if( x<x_table(0) || x>x_table(n) )
      raise_error( "coordinates not in range" ) ;
   
   if( shift )
   {
      x += 0.499999*( x_table(1)-x_table(0) ) ;
      if( x>x_table(n) )
      {
         x -= x_table(n)-x_table(0) ;
      }
   }
   size_t i=0 ;
   size_t j=n ;
   while( j-i>1 )
   {
      size_t k=(i+j)/2 ;
      double xx= x_table(k) ;
      if( x>=xx )
      {
         i = k ;
      }
      else
      {
         j = k ;
      }
   }
   int const result = i ;
   
   PEL_CHECK_POST( result>=0 && result<(int) x_table.size()-1 ) ;
   PEL_CHECK_POST( x>=x_table(result) && x<=x_table(result+1) ) ;
   return( result ) ;
}

//internal--------------------------------------------------------------
void 
PEL_GroupExp_ERROR:: n0( std::string const& f_name,
                         std::string const& op_name )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** PEL_GroupExp::"+f_name+"\n" ;
   mesg += "    operation "+op_name+" not implemented." ;
   PEL_Error::object()->raise_internal( mesg ) ;
}
