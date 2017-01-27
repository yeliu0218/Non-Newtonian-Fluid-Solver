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

#include <PEL_MembershipExp.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>
#include <doubleVector.hh>
#include <intVector.hh>

PEL_MembershipExp const*  
PEL_MembershipExp::PROTOTYPE_in_range = 
                             new PEL_MembershipExp( in_range, "in_range" ) ;

PEL_MembershipExp const*  
PEL_MembershipExp::PROTOTYPE_in_box = 
                             new PEL_MembershipExp( in_box, "in_box" ) ;

//----------------------------------------------------------------------
PEL_MembershipExp:: PEL_MembershipExp( MembExp exp_id, 
                                       std::string const& a_name  ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , OP( exp_id )
   , ARG0( 0 )
   , ARG1( 0 )
   , ARG2( 0 )
{
   PEL_LABEL( "PEL_MembershipExp:: PEL_MembershipExp" ) ;
}

//----------------------------------------------------------------------
PEL_MembershipExp*
PEL_MembershipExp:: create_replica( PEL_Object* a_owner,
                                    PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MembershipExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_MembershipExp* result = new PEL_MembershipExp( a_owner, 
                                                      OP,
                                                      name(),
                                                      argument_list ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_MembershipExp:: PEL_MembershipExp( PEL_Object* a_owner,
                                       MembExp exp_id,
                                       std::string const& a_name,
                                       PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , OP( exp_id )
   , ARG0( arg(0) )
   , ARG1( arg(1) )
   , ARG2( exp_id == in_box ? arg(2) : 0 )
{
   PEL_LABEL( "PEL_MembershipExp:: PEL_MembershipExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_MembershipExp:: ~PEL_MembershipExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MembershipExp:: ~PEL_MembershipExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
std::string const& 
PEL_MembershipExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   switch( OP )
   {
      case in_range :
         result =
           "in_range(DS|IS,DV|IV)" ;
	 break ;
      case in_box :
         result = 
           "in_box(DV,DV,DV)" ;
	 break ;
   }
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_MembershipExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MembershipExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = false ;
   switch( OP )
   {
      case in_range :
         result = some_arguments->count()==2 ;
	 if( result )
	 {
	    Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
	    Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
	    result = result && ( ( t0==Double && t1==DoubleVector ) || 
                                 ( t0==Int    && t1==IntVector    ) ) ;
	 }
	 break ;
      case in_box :
         result = ( some_arguments->count() == 3 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
            Type t2 = extract_arg( some_arguments, 2 )->data_type() ;
            result = result && ( t0 == DoubleVector && 
                                 t1 == DoubleVector && 
                                 t2 == DoubleVector ) ;
         }
	 break ;
      default : result = false ;
   }
   return result ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_MembershipExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   return PEL_Data::Bool ;
}

//----------------------------------------------------------------------
bool
PEL_MembershipExp:: to_bool( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MembershipExp:: to_bool" ) ;
   PEL_CHECK_PRE( to_bool_PRE( ct ) ) ;

   bool result = false ;
   switch( OP )
   {
      case in_range :
         double x ;
         double x_min ;
         double x_max ;
         if( ARG0->data_type()==Double )
         {
            x = ARG0->to_double( ct ) ;
            if( ARG1->to_double_vector( ct ).size()!=2 )
               raise_error( "second argument: two components expected" ) ;
            x_min = ARG1->to_double_vector( ct )( 0 ) ;
            x_max = ARG1->to_double_vector( ct )( 1 ) ;
         }
         else
         {
            x = ARG0->to_int( ct ) ;
            if( ARG1->to_int_vector( ct ).size()!=2 )
               raise_error( "second argument: two components expected" ) ;
            x_min = ARG1->to_int_vector( ct )(0) ;
            x_max = ARG1->to_int_vector( ct )(1) ;
         }
         if( x_min >= x_max )
            raise_error( "second argument: increasing values expected" ) ;
         result = ( x_min<=x  && x<=x_max ) ;
         break ;
      case in_box :
         doubleVector const& xx     = ARG0->to_double_vector( ct ) ;
         doubleVector const& xx_min = ARG1->to_double_vector( ct ) ;
         doubleVector const& xx_max = ARG2->to_double_vector( ct ) ;
         if( xx_min.size() != xx.size() )
            raise_error( "first and second arguments: incompatible sizes" ) ;
         if( xx_max.size() != xx.size() )
            raise_error( "first and third arguments: incompatible sizes" ) ;
         result = true ;
         for( size_t i=0 ; i<xx.size() ; ++i )
         {
            if( xx_min(i) >= xx_max(i) )
               raise_error(
                  "components of the third argument should be greater"
                  "\n    than those of the second one" ) ;
            result = result && ( xx_min(i) <= xx(i)     ) ;
            result = result && ( xx(i)     <= xx_max(i) ) ;
         }
         break ;
   }
   return result ;
}
