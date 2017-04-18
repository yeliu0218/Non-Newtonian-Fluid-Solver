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

#include <PEL_InterpolExp.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <iostream>
#include <fstream>
#include <sstream>

PEL_InterpolExp const* 
PEL_InterpolExp::PROTOTYPE_LIN_1D = new PEL_InterpolExp(
                                                "interpol", lin_inter_1D ) ;

struct PEL_InterpolExp_ERROR
{
   static void n0( std::string const& f_name, std::string const& op_name ) ;
   static void n1( std::string const& f_name,
                   std::string const& filename,
                   size_t line_nb,
                   std::string const& line ) ;
   static void n2( std::string const& f_name, std::string const& filename ) ;
   static void n3( std::string const& f_name,
                   std::string const& filename,
                   size_t line_nb,
                   double xn,
                   double xn1 ) ;
   static void n4( std::string const& f_name,
                   std::string const& filename ) ;
} ;

//----------------------------------------------------------------------
PEL_InterpolExp:: PEL_InterpolExp( std::string const& a_name,
                                 PEL_InterpolOp a_op ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , OP( a_op )
   , FROM_FILE( false )
   , X1_IS_SET( false )
   , X1(0)
   , F_IS_SET( false )
   , FX1(0)
   , CHECK( false )
{
   PEL_LABEL( "PEL_InterpolExp:: PEL_InterpolExp" ) ;
}

//----------------------------------------------------------------------
PEL_Expression*
PEL_InterpolExp:: create_replica( PEL_Object* a_owner,
                                 PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_InterpolExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_InterpolExp* result = new PEL_InterpolExp( a_owner,
                                                name(), 
                                                argument_list,
                                                OP ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_InterpolExp:: PEL_InterpolExp( PEL_Object* a_owner,
                                   std::string const& a_name,
                                   PEL_Sequence const* argument_list,
                                   PEL_InterpolOp a_op ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , OP( a_op )
   , FROM_FILE( false )
   , X1_IS_SET( false )
   , X1(0)
   , F_IS_SET( false )
   , FX1(0)
   , CHECK( false )
{
   PEL_LABEL( "PEL_InterpolExp:: PEL_InterpolExp" ) ;

   switch( OP )
   {
      case lin_inter_1D :
         FROM_FILE = ( arg(0)->data_type() == String ) ;
         break ;
      default:
         PEL_InterpolExp_ERROR::n0( "PEL_InterpolExp::", name() ) ;
         break ;
   }
   
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_InterpolExp:: ~PEL_InterpolExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_InterpolExp:: ~PEL_InterpolExp" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( PROTOTYPE_LIN_1D == this )
   {
      PROTOTYPE_LIN_1D = 0 ;
   }
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_InterpolExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_InterpolExp:: data_type" ) ;

   PEL_Data::Type result = PEL_Data::Undefined ;
   switch( OP )
   {
      case lin_inter_1D :
         result = PEL_Data::Double ;
         break ;
      default:
         PEL_InterpolExp_ERROR::n0( "data_type", name() ) ;
         break ;
   }     
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
PEL_InterpolExp:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_InterpolExp:: usage" ) ;

   static std::string result = "undefined" ;
   switch( OP )
   {
      case lin_inter_1D :
         result = name() + "(DV,DV,DV) or "+ name() + "(IS,DV) " ;
         break ;
      default:
         PEL_InterpolExp_ERROR::n0( "usage", name() ) ;
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_InterpolExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_InterpolExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = false ;
   switch( OP )
   {
      case lin_inter_1D :
         result = ( some_arguments->count() > 0 ) ;
         if( result )
         {
            PEL_Data::Type k0 =
                             extract_arg( some_arguments, 0 )->data_type() ;
            if( k0 == String )
            {
               result = ( some_arguments->count() == 2 ) ;
               if( result )
               {
                  PEL_Data::Type k1 =
                             extract_arg( some_arguments, 1 )->data_type() ;
                  result = ( k1 == Double ) ;
               }
            }
            else if( k0 == DoubleVector )
            {
               result = ( some_arguments->count() == 3 ) ;
               if( result )
               {
                  PEL_Data::Type k1 =
                             extract_arg( some_arguments, 1 )->data_type() ;
                  result = ( k1 == DoubleVector ) ;
                  PEL_Data::Type k2 =
                             extract_arg( some_arguments, 2 )->data_type() ;
                  result &= ( k2 == Double ) ;
               }
            }
            else
            {
               result = false ;
            }
         }
         break ;
      default:
         PEL_InterpolExp_ERROR::n0( "valid_arguments", name() ) ;
         break ;
   }     
   return( result ) ;
}

//----------------------------------------------------------------------
double
PEL_InterpolExp:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_InterpolExp:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE(ct) ) ;

   double result = PEL::bad_double() ;

   switch( OP )
   {
      case lin_inter_1D :
         {
            double x = PEL::bad_double() ;
            if( FROM_FILE )
            {
               if( !X1_IS_SET )
               {
                  read_tables_1D( arg(0)->to_string( ct ), X1, FX1 ) ;
                  X1_IS_SET = true ;
                  F_IS_SET = true ;
                  CHECK = true ;
               }
               x = arg(1)->to_double( ct ) ;
            }
            else
            {
               if( !X1_IS_SET )
               {
                  X1 = arg(0)->to_double_vector(ct) ;
               }
               if( !F_IS_SET )
               {
                  FX1 = arg(1)->to_double_vector(ct) ;
               }
               if( !X1_IS_SET || !F_IS_SET )
               {
                  check_tables_1D( X1, FX1 ) ;
               }
               X1_IS_SET = true ;
               F_IS_SET = true ;
               CHECK = true ;
               x = arg(2)->to_double( ct ) ;
            }
            result = linear_interpol_1D( X1, FX1, x ) ;
         }
         break ;
      default:
         PEL_InterpolExp_ERROR::n0( "to_double", name() ) ;
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_InterpolExp:: read_tables_1D( std::string const& filename,
                                  doubleVector& X_table,
                                  doubleVector& FX_table ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_InterpolExp:: read_tables_1D" ) ;
   
   std::ifstream in( filename.c_str() ) ;
   if( !in ) PEL_InterpolExp_ERROR::n2( name(), filename ) ;

   size_t n = 0 ;
   double x, y ;
   std::string line ;
   size_t nb = 0 ;
   
   while( getline( in, line ) )
   {
      std::istringstream sin( line ) ;
      
      while( sin >> x )
      {
         if( ! ( sin >> y ) )
         {
            PEL_InterpolExp_ERROR::n1( name(), filename, nb+1, line ) ;
         }
         X_table.append( x ) ;
         FX_table.append( y ) ;
         if( n!=0 &&  X_table(n)<=X_table(n-1) )
         {
            PEL_InterpolExp_ERROR::n3(
                  name(), filename, nb+1, X_table(n), X_table(n-1) ) ;
         }
         n++ ;
      }
      if( !sin.eof() )
      {
         PEL_InterpolExp_ERROR::n1( name(), filename, nb+1, line ) ;
      }
      nb++ ;
   }
   if( n == 0 ) PEL_InterpolExp_ERROR::n4( name(), filename ) ;
}

//----------------------------------------------------------------------
double
PEL_InterpolExp:: linear_interpol_1D( doubleVector const& X_table,
                                      doubleVector const& FX_table,
                                      double x ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_InterpolExp:: linear_interpol_1D" ) ;
   PEL_CHECK( X_table.size() == FX_table.size() ) ;

   double result = FX_table(0) ;
   if( x>X_table(0) )
   {
      size_t i=1 ;
      for( ; i<X_table.size() ; i++ )
      {
         if( x<X_table(i) )
         {
            break ;
         }
      }
      if( i<X_table.size() )
      {
         double const dy = FX_table(i)-FX_table(i-1) ;
         double const dx = X_table(i)-X_table(i-1) ;
         result = FX_table(i-1) + ( x-X_table(i-1) )*dy/dx ;
      }
      else
      {
         result = FX_table(i-1) ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_InterpolExp:: check_tables_1D( doubleVector const& X_table,
                                   doubleVector const& FX_table ) const
//----------------------------------------------------------------------
{
   if( X_table.size() != FX_table.size() )
   {
      raise_error( "Same size expected for abscissea and ordinates tables" ) ;
   }
   for( size_t i=1 ; i<X_table.size() ; ++i )
   {
      if( X_table(i-1)>=X_table(i) )
      {
         raise_error( "Increasing first table expected" ) ;
      }
   }
}

//internal--------------------------------------------------------------
void 
PEL_InterpolExp_ERROR:: n0( std::string const& f_name,
                            std::string const& op_name )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** PEL_InterpolExp::"+f_name+"\n" ;
   mesg += "    operation "+op_name+" not implemented." ;
   PEL_Error::object()->raise_internal( mesg ) ;
}

//internal--------------------------------------------------------------
void 
PEL_InterpolExp_ERROR:: n1( std::string const& f_name,
                            std::string const& filename,
                            size_t line_nb,
                            std::string const& line )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** PEL_InterpolExp::"+f_name+" error\n" ;
   mesg << "    Syntax error reading file : \""+filename+"\"\n" ;
   mesg << "    Error at line " ;
   mesg << line_nb << " : \"" ;
   mesg << line << "\"" << std::endl ;
   
   
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PEL_InterpolExp_ERROR:: n2( std::string const& f_name,
                            std::string const& filename )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** PEL_InterpolExp::"+f_name+" error\n" ;
   mesg += "    Unable to open \""+filename+"\"." ;
   PEL_Error::object()->raise_plain( mesg ) ;
}

//internal--------------------------------------------------------------
void 
PEL_InterpolExp_ERROR:: n3( std::string const& f_name,
                            std::string const& filename,
                            size_t line_nb,
                            double xn,
                            double xn1 )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   
   mesg << "*** PEL_InterpolExp::"+f_name+" error\n" ;
   mesg << "    Bad file \""+filename+"\"\n" ;
   mesg << "    At line ";
   mesg << line_nb ;
   mesg << " increasing first table expected but: "  ;
   PEL::print_double( mesg, xn ) ;
   mesg << " < " ;
   PEL::print_double( mesg, xn1 ) ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
PEL_InterpolExp_ERROR:: n4( std::string const& f_name,
                            std::string const& filename )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** PEL_InterpolExp::"+f_name+" error\n" ;
   mesg += "    syntax error reading file \""+filename+"\"" ;
   mesg += " : empty file." ;
    
   PEL_Error::object()->raise_plain( mesg ) ;
}
