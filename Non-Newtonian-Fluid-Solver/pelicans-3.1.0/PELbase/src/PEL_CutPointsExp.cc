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

#include <PEL_CutPointsExp.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>
#include <doubleVector.hh>
#include <doubleArray2D.hh>

PEL_CutPointsExp const*
PEL_CutPointsExp::PROTOTYPE_MP = new PEL_CutPointsExp(
                         "middle_point", PEL_CutPointsExp::mid_point ) ;
PEL_CutPointsExp const*
PEL_CutPointsExp::PROTOTYPE_MPS = new PEL_CutPointsExp(
                       "middle_points", PEL_CutPointsExp::mid_points ) ;
PEL_CutPointsExp const*
PEL_CutPointsExp::PROTOTYPE_X = new PEL_CutPointsExp(
                              "x_cut_points", PEL_CutPointsExp::xcut ) ;
PEL_CutPointsExp const*
PEL_CutPointsExp::PROTOTYPE_Y = new PEL_CutPointsExp(
                              "y_cut_points", PEL_CutPointsExp::ycut ) ;
PEL_CutPointsExp const*
PEL_CutPointsExp::PROTOTYPE_Z = new PEL_CutPointsExp(
                              "z_cut_points", PEL_CutPointsExp::zcut ) ;

struct PEL_CutPointsExp_ERROR
{
   static void n0( std::string const& f_name, std::string const& op_name ) ;
} ;

//----------------------------------------------------------------------
PEL_CutPointsExp:: PEL_CutPointsExp( std::string const& a_name,
                                     IS_CutOp a_op ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , OP( a_op )
{
   PEL_LABEL( "PEL_CutPointsExp:: PEL_CutPointsExp" ) ;
}

//----------------------------------------------------------------------
PEL_CutPointsExp:: PEL_CutPointsExp( PEL_Object* a_owner,
                                     std::string const& a_name,
                                     PEL_Sequence const* argument_list,
                                     IS_CutOp a_op ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , OP( a_op )
{
   PEL_LABEL( "PEL_CutPointsExp:: PEL_CutPointsExp" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_CutPointsExp:: ~PEL_CutPointsExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CutPointsExp:: ~PEL_CutPointsExp" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( this == PROTOTYPE_X ) PROTOTYPE_X = 0 ;
   if( this == PROTOTYPE_Y ) PROTOTYPE_Y = 0 ;
   if( this == PROTOTYPE_Z ) PROTOTYPE_Z = 0 ;
   if( this == PROTOTYPE_MP ) PROTOTYPE_MP = 0 ;
   if( this == PROTOTYPE_MPS ) PROTOTYPE_MPS = 0 ;
   
}

//----------------------------------------------------------------------
PEL_CutPointsExp*
PEL_CutPointsExp:: create_replica( PEL_Object* a_owner,
                                   PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CutPointsExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_CutPointsExp* result =
           new PEL_CutPointsExp( a_owner, name(), argument_list, OP ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
PEL_CutPointsExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result = "unspecified" ;
   if( OP == mid_point )
   {
      result = name()+"(DS,DV)" ;
   }
   else if( OP == mid_points )
   {
      result = name()+"([DV],DV)" ;
   }
   else if( OP == xcut )
   {
      result = name()+"(DV,DS[,DS])" ;
   }
   else if( OP == ycut )
   {
      result = name()+"(DS,DV[,DS] )" ;
   }
   else if( OP == zcut )
   {
      result = name()+"(DS,DS,DV)" ;
   }
   else
   {
      PEL_CutPointsExp_ERROR::n0( "usage", name() ) ;
   }
   
   return( result ) ;
}

//----------------------------------------------------------------------
bool
PEL_CutPointsExp:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CutPointsExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = false ;
   if( OP == mid_point )
   {
      result = ( some_arguments->count() == 2 ) ;
      if( result )
      {
         PEL_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
         result &= ( k0 == Double ) ;
         PEL_Data::Type k1 = extract_arg( some_arguments, 1 )->data_type() ;
         result &= ( k1 == DoubleVector ) ;
      }
   }
   else if( OP == mid_points )
   {
      result = ( some_arguments->count() == 1 || some_arguments->count() == 2 ) ;
      if( result )
      {
         PEL_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
         result &= ( k0 == DoubleVector ) ;
         if( some_arguments->count() == 2 )
         {
            PEL_Data::Type k1 = extract_arg( some_arguments, 1 )->data_type() ;
            result &= ( k1 == DoubleVector ) ;
         }
      }
   }
   else if( OP == xcut )
   {
      result = ( some_arguments->count() == 2 || some_arguments->count() == 3 ) ;
      if( result )
      {
         PEL_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
         result &= ( k0 == DoubleVector ) ;
         PEL_Data::Type k1 = extract_arg( some_arguments, 1 )->data_type() ;
         result &= ( k1 == Double ) ;
         if( some_arguments->count() == 3 )
         {
            PEL_Data::Type k2 = extract_arg( some_arguments, 2 )->data_type() ;
            result &= ( k2 == Double ) ;
         }
      }
   }
   else if( OP == ycut )
   {
      result = ( some_arguments->count() == 2 || some_arguments->count() == 3 ) ;
      if( result )
      {
         PEL_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
         result &= ( k0 == Double ) ;
         PEL_Data::Type k1 = extract_arg( some_arguments, 1 )->data_type() ;
         result &= ( k1 == DoubleVector ) ;
         if( some_arguments->count() == 3 )
         {
            PEL_Data::Type k2 = extract_arg( some_arguments, 2 )->data_type() ;
            result &= ( k2 == Double ) ;
         }
      }
   }
   else if( OP == zcut )
   {
      result = ( some_arguments->count() == 3 ) ;
      if( result )
      {
         PEL_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
         result &= ( k0 == Double ) ;
         PEL_Data::Type k1 = extract_arg( some_arguments, 1 )->data_type() ;
         result &= ( k1 == Double ) ;
         PEL_Data::Type k2 = extract_arg( some_arguments, 2 )->data_type() ;
         result &= ( k2 == DoubleVector ) ;
      }
   }
   else
   {
      PEL_CutPointsExp_ERROR::n0( "valid_arguments", name() ) ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_CutPointsExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_Data::Type result = Undefined ;
   if( OP == mid_point )
   {
      result = Double ;
   }
   else if( OP == mid_points )
   {
      result = DoubleVector ;
   }
   else if( OP == xcut ||  OP == ycut || OP == zcut )
   {
      result = DoubleArray2D ;
   }
   else
   {
      PEL_CutPointsExp_ERROR::n0( "data_type", name() ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
PEL_CutPointsExp:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CutPointsExp:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE( ct ) ) ;

   double result = PEL::bad_double() ;

   if( OP == mid_point )
   {
      double const x = arg(0)->to_double( ct ) ;
      doubleVector const& verts_table = arg(1)->to_double_vector( ct ) ;
      result = m_pt( x, verts_table ) ;
   }
   else
   {
      PEL_CutPointsExp_ERROR::n0( "to_double", name() ) ;
   }
   
   return( result ) ;
}

//----------------------------------------------------------------------
doubleVector const&
PEL_CutPointsExp:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CutPointsExp:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ;

   static doubleVector result(0) ;

   if( OP == mid_points )
   {
      if( nb_arguments()==1 )
      {
         doubleVector const& verts_table = arg(0)->to_double_vector( ct ) ;
         build_coords( verts_table, result ) ;
      }
      else if( nb_arguments()==2 )
      {
         doubleVector const& pt_table = arg(0)->to_double_vector( ct ) ;
         doubleVector const& verts_table = arg(1)->to_double_vector( ct ) ;
         result.re_initialize( pt_table.size() ) ;
         for( size_t i=0 ; i<pt_table.size() ; ++i )
         {
            result(i) = m_pt( pt_table(i), verts_table ) ;
         }
      }
      else
      {
         PEL_CutPointsExp_ERROR::n0( "to_double_vector", name() ) ;
      }
   }
   else
   {
      PEL_CutPointsExp_ERROR::n0( "to_double_vector", name() ) ;
   }
   
   return( result ) ;
}

//----------------------------------------------------------------------
doubleArray2D const&
PEL_CutPointsExp:: to_double_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CutPointsExp:: to_double_array2D" ) ;
   PEL_CHECK_PRE( to_double_array2D_PRE( ct ) ) ;

   static doubleArray2D result(0,0) ;
   
   size_t const nb_dims = nb_arguments() ;
   if( OP == xcut )
   {
      doubleVector const& verts_table = arg(0)->to_double_vector( ct ) ;
      doubleVector coords_table(0) ;
      build_coords( verts_table, coords_table ) ;
      result.re_initialize( coords_table.size(), nb_dims ) ;
      for( size_t i=0 ; i<coords_table.size() ; ++i )
      {
         result(i,0) = coords_table(i) ;
         result(i,1) = arg(1)->to_double( ct ) ;
         if( nb_dims==3 ) result(i,2) = arg(2)->to_double( ct ) ;
      }
   }
   else if( OP == ycut )
   {
      doubleVector const& verts_table = arg(1)->to_double_vector( ct ) ;
      doubleVector coords_table(0) ;
      build_coords( verts_table, coords_table ) ;
      result.re_initialize( coords_table.size(), nb_dims ) ;
      for( size_t i=0 ; i<coords_table.size() ; ++i )
      {
         result(i,0) = arg(0)->to_double( ct ) ;
         result(i,1) = coords_table(i) ;
         if( nb_dims==3 ) result(i,2) = arg(2)->to_double( ct ) ;
      }
      
   }
   else if( OP == zcut )
   {
      doubleVector const& verts_table = arg(2)->to_double_vector( ct ) ;
      doubleVector coords_table(0) ;
      build_coords( verts_table, coords_table ) ;
      result.re_initialize( coords_table.size(), nb_dims ) ;
      for( size_t i=0 ; i<coords_table.size() ; ++i )
      {
         result(i,0) = arg(0)->to_double( ct ) ;
         result(i,1) = arg(1)->to_double( ct ) ;
         result(i,2) = coords_table(i) ;
      }
   }
   else
   {
      PEL_CutPointsExp_ERROR::n0( "to_double_array2D", name() ) ;
   }
   
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_CutPointsExp:: check_table( doubleVector const& verts_table ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CutPointsExp:: check_table" ) ;

   if( verts_table.size() <= 1 )
   {
      raise_error( "cutline vector should have at least two elements" ) ;
   }

   for( size_t i=0 ; i<verts_table.size()-1 ; ++i )
   {
      if( verts_table(i)>=verts_table(i+1) )
      {
         raise_error( "cutline vector should be increasing" ) ;
      }
   }
}

//----------------------------------------------------------------------
double
PEL_CutPointsExp:: m_pt( double x,
                         doubleVector const& verts_table ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CutPointsExp:: m_pt" ) ;

   check_table( verts_table ) ;

   double result = PEL::bad_double() ;
   for( size_t i=0 ; i<verts_table.size()-1 ; ++i )
   {
      PEL_CHECK( verts_table(i)<verts_table(i+1) ) ;
      result = 0.5*( verts_table(i)+verts_table(i+1) ) ;
      if( x<verts_table(i+1) ) break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_CutPointsExp:: build_coords( doubleVector const& verts_table,
                                 doubleVector& coords_table ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_CutPointsExp:: build_coords" ) ;

   check_table( verts_table ) ;
   
   coords_table.re_initialize( verts_table.size()-1 ) ;
   for( size_t i=0 ; i<verts_table.size()-1 ; ++i )
   {
      coords_table(i) = 0.5*( verts_table(i)+verts_table(i+1) ) ;
   }
}

//internal--------------------------------------------------------------
void 
PEL_CutPointsExp_ERROR:: n0( std::string const& f_name,
                             std::string const& op_name )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** PEL_CutPointsExp::" + f_name +"\n" ;
   mesg += "    operation " + op_name + " not implemented." ;
   PEL_Error::object()->raise_internal( mesg ) ;
}
