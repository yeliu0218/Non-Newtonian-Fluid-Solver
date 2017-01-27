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

#include <GE_PerturbatedMeshingExp.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Randomizer.hh>
#include <PEL_Sequence.hh>

#include <doubleVector.hh>

#include <sstream>

GE_PerturbatedMeshingExp const*  
GE_PerturbatedMeshingExp::PROTOTYPE_PERT_COORDS =
   new GE_PerturbatedMeshingExp( GE_PerturbatedMeshingExp::pert_coords,
                                 "perturbated_coordinates" ) ;

struct GE_PerturbatedMeshingExp_ERROR
{
   static void n0( void ) ;
   static void n1( double coeff ) ;
   static void n2( double h_min ) ;
} ;

//----------------------------------------------------------------------
GE_PerturbatedMeshingExp:: GE_PerturbatedMeshingExp(
                   PerturbatedMeshingExp id, std::string const& a_name ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , ID( id )
   , RAND( 0 )
{
}

//----------------------------------------------------------------------
GE_PerturbatedMeshingExp*
GE_PerturbatedMeshingExp:: create_replica( PEL_Object* a_owner,
                                           PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerturbatedMeshingExp:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   GE_PerturbatedMeshingExp* result =
      new GE_PerturbatedMeshingExp( a_owner, ID, name(), argument_list ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
GE_PerturbatedMeshingExp:: GE_PerturbatedMeshingExp(
                                PEL_Object* a_owner, 
                                PerturbatedMeshingExp id,
                                std::string const& a_name,
                                PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , ID( id )
   , RAND( PEL_Randomizer::create( this, 14021972 ) )
{
   PEL_LABEL( "GE_PerturbatedMeshingExp:: GE_PerturbatedMeshingExp" ) ;
   PEL_CHECK_INV( invariant() ) ;

   RAND->start() ;
}

//----------------------------------------------------------------------
GE_PerturbatedMeshingExp:: ~GE_PerturbatedMeshingExp( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerturbatedMeshingExp:: ~GE_PerturbatedMeshingExp" ) ;
   PEL_CHECK_INV( invariant() ) ;

   if( is_a_prototype() )
   {
      if( ID == pert_coords ) PROTOTYPE_PERT_COORDS = 0 ;
   }
}

//----------------------------------------------------------------------
std::string const& 
GE_PerturbatedMeshingExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result = "undefined" ;
   switch( ID )
   {
      case pert_coords :
         result = name()+"(DS,DV,DS,BS)" ;
	 break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
GE_PerturbatedMeshingExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   return( PEL_Data::DoubleVector ) ;
}

//----------------------------------------------------------------------
doubleVector const&
GE_PerturbatedMeshingExp:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerturbatedMeshingExp:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ;

   static doubleVector result( 0 ) ;

   switch( ID )
   {
      case pert_coords :
         result = arg(1)->to_double_vector( ct ) ;
         size_t nb_sp_dims = result.size() ;
         if( nb_sp_dims != 1 && nb_sp_dims != 2 && nb_sp_dims != 3 )
         {
            GE_PerturbatedMeshingExp_ERROR::n0() ;
         }
         if( ! arg(3)->to_bool( ct ) )
         {
            double const coeff = arg(0)->to_double( ct ) ;
            if( coeff<0. || coeff>=0.5 )
            {
               GE_PerturbatedMeshingExp_ERROR::n1( coeff ) ;
            }
            double const h_min = arg(2)->to_double( ct ) ;
            if( h_min<=0. )
            {
               GE_PerturbatedMeshingExp_ERROR::n2( h_min ) ;
            }
            
            double const phi = 2. * PEL::pi() * random_value() ;
            double const radius = coeff*h_min ;
            if( nb_sp_dims == 1 )
            {
               double const s = ( phi<PEL::pi() ? 1. : -1. ) ;
               result(0) += radius*s ;
            }
            else if( nb_sp_dims == 2 )
            {
               result( 0 ) += radius*PEL::cos( phi ) ;
               result( 1 ) += radius*PEL::sin( phi ) ;
            }
            else if( nb_sp_dims == 3 )
            {
               double const theta = 2. * PEL::pi() * random_value() ;
               result( 0 ) += radius*PEL::sin( theta )*PEL::cos( phi ) ;
               result( 1 ) += radius*PEL::sin( theta )*PEL::sin( phi ) ;
               result( 2 ) += radius*PEL::cos( theta ) ;
            }
         }
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
GE_PerturbatedMeshingExp:: valid_arguments(
                             PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerturbatedMeshingExp:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = false ;
   switch( ID )
   {
      case pert_coords :
         result =
            some_arguments->count() == 4 &&
            extract_arg( some_arguments, 0 )->data_type() == PEL_Data::Double &&
            extract_arg( some_arguments, 1 )->data_type() == PEL_Data::DoubleVector &&
            extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double &&
            extract_arg( some_arguments, 3 )->data_type() == PEL_Data::Bool ;
	 break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
GE_PerturbatedMeshingExp:: random_value( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "GE_PerturbatedMeshingExp:: random_value" ) ;
   
   double result = RAND->item() ;
   RAND->go_next() ;
   
   PEL_CHECK_POST( result >=0.0 && result <= 1.0 ) ;
   return( result ) ;
}

//internal--------------------------------------------------------------
void 
GE_PerturbatedMeshingExp_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** GE_PerturbatedMeshingExp error\n" ;
   mesg += "    the position vector (second argument) should have\n" ;
   mesg += "    1, 2 or 3 dimensions (space number of dimensions)" ;
   PEL_Error::object()->raise_plain( mesg ) ;
}

//internal--------------------------------------------------------------
void 
GE_PerturbatedMeshingExp_ERROR:: n1( double coeff )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** GE_PerturbatedMeshingExp error\n"
        << "    the perturbation coefficient (first argument) should be\n"
        << "    greater than 0. and strickly lower than 0.5\n"
        << "        current value: " << coeff ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
GE_PerturbatedMeshingExp_ERROR:: n2( double h_min )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** GE_PerturbatedMeshingExp error\n"
        << "    the reference distance (fourth argument) should be\n"
        << "    greater than 0.\n"
        << "        current value: " << h_min ;
   PEL_Error::object()->raise_plain( mesg.str() ) ;
}
