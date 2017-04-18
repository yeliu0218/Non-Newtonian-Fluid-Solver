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

#include <RS_TransientDiffusion1EXP.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <iostream>

RS_TransientDiffusion1EXP const* 
RS_TransientDiffusion1EXP::PROTOTYPE = new RS_TransientDiffusion1EXP() ;

//----------------------------------------------------------------------
RS_TransientDiffusion1EXP:: RS_TransientDiffusion1EXP( void ) 
//----------------------------------------------------------------------
   : PEL_Expression( "RS_TransientDiffusion1EXP" )
{
}

//----------------------------------------------------------------------
RS_TransientDiffusion1EXP*
RS_TransientDiffusion1EXP:: create_replica(
          PEL_Object* a_owner, PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_TransientDiffusion1EXP:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   RS_TransientDiffusion1EXP* result =
               new RS_TransientDiffusion1EXP( a_owner, argument_list ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
RS_TransientDiffusion1EXP:: RS_TransientDiffusion1EXP(
                 PEL_Object* a_owner, PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, "RS_TransientDiffusion1EXP", argument_list )
{
   PEL_LABEL( "RS_TransientDiffusion1EXP:: RS_TransientDiffusion1EXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
RS_TransientDiffusion1EXP:: ~RS_TransientDiffusion1EXP( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_TransientDiffusion1EXP:: ~RS_TransientDiffusion1EXP" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
RS_TransientDiffusion1EXP:: data_type( void ) const
//----------------------------------------------------------------------
{
   return( Double ) ;
}

//----------------------------------------------------------------------
double
RS_TransientDiffusion1EXP:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_TransientDiffusion1EXP:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE( ct ) ) ; 

   double const epsilon = 1.e-14 ;
   size_t const iter_max = 10000 ;
   
   double const time = arg(0)->to_double(ct) ;

   double result = 0. ;
   if( time > 1.e-8 )
   {
      double const ll = 1.0 ;
      double const a  = 0.25 ;

      double const lsdat = ll/2.0/PEL::sqrt(a*time) ;
   
      double const temp_0   = 0.0 ;
      double const temp_inf = 1.0 ;

      double sgn = -1.0 ;
      bool cv = false ;
      for( size_t i=0 ; !cv && i<iter_max ; ++i )
      {
         sgn = -sgn ;
         double const uu = sgn * PEL::erfc( (2.0*i+1.0)*lsdat ) ;
         result +=  uu ;
         cv = ( PEL::abs(uu) < epsilon ) ;
      }
      if( !cv )
      {
         PEL_Error::object()->raise_plain(
            "*** RS_TransientDiffusion1EXP error\n"
            "    no convergence" ) ;
      }
      result = temp_0 + 2.0*(temp_inf-temp_0)*result ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
RS_TransientDiffusion1EXP:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_TransientDiffusion1EXP:: usage" ) ;

   static std::string result = "TransientDiffusion1($DS_T)" ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
RS_TransientDiffusion1EXP:: valid_arguments(
                              PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_TransientDiffusion1EXP:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result =
      ( some_arguments->count() == 1 ) &&
      ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::Double ) ;
   return( result ) ;
}
