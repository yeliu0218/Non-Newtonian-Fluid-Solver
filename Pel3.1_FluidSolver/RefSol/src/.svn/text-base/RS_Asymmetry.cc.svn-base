/*
 *  Copyright : 
 *    "Institut de Radioprotection et de S�ret� Nucl�aire - IRSN" (1995-2008)
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

#include <RS_Asymmetry.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <doubleVector.hh>
#include <math.h>
#include <iostream>

RS_Asymmetry const* 
RS_Asymmetry:: PROTOTYPE_U = 
   new RS_Asymmetry( "Asym_Velocity_BC", U ) ;


RS_Asymmetry const* 
RS_Asymmetry:: PROTOTYPE_CC = 
   new RS_Asymmetry( "Asym_Concentration_BC", CC ) ;




//----------------------------------------------------------------------
RS_Asymmetry:: RS_Asymmetry( 
                                                   std::string const& a_name,
                			           Func an_exp ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , EXPR( an_exp )
   , DV_result_1( 1 ) 
   , DV_result_2( 2 ) 
   
{
 PEL_LABEL( "RS_Asymmetry" ) ;
 //PEL::out() << "Constructor" << std::endl;
 //PEL::out() << DV_result_1.size() << std::endl;

}

//----------------------------------------------------------------------
RS_Asymmetry*
RS_Asymmetry:: create_replica( 
                                PEL_Object* a_owner,
				PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Asymmetry:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   RS_Asymmetry* result = new RS_Asymmetry( 
                                                      a_owner, 
						      name(), 
						      argument_list, 
						      EXPR ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
RS_Asymmetry:: RS_Asymmetry(
                                   PEL_Object* a_owner,
   			           std::string const& a_name,
			           PEL_Sequence const* argument_list,
			           Func an_exp ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , EXPR( an_exp )
   , DV_result_1( 1 ) 
   , DV_result_2( 2 ) 
  
{
   PEL_LABEL( "RS_Asymmetry:: RS_Asymmetry" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
RS_Asymmetry:: ~RS_Asymmetry( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Asymmetry:: ~RS_Asymmetry" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
RS_Asymmetry:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Asymmetry:: data_type" ) ;

   PEL_Data::Type result = Undefined ;
   switch( EXPR )
   {
      case U : result = DoubleVector ; 
         break ;
      case CC : result = DoubleVector ;
	 break ;
    }
   return result ;      
}

//----------------------------------------------------------------------
doubleVector const&
RS_Asymmetry:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Asymmetry:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ; 

//	PEL::out() << "in RS_Move_B START" << std::endl;
   doubleVector& result = ( ( EXPR == CC ) ? 
                               DV_result_1 
                            : DV_result_2 ) ;

   doubleVector const& res = arg(0)->to_double_vector(ct) ;
	//PEL::out() << "in RS_Move_B - 1" << std::endl;
   double Rl = arg(1)->to_double(ct) ;
   double Rr = arg(2)->to_double(ct) ;
   double Um = arg(3)->to_double(ct) ;
   double Ul = arg(4)->to_double(ct) ;
   double Ur = arg(5)->to_double(ct) ;
   double Vm = arg(6)->to_double(ct) ;
   double Vl = arg(7)->to_double(ct) ;
   double Vr = arg(8)->to_double(ct) ;
   double Cm = arg(9)->to_double(ct) ;
   double Cl = arg(10)->to_double(ct) ;
   double Cr = arg(11)->to_double(ct) ;
	//PEL::out() << "in RS_Move_B - 2" << std::endl;
   if( EXPR == U )
   {
	//PEL::out() << "*** change U ***" << std::endl;
    double r = res( 0 ) ; 
    //double z = res( 1 ) ;
   
      
     if (r >=-1.0 && r<=Rl)
    {
     result( 0 ) = Vl;
     result( 1 ) = Ul;
    }
    if (r >Rl&& r<Rr)
    {
     result( 0 ) = Vm;
     result( 1 ) = Um;
    }
      if (r >=Rr&& r<=1.0)
    {
     result( 0 ) = Vr;
     result( 1 ) = Ur;
    }
  
   }

 else if ( EXPR == CC )
   { 
	//PEL::out() << "*** change C ***" << std::endl;
     double r = res( 0 ) ;
   
    
     if( r >=-1.0 && r<=Rl )
    {
      result( 0 ) = Cl ;
     }
     if( r >Rl && r<Rr )
    {
      result( 0 ) = Cm ;
     }
     if( r >=Rr && r<=1.0 )
    {
      result( 0 ) = Cr ;
    
    }
   }
	//else	{
	//PEL::out() << "EXPR =" << EXPR << std::endl;
	//PEL::out() << "NOT IMPLEMENTED " << std::endl;
	//}
   return result ;


}



//----------------------------------------------------------------------
std::string const& 
RS_Asymmetry:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Asymmetry:: usage" ) ;

   static std::string result ;
   switch( EXPR )
   {
      case U : result = "Asym_Velocity_BC($DV_X,$DS_Rl,$DS_Rr,$DS_Um,$DS_Ul,$DS_Ur,$DS_Vm,$DS_Vl,$DS_Vr,$DS_Cm,$DS_Cl,$DS_Cr)"; 
         break ;
      
      case CC : result = "Asym_Concentration_BC($DV_X,$DS_Rl,$DS_Rr,$DS_Um,$DS_Ul,$DS_Ur,$DS_Vm,$DS_Vl,$DS_Vr,$DS_Cm,$DS_Cl,$DS_Cr)";
	 break ;
     
        }   
   return result ;
}

//----------------------------------------------------------------------
bool
RS_Asymmetry:: valid_arguments( 
                                   PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Asymmetry:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result =  
 ( some_arguments->count() == 12 ) &&
 ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::DoubleVector ) &&
 ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 3 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 4 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 5 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 6 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 7 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 8 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 9 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 10 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 11 )->data_type() == PEL_Data::Double ) ;
   return result ;
}
