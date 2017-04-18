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

#include <RS_Move_BC.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <doubleVector.hh>
#include <math.h>
#include <iostream>

RS_Move_BC const* 
RS_Move_BC:: PROTOTYPE_U = 
   new RS_Move_BC( "Displace_Velocity_BC", U ) ;


RS_Move_BC const* 
RS_Move_BC:: PROTOTYPE_CC = 
   new RS_Move_BC( "Displace_Concentration_BC", CC ) ;




//----------------------------------------------------------------------
RS_Move_BC:: RS_Move_BC( 
                                                   std::string const& a_name,
                			           Func an_exp ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , EXPR( an_exp )
   , DV_result_1( 1 ) 
   , DV_result_2( 2 ) 
   
{
 PEL_LABEL( "RS_Move_BC" ) ;
 //PEL::out() << "Constructor" << std::endl;
 //PEL::out() << DV_result_1.size() << std::endl;

}

//----------------------------------------------------------------------
RS_Move_BC*
RS_Move_BC:: create_replica( 
                                PEL_Object* a_owner,
				PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Move_BC:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   RS_Move_BC* result = new RS_Move_BC( 
                                                      a_owner, 
						      name(), 
						      argument_list, 
						      EXPR ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
RS_Move_BC:: RS_Move_BC(
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
   PEL_LABEL( "RS_Move_BC:: RS_Move_BC" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
RS_Move_BC:: ~RS_Move_BC( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Move_BC:: ~RS_Move_BC" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
RS_Move_BC:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Move_BC:: data_type" ) ;

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
RS_Move_BC:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Move_BC:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ; 

//	PEL::out() << "in RS_Move_B START" << std::endl;
   doubleVector& result = ( ( EXPR == CC ) ? 
                               DV_result_1 
                            : DV_result_2 ) ;

   doubleVector const& res = arg(0)->to_double_vector(ct) ;
	//PEL::out() << "in RS_Move_B - 1" << std::endl;
   double Rl1     = arg(1)->to_double(ct) ;
   double Rr1     = arg(2)->to_double(ct) ;
   double Rl2     = arg(3)->to_double(ct) ;
   double Rr2     = arg(4)->to_double(ct) ;
   double Um1  = arg(5)->to_double(ct) ;
   double Ul1 = arg(6)->to_double(ct) ;
   double Ur1  = arg(7)->to_double(ct) ;
   double Um2  = arg(8)->to_double(ct) ;
   double Ul2 = arg(9)->to_double(ct) ;
   double Ur2  = arg(10)->to_double(ct) ;
   double T  = arg(11)->to_double(ct) ;
   double T1 = arg(12)->to_double(ct) ;
   double Vm1  = arg(13)->to_double(ct) ;
   double Vl1 = arg(14)->to_double(ct) ;
   double Vr1  = arg(15)->to_double(ct) ;
   double Vm2  = arg(16)->to_double(ct) ;
   double Vl2 = arg(17)->to_double(ct) ;
   double Vr2  = arg(18)->to_double(ct) ;
   double Cl  = arg(19)->to_double(ct) ;
   double Cm = arg(20)->to_double(ct) ;
   double Cr  = arg(21)->to_double(ct) ;
 
	//PEL::out() << "in RS_Move_B - 2" << std::endl;
   if( EXPR == U )
   {
	//PEL::out() << "*** change U ***" << std::endl;
    double r = res( 0 ) ; 
    //double z = res( 1 ) ;
    if ( T<=T1)
    {
      
     if (r >=-1.0 && r<=Rl1)
    {
     result( 0 ) = Vl1;
     result( 1 ) = Ul1;
    }
    if (r >Rl1 && r<Rr1)
    {
     result( 0 ) = Vm1;
     result( 1 ) = Um1;
    }
      if (r >=Rr1 && r<=1.0)
    {
     result( 0 ) = Vr1;
     result( 1 ) = Ur1;
    }
   }
   if ( T>T1 && T<(2.0*T1))
    {
      
     if (r >=-1.0 && r<=Rl2)
    {
     result( 0 ) = Vl2;
     result( 1 ) = Ul2;
    }
    if (r >Rl2 && r<Rr2)
    {
     result( 0 ) = Vm2;
     result( 1 ) = Um2;
    }
      if (r >=Rr2 && r<=1.0)
    {
     result( 0 ) = Vr2;
     result( 1 ) = Ur2;
    }
   }
   if ( T>=(2.0*T1))
    {
      
     if (r >=-1.0 && r<=Rl1)
    {
     result( 0 ) = Vl1;
     result( 1 ) = Ul1;
    }
    if (r >Rl1 && r<Rr1)
    {
     result( 0 ) = Vm1;
     result( 1 ) = Um1;
    }
      if (r >=Rr1 && r<=1.0)
    {
     result( 0 ) = Vr1;
     result( 1 ) = Ur1;
    }
   }
   } else if ( EXPR == CC )
   { 
	//PEL::out() << "*** change C ***" << std::endl;
     double r = res( 0 ) ;
   
     if (T<=T1)
     {
     if( r >=-1.0 && r<=Rl1 )
    {
      result( 0 ) = Cl ;
     }
     if( r >Rl1 && r<Rr1 )
    {
      result( 0 ) = Cm ;
     }
     if( r >=Rr1 && r<=1.0 )
    {
      result( 0 ) = Cr ;
     }
    }
       if ( T>T1 && T<(2.0*T1))
     {
     if( r >=-1.0 && r<=Rl2 )
    {
      result( 0 ) = Cl ;
     }
     if( r >Rl2 && r<Rr2 )
    {
      result( 0 ) = Cm ;
     }
     if( r >=Rr2 && r<=1.0 )
    {
      result( 0 ) = Cr ;
     }
    }
     if (T>=(2.0*T1))
     {
     if( r >=-1.0 && r<=Rl1 )
    {
      result( 0 ) = Cl ;
     }
     if( r >Rl1 && r<Rr1 )
    {
      result( 0 ) = Cm ;
     }
     if( r >=Rr1 && r<=1.0 )
    {
      result( 0 ) = Cr ;
     }
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
RS_Move_BC:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Move_BC:: usage" ) ;

   static std::string result ;
   switch( EXPR )
   {
      case U : result = "Displace_Velocity_BC($DV_X,$DS_Rl1,$DS_Rr1,$DS_Rl2,$DS_Rr2,$DS_Um1,$DS_Ul1,$DS_Ur1,$DS_Um2,$DS_Ul2,$DS_Ur2,$DS_T,$DS_T1,$DS_Vm1,$DS_Vl1,$DS_Vr1,$DS_Vm2,$DS_Vl2,$DS_Vr2,$DS_Cl,$DS_Cm,$DS_Cr)"; 
         break ;
      
      case CC : result = "Displace_Concentration_BC($DV_X,$DS_Rl1,$DS_Rr1,$DS_Rl2,$DS_Rr2,$DS_Um1,$DS_Ul1,$DS_Ur1,$DS_Um2,$DS_Ul2,$DS_Ur2,$DS_T,$DS_T1,$DS_Vm1,$DS_Vl1,$DS_Vr1,$DS_Vm2,$DS_Vl2,$DS_Vr2,$DS_Cl,$DS_Cm,$DS_Cr)";
	 break ;
     
        }   
   return result ;
}

//----------------------------------------------------------------------
bool
RS_Move_BC:: valid_arguments( 
                                   PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Move_BC:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result =  
 ( some_arguments->count() == 22 ) &&
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
 ( extract_arg( some_arguments, 11 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 12 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 13 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 14 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 15 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 16 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 17 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 18 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 19 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 20 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 21 )->data_type() == PEL_Data::Double ) ;
   return result ;
}
