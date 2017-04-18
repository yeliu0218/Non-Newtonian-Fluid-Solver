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

#include <RS_Multi_Layer5.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <doubleVector.hh>
#include <math.h>
#include <iostream>

RS_Multi_Layer5 const* 
RS_Multi_Layer5:: PROTOTYPE_U = 
   new RS_Multi_Layer5( "Multi_Layer_Velocity_BC5", U ) ;


RS_Multi_Layer5 const* 
RS_Multi_Layer5:: PROTOTYPE_CC = 
   new RS_Multi_Layer5( "Multi_Layer_Concentration_BC5", CC ) ;




//----------------------------------------------------------------------
RS_Multi_Layer5:: RS_Multi_Layer5( 
                                                   std::string const& a_name,
                			           Func an_exp ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , EXPR( an_exp )
   , DV_result_1( 1 ) 
   , DV_result_2( 2 ) 
   
{
 PEL_LABEL( "RS_Multi_Layer5" ) ;
 //PEL::out() << "Constructor" << std::endl;
 //PEL::out() << DV_result_1.size() << std::endl;

}

//----------------------------------------------------------------------
RS_Multi_Layer5*
RS_Multi_Layer5:: create_replica( 
                                PEL_Object* a_owner,
				PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Multi_Layer5:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   RS_Multi_Layer5* result = new RS_Multi_Layer5( 
                                                      a_owner, 
						      name(), 
						      argument_list, 
						      EXPR ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
RS_Multi_Layer5:: RS_Multi_Layer5(
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
   PEL_LABEL( "RS_Multi_Layer5:: RS_Multi_Layer5" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
RS_Multi_Layer5:: ~RS_Multi_Layer5( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Multi_Layer5:: ~RS_Multi_Layer5" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
RS_Multi_Layer5:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Multi_Layer5:: data_type" ) ;

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
RS_Multi_Layer5:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Multi_Layer5:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ; 

//	PEL::out() << "in RS_Move_B START" << std::endl;
   doubleVector& result = ( ( EXPR == CC ) ? 
                               DV_result_1 
                            : DV_result_2 ) ;

   doubleVector const& res = arg(0)->to_double_vector(ct) ;
	//PEL::out() << "in RS_Move_B - 1" << std::endl;
   //double lwlf     = arg(1)->to_double(ct) ;
   double rwlf     = arg(1)->to_double(ct) ;
   double lwmf     = arg(2)->to_double(ct) ;
   double rwmf     = arg(3)->to_double(ct) ;
   double lwrf     = arg(4)->to_double(ct) ;
  // double rwrf     = arg(6)->to_double(ct) ;
   double Um0      = arg(5)->to_double(ct) ;
   double Ul0      = arg(6)->to_double(ct) ;
   double Ur0      = arg(7)->to_double(ct) ;
   double Um1      = arg(8)->to_double(ct) ;
   double Ul1      = arg(9)->to_double(ct) ;
   double Ul2      = arg(10)->to_double(ct) ;
   //double Ul3      = arg(11)->to_double(ct) ;
   double Ur1      = arg(11)->to_double(ct) ;
   double Ur2      = arg(12)->to_double(ct) ;
  // double Ur3      = arg(14)->to_double(ct) ;
   double CN       = arg(13)->to_double(ct) ;
   double CB       = arg(14)->to_double(ct) ;
   double T        = arg(15)->to_double(ct) ;
   double T1         = arg(16)->to_double(ct) ;
  
 
	//PEL::out() << "in RS_Move_B - 2" << std::endl;
   if( EXPR == U )
   {
	//PEL::out() << "*** change U ***" << std::endl;
    double r = res( 0 ) ; 
    //double z = res( 1 ) ;
    if ( T<=T1)
    {
      
     if (r >=-1.0 && r<=lwmf)
    {
     result( 0 ) = 0.0;
     result( 1 ) = Ul0;
    }
    if (r >lwmf && r<rwmf)
    {
     result( 0 ) = 0.0;
     result( 1 ) = Um0;
    }
      if (r >=rwmf && r<=1.0)
    {
     result( 0 ) = 0.0;
     result( 1 ) = Ur0;
    }
   }
   if ( T>T1 && T<(2.0*T1))
   //if ( T>T1)
    {
      
     //if (r >=-1.0 && r<=lwlf)
   // {
    // result( 0 ) = 0.0;
   //  result( 1 ) = Ul3;
  //  }
    if (r >=-1.0 && r<rwlf)
    {
     result( 0 ) = 0.0;
     result( 1 ) = Ul2;
    }
    if (r >=rwlf && r<=lwmf)
    {
     result( 0 ) = 0.0;
     result( 1 ) = Ul1;
    }
    if (r >lwmf && r<rwmf)
    {
     result( 0 ) = 0.0;
     result( 1 ) = Um1;
    }
    if (r >=rwmf && r<=lwrf)
    {
     result( 0 ) = 0.0;
     result( 1 ) = Ur1;
    }
    if (r >lwrf && r<=1.0)
    {
     result( 0 ) = 0.0;
     result( 1 ) = Ur2;
    }
    //if (r >=rwrf && r<=1.0)
    //{
    // result( 0 ) = 0.0;
   //  result( 1 ) = Ur3;
   // }
   }
   if ( T>=(3.0*T1))
   {
      
    if (r >=-1.0 && r<=lwmf)
    {
     result( 0 ) = 0.0;
     result( 1 ) = Ul0;
    }
    if (r >lwmf && r<rwmf)
    {
     result( 0 ) = 0.0;
     result( 1 ) = Um0;
    }
      if (r >=rwmf && r<=1.0)
    {
     result( 0 ) = 0.0;
     result( 1 ) = Ur0;
    }
   }
   } else if ( EXPR == CC )
   { 
	//PEL::out() << "*** change C ***" << std::endl;
     double r = res( 0 ) ;
   
     if ( T<=T1)
    {
      
     if (r >=-1.0 && r<=lwmf)
    {
     result( 0 ) = CB;
     
    }
    if (r >lwmf && r<rwmf)
    {
     result( 0 ) = CN;
     
    }
      if (r >=rwmf && r<=1.0)
    {
     result( 0 ) = CB;
     
    }
   }
   if ( T>T1 && T<(2.0*T1))
    {
      
     //if (r >=-1.0 && r<=lwlf)
   // {
   //  result( 0 ) = CB;
    
  //  }
    if (r >=-1.0 && r<rwlf)
    {
     result( 0 ) = CN;
     
    }
    if (r >=rwlf && r<=lwmf)
    {
     result( 0 ) = CB;
     
    }
    if (r >lwmf && r<rwmf)
    {
     result( 0 ) = CN;
     
    }
    if (r >=rwmf && r<=lwrf)
    {
     result( 0 ) = CB;
     
    }
    if (r >lwrf && r<=1.0)
    {
     result( 0 ) = CN;
     
    }
   // if (r >=rwrf && r<=1.0)
    //{
    // result( 0 ) = CB;
     
  //  }
   }
   if ( T>=(2.0*T1))
    {
      
    if (r >=-1.0 && r<=lwmf)
    {
     result( 0 ) = CB;
     
    }
    if (r >lwmf && r<rwmf)
    {
     result( 0 ) = CN;
    
    }
      if (r >=rwmf && r<=1.0)
    {
     result( 0 ) = CB;
     
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
RS_Multi_Layer5:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Multi_Layer5:: usage" ) ;

   static std::string result ;
   switch( EXPR )
   {
      case U : result = "Multi_Layer_Velocity_BC5($DV_X,$DS_rwlf,$DS_lwmf,$DS_rwmf,$DS_lwrf,$DS_Um0,$DS_Ul0,$DS_Ur0,$DS_Um1,$DS_Ul1,$DS_Ul2,$DS_Ur1,$DS_Ur2,$DS_CN,$DS_CB,$DS_T,$DS_T1)"; 
         break ;
      
      case CC : result = "Multi_Layer_Concentration_BC5($DV_X,$DS_lwlf,$DS_rwlf,$DS_lwmf,$DS_rwmf,$DS_lwrf,$DS_rwrf,$DS_Um0,$DS_Ul0,$DS_Ur0,$DS_Um1,$DS_Ul1,$DS_Ul2,$DS_Ur1,$DS_Ur2,$DS_CN,$DS_CB,$DS_T,$DS_T1)";
	 break ;
     
        }   
   return result ;
}

//----------------------------------------------------------------------
bool
RS_Multi_Layer5:: valid_arguments( 
                                   PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Multi_Layer5:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result =  
 ( some_arguments->count() == 17 ) &&
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
 ( extract_arg( some_arguments, 16 )->data_type() == PEL_Data::Double ) ;
// ( extract_arg( some_arguments, 17 )->data_type() == PEL_Data::Double ) &&
 //( extract_arg( some_arguments, 18 )->data_type() == PEL_Data::Double ) ;
   return result ;
}
