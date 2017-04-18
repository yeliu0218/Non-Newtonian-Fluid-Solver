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

#include <RS_Perturbation.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <doubleVector.hh>
#include <math.h>
#include <iostream>

RS_Perturbation const* 
RS_Perturbation:: PROTOTYPE_U = 
   new RS_Perturbation( "Perturbed_And_Mean_Velocity", U ) ;

RS_Perturbation const* 
RS_Perturbation:: PROTOTYPE_UD = 
   new RS_Perturbation( "Perturbed_grad_velocity", UD ) ;

RS_Perturbation const* 
RS_Perturbation:: PROTOTYPE_P = 
   new RS_Perturbation( "Perturbed_pressure", P ) ;


RS_Perturbation const* 
RS_Perturbation:: PROTOTYPE_RHO = 
   new RS_Perturbation( "Perturbed_density", RHO ) ;

RS_Perturbation const* 
RS_Perturbation:: PROTOTYPE_CC = 
   new RS_Perturbation( "Mean_Concentration", CC ) ;

RS_Perturbation const* 
RS_Perturbation:: PROTOTYPE_Uave = 
   new RS_Perturbation( "Mean_Velocity", Uave ) ;





//----------------------------------------------------------------------
RS_Perturbation:: RS_Perturbation( 
                                                   std::string const& a_name,
                			           Func an_exp ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , EXPR( an_exp )
   , DV_result_1( 1 ) 
   , DV_result_2( 2 ) 
   , doubleArray2D_result( 2, 2 )
{
 PEL_LABEL( "RS_Perturbation" ) ;
}

//----------------------------------------------------------------------
RS_Perturbation*
RS_Perturbation:: create_replica( 
                                PEL_Object* a_owner,
				PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Perturbation:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;
   
   RS_Perturbation* result = new RS_Perturbation( 
                                                      a_owner, 
						      name(), 
						      argument_list, 
						      EXPR ) ;
   
   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
RS_Perturbation:: RS_Perturbation(
                                   PEL_Object* a_owner,
   			           std::string const& a_name,
			           PEL_Sequence const* argument_list,
			           Func an_exp ) 
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , EXPR( an_exp )
   , DV_result_1( 1 ) 
   , DV_result_2( 2 ) 
   , doubleArray2D_result( 2, 2 )
{
   PEL_LABEL( "RS_Perturbation:: RS_Perturbation" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
RS_Perturbation:: ~RS_Perturbation( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Perturbation:: ~RS_Perturbation" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
RS_Perturbation:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Perturbation:: data_type" ) ;

   PEL_Data::Type result = Undefined ;
   switch( EXPR )
   {
      case U : result = DoubleVector ; 
         break ;
      case P : result = DoubleVector ;
	 break ;
      case UD : result = DoubleArray2D ;
         break ;
      case RHO : result= DoubleVector ;
         break ;
      case CC : result = DoubleVector ;
	 break ;
     case Uave : result = DoubleVector ;
	 break ;
   }
   return result ;      
}

//----------------------------------------------------------------------
doubleVector const&
RS_Perturbation:: to_double_vector( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Perturbation:: to_double_vector" ) ;
   PEL_CHECK_PRE( to_double_vector_PRE( ct ) ) ; 

   doubleVector& result = ( ( EXPR == P || EXPR == RHO || EXPR == CC ) ? 
                               DV_result_1 
                            : DV_result_2 ) ;

   doubleVector const& res = arg(0)->to_double_vector(ct) ;
   double ri     = arg(1)->to_double(ct) ;
   double ry = arg(2)->to_double(ct) ;
   double B  = arg(3)->to_double(ct) ;
   double m = arg(4)->to_double(ct) ;
   double alpha  = arg(5)->to_double(ct) ;
   double amplitude  = arg(6)->to_double(ct) ;
   double chanorpipe  = arg(7)->to_double(ct) ;
   if ( chanorpipe==2.0 )
{
   if( EXPR == U )
   {
      double r = res( 0 ) ; 
      double z = res( 1 ) ;
     if (r >=0.0 && r<=ri)
    {
     result( 0 ) = -amplitude*(pow(r,3)*pow((r-ri),2)*alpha*(PEL::cos(alpha*z)));
     result( 1 ) = (B/(2.*ry))*((1./m)*(pow(ri,2)-pow(r,2))+pow((1.0-ry),2))+(4.0*pow(r,2)*pow((r-ri),2)+2.0*pow(r,3)*(r-ri))*(PEL::sin(alpha*z))*amplitude;
    }
    if (r >ri && r<=ry)
    {
     result( 0 ) = 0.0;
     result( 1 ) = (B/(2.*ry))*pow((1.0-ry),2);
    }
      if (r >ry && r<=1.0)
    {
     result( 0 ) = -(pow((r-ry),2)*pow((r-1.0),2)*amplitude*(alpha*(PEL::cos(alpha*z))))/r;
     result( 1 ) = (B/(2.*ry))*(pow((1.0-ry),2)-pow((r-ry),2))+(PEL::sin(alpha*z))*amplitude*(((2.0*(r-1.0)*pow((r-ry),2))/r)+((2.0*(r-ry)*pow((r-1.0),2))/r));
    }
   }
  else if( EXPR == Uave)
   {
      double r = res( 0 ) ; 
      //double z = res( 1 ) ;
     if (r >=0.0 && r<=ri)
    {
     result( 0 ) = 0.0;
     result( 1 ) = (B/(2.*ry))*((1./m)*(pow(ri,2)-pow(r,2))+pow((1.0-ry),2));
    }
    if (r >ri && r<=ry)
    {
     result( 0 ) = 0.0;
     result( 1 ) = (B/(2.*ry))*pow((1.0-ry),2);
    }
      if (r >ry && r<=1.0)
    {
     result( 0 ) = 0.0;
     result( 1 ) = (B/(2.*ry))*(pow((1.0-ry),2)-pow((r-ry),2));
    }
   }
   else if( EXPR == P )
   { 
      
      result( 0 ) = 0.0 ;
   }
   else if( EXPR == RHO )
   { 
     
      result( 0 ) = 0.0 ;
   }
   else if( EXPR == CC )
   { 
     double r = res( 0 ) ;
     if( r >=0.0 && r<=ri )
    {
      result( 0 ) = 0.0 ;
     }
     if( r >ri && r<=1.0 )
    {
      result( 0 ) = 1.0 ;
     }
   }
   return result ;
}

   if ( chanorpipe==1.0 )
{
   if( EXPR == U )
   {
      double r = res( 0 ) ; 
      double z = res( 1 ) ;
     if (r >=-ri && r<=ri)
    {
     
     result( 0 ) = -amplitude*(pow(PEL::abs(r),3)*pow((PEL::abs(r)-ri),2)*alpha*(PEL::cos(alpha*z)));
     result( 1 ) = (B/(2.*ry))*((1./m)*(pow(ri,2)-pow(PEL::abs(r),2))+pow((1.0-ry),2))+(3.0*pow(PEL::abs(r),2)*pow((PEL::abs(r)-ri),2)+2.0*pow(PEL::abs(r),3)*(PEL::abs(r)-ri))*(PEL::sin(alpha*z))*amplitude;
    }
    if ( (r >ri && r<=ry) || (r >=-ry && r<=-ri) )
    {
     result( 0 ) = 0.0;
     result( 1 ) = (B/(2.*ry))*pow((1.0-ry),2);
    }
      if ( (r >ry && r<=1.0) || (r>=-1.0 && r<-ry) )
    {
    
     result( 0 ) = -(pow((PEL::abs(r)-ry),2)*pow((PEL::abs(r)-1.0),2)*amplitude*(alpha*(PEL::cos(alpha*z))));
     result( 1 ) = (B/(2.*ry))*(pow((1.0-ry),2)-pow((PEL::abs(r)-ry),2))+(PEL::sin(alpha*z))*amplitude*(((2.0*(PEL::abs(r)-1.0)*pow((PEL::abs(r)-ry),2)))+((2.0*(PEL::abs(r)-ry)*pow((PEL::abs(r)-1.0),2))));
    }
   }
  else if( EXPR == Uave)
   {
      double r = res( 0 ) ; 
      //double z = res( 1 ) ;
     if (r >=-ri && r<=ri)
    {
     
     result( 0 ) = 0.0;
     result( 1 ) = (B/(2.*ry))*((1./m)*(pow(ri,2)-pow(PEL::abs(r),2))+pow((1.0-ry),2));
    }
    if ( (r >ri && r<=ry) || (r >=-ry && r<=-ri) )
    {
     result( 0 ) = 0.0;
     result( 1 ) = (B/(2.*ry))*pow((1.0-ry),2);
    }
      if ( (r >ry && r<=1.0) || (r>=-1.0 && r<-ry) )
    {
     
     result( 0 ) = 0.0;
     result( 1 ) = (B/(2.*ry))*(pow((1.0-ry),2)-pow((PEL::abs(r)-ry),2));
    }
   }
   else if( EXPR == P )
   { 
      
      result( 0 ) = 0.0 ;
   }
   else if( EXPR == RHO )
   { 
     
      result( 0 ) = 0.0 ;
   }
   else if( EXPR == CC )
   { 
     double r = res( 0 ) ;
     if (r >=-ri && r<=ri)
    {
      result( 0 ) = 0.0 ;
     }
     if( (r >ri && r<=1.0) || (r>=-1.0 && r<-ri) )
    {
      result( 0 ) = 1.0 ;
     }
   }
   return result ;
}
   if ( chanorpipe==3.0 )
{
   if( EXPR == U )
   {
      double r = res( 0 ) ; 
      double z = res( 1 ) ;
     if (r >=0.0 && r<=ri)
    {
     result( 0 ) = -amplitude*(pow(r,3)*pow((r-ri),3)*alpha*(PEL::cos(alpha*z)));
     result( 1 ) = (B/(2.*ry))*((1./m)*(pow(ri,2)-pow(r,2))+pow((1.0-ry),2))+((4.0*pow(r,2)*pow((r-ri),3)+3.0*pow(r,3)*pow((r-ri),2.0))*(PEL::sin(alpha*z))+(4.0*pow(r,2)*pow((r-ri),2)+2.0*pow(r,3)*pow((r-ri),1.0)))*amplitude;
    }
    if (r >ri && r<=ry)
    {
     result( 0 ) =-(pow((r-ri),3)*pow((r-1.0),2)*amplitude*(alpha*(PEL::cos(alpha*z))))/r;
     result( 1 ) = (B/(2.*ry))*pow((1.0-ry),2)+amplitude*(((2*(r-1.0)*pow((r-ri),3)+3.0*pow((r-1.0),2)*pow((r-ri),2.0))*(PEL::sin(alpha*z)))+((2.0*(r-1.0)*pow((r-ri),2)+2.0*pow((r-1.0),2)*(r-ri)))*((2.0*m*pow(ri,4)-(B*pow(ri,2))/ry-B*ri)/(2*pow((ri-1.0),2))))/r;
    }
      if (r >ry && r<=1.0)
    {
     result( 0 ) = -(pow((r-ri),3)*pow((r-1.0),2)*amplitude*(alpha*(PEL::cos(alpha*z))))/r;
     result( 1 ) = (B/(2.*ry))*(pow((1.0-ry),2)-pow((r-ry),2))+amplitude*(((2*(r-1.0)*pow((r-ri),3)+3.0*pow((r-1.0),2)*pow((r-ri),2.0))*(PEL::sin(alpha*z)))+((2.0*(r-1.0)*pow((r-ri),2)+2.0*pow((r-1.0),2)*(r-ri)))*((2.0*m*pow(ri,4)-(B*pow(ri,2))/ry-B*ri)/(2*pow((ri-1.0),2))))/r;
    }
   }
  else if( EXPR == Uave)
   {
      double r = res( 0 ) ; 
      //double z = res( 1 ) ;
     if (r >=0.0 && r<=ri)
    {
     result( 0 ) = 0.0;
     result( 1 ) = (B/(2.*ry))*((1./m)*(pow(ri,2)-pow(r,2))+pow((1.0-ry),2));
    }
    if (r >ri && r<=ry)
    {
     result( 0 ) = 0.0;
     result( 1 ) = (B/(2.*ry))*pow((1.0-ry),2);
    }
      if (r >ry && r<=1.0)
    {
     result( 0 ) = 0.0;
     result( 1 ) = (B/(2.*ry))*(pow((1.0-ry),2)-pow((r-ry),2));
    }
   }
   else if( EXPR == P )
   { 
      
      result( 0 ) = 0.0 ;
   }
   else if( EXPR == RHO )
   { 
     
      result( 0 ) = 0.0 ;
   }
   else if( EXPR == CC )
   { 
     double r = res( 0 ) ;
     if( r >=0.0 && r<=ri )
    {
      result( 0 ) = 0.0 ;
     }
     if( r >ri && r<=1.0 )
    {
      result( 0 ) = 1.0 ;
     }
   }
   return result ;
}
 if ( chanorpipe==4.0 )
 {
   if( EXPR == U )
   {
      double r = res( 0 ) ; 
      double z = res( 1 ) ;
     if (r >=-ri && r<=ri)
    {
     
     result( 0 ) = -amplitude*(pow(PEL::abs(r),3)*pow((PEL::abs(r)-ri),3)*alpha*(PEL::cos(alpha*z)));
     result( 1 ) = (B/(2.*ry))*((1./m)*(pow(ri,2)-pow(PEL::abs(r),2))+pow((1.0-ry),2))+((3.0*pow(PEL::abs(r),2)*pow((PEL::abs(r)-ri),3)+3.0*pow(PEL::abs(r),3)*pow((PEL::abs(r)-ri),2))*(PEL::sin(alpha*z))+(3.0*pow(PEL::abs(r),2)*pow((PEL::abs(r)-ri),2)+2.0*pow(PEL::abs(r),3)*(PEL::abs(r)-ri)))*amplitude;
    }
    if ((r >ri && r<=ry) || (r >=-ry && r<=-ri) )
    {
     result( 0 ) = -(pow((PEL::abs(r)-ri),3)*pow((PEL::abs(r)-1.0),2)*amplitude*(alpha*(PEL::cos(alpha*z))));
     result( 1 ) = (B/(2.*ry))*pow((1.0-ry),2)+amplitude*((3.0*pow((PEL::abs(r)-ri),2)*pow((PEL::abs(r)-1.0),2)+2.0*pow((PEL::abs(r)-ri),3)*pow((PEL::abs(r)-1.0),1))*(PEL::sin(alpha*z))+((2.0*pow((PEL::abs(r)-ri),2)*pow((PEL::abs(r)-1.0),1)+2.0*pow((PEL::abs(r)-ri),1)*pow((PEL::abs(r)-1.0),2))*((2*m*pow(ri,3.0)-B*ri/ry-B)/(2*pow((ri-1.0),2.0)))));
    }
      if ((r >ry && r<=1.0) || (r>=-1.0 && r<-ry))
    {
    
     result( 0 ) = -(pow((PEL::abs(r)-ri),3)*pow((PEL::abs(r)-1.0),2)*amplitude*(alpha*(PEL::cos(alpha*z))));
     result( 1 ) = (B/(2.*ry))*(pow((1.0-ry),2)-pow((PEL::abs(r)-ry),2))+amplitude*((3.0*pow((PEL::abs(r)-ri),2)*pow((PEL::abs(r)-1.0),2)+2.0*pow((PEL::abs(r)-ri),3)*pow((PEL::abs(r)-1.0),1))*(PEL::sin(alpha*z))+((2.0*pow((PEL::abs(r)-ri),2)*pow((PEL::abs(r)-1.0),1)+2.0*pow((PEL::abs(r)-ri),1)*pow((PEL::abs(r)-1.0),2))*((2*m*pow(ri,3.0)-B*ri/ry-B)/(2*pow((ri-1.0),2.0)))));
    }
   }
  else if( EXPR == Uave)
   {
      double r = res( 0 ) ; 
      //double z = res( 1 ) ;
     if (r >=-ri && r<=ri)
    {    
     result( 0 ) = 0.0;
     result( 1 ) = (B/(2.*ry))*((1./m)*(pow(ri,2)-pow(PEL::abs(r),2))+pow((1.0-ry),2));
    }
    if ((r >ri && r<=ry) || (r >=-ry && r<=-ri) )
    {
     result( 0 ) = 0.0;
     result( 1 ) = (B/(2.*ry))*pow((1.0-ry),2);
    }
      if ((r >ry && r<=1.0) || (r>=-1.0 && r<-ry))
    {    
     result( 0 ) = 0.0;
     result( 1 ) = (B/(2.*ry))*(pow((1.0-ry),2)-pow((PEL::abs(r)-ry),2));
    }
   }
   else if( EXPR == P )
   {       
      result( 0 ) = 0.0 ;
   }
   else if( EXPR == RHO )
   {      
      result( 0 ) = 0.0 ;
   }
   else if( EXPR == CC )
   { 
     double r = res( 0 ) ;
     if (r >=-ri && r<=ri)
     {
      result( 0 ) = 0.0 ;
     }
     if( (r >ri && r<=1.0) || (r>=-1.0 && r<-ri) )
     {
      result( 0 ) = 1.0 ;
     }
   }
   return result ;
}
 // Should give an error if chanorpipe isn't 1., 3. or 4.
 // Better: if( chanorpipe==1.0 )
 //         else if( chanorpipe==3.0 )
 //         else if( chanorpipe==4.0 )
 //         else ERROR
 // and then RETURN statement (better would be only one return statement)
 return result; // is only reached if chanorpipe isn't 1., 3. or 4.
}

//----------------------------------------------------------------------
doubleArray2D  const&
RS_Perturbation:: to_double_array2D( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Perturbation:: to_double_array2D" ) ;
   PEL_CHECK_PRE( to_double_array2D_PRE( ct ) ) ;

   doubleArray2D& result = doubleArray2D_result ;
   //doubleVector const& res = arg(0)->to_double_vector(ct) ;
   //double ri     = arg(1)->to_double(ct) ;
   //double ry = arg(2)->to_double(ct) ;
   //double B  = arg(3)->to_double(ct) ;
   //double m = arg(4)->to_double(ct) ;
   //double alpha  = arg(5)->to_double(ct) ;
   //double amplitude  = arg(6)->to_double(ct) ;
   //double chanorpipe  = arg(7)->to_double(ct) ;
   
   if( EXPR == UD )
   {
      result( 0, 0 ) = 0.0 ;
      result( 0, 1 ) = 0.0 ;
      result( 1, 0 ) = 0.0 ;
      result( 1, 1 ) = 0.0 ;
   }
   
   return result ;
}


//----------------------------------------------------------------------
std::string const& 
RS_Perturbation:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Perturbation:: usage" ) ;

   static std::string result ;
   switch( EXPR )
   {
      case U : result = "Perturbed_And_Mean_Velocity($DV_X,$DS_ri,$DS_ry,$DS_Bn2,$DS_kappa1,$DS_ALPHA1,$DS_AMPLITUDE,$DS_CHANORPIPE)"; 
         break ;
      case P : result = "Perturbed_pressure($DV_X,$DS_ri,$DS_ry,$DS_Bn2,$DS_kappa1,$DS_ALPHA1,$DS_AMPLITUDE,$DS_CHANORPIPE)";
	 break ;
      case RHO : result = "Perturbed_density($DV_X,$DS_ri,$DS_ry,$DS_Bn2,$DS_kappa1,$DS_ALPHA1,$DS_AMPLITUDE,$DS_CHANORPIPE)";
	 break ;
      case UD : result = "Perturbed_grad_velocity($DV_X,$DS_ri,$DS_ry,$DS_Bn2,$DS_kappa1,$DS_ALPHA1,$DS_AMPLITUDE,$DS_CHANORPIPE)";
         break ;
      case CC : result = "Mean_Concentration($DV_X,$DS_ri,$DS_ry,$DS_Bn2,$DS_kappa1,$DS_ALPHA1,$DS_AMPLITUDE,$DS_CHANORPIPE)";
	 break ;
      case Uave : result = "Mean_Velocity($DV_X,$DS_ri,$DS_ry,$DS_Bn2,$DS_kappa1,$DS_ALPHA1,$DS_AMPLITUDE,$DS_CHANORPIPE)";
	 break ;
        }   
   return result ;
}

//----------------------------------------------------------------------
bool
RS_Perturbation:: valid_arguments( 
                                   PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Perturbation:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result =  
 ( some_arguments->count() == 8 ) &&
 ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::DoubleVector ) &&
 ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 3 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 4 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 5 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 6 )->data_type() == PEL_Data::Double ) &&
 ( extract_arg( some_arguments, 7 )->data_type() == PEL_Data::Double ) ;
   return result ;
}
