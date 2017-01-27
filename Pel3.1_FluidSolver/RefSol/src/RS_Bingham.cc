/*
 *  Copyright :
 *    "Institut de Radioprotection et de Sret�Nucl�ire - IRSN" (1995-2008)
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

#include <RS_Bingham.hh>

#include <PEL_assertions.hh>
#include <PEL_Error.hh>
#include <PEL_Sequence.hh>
#include <PEL.hh>

#include <math.h>
#include <iostream>

RS_Bingham const*
RS_Bingham::PROTOTYPE_ONE
    = new RS_Bingham( "Bingham_YieldSurface", ONE ) ;

RS_Bingham const*
RS_Bingham::PROTOTYPE_TWO
    = new RS_Bingham( "Bingham_VelocityProfile", TWO ) ;

RS_Bingham const*
	RS_Bingham::PROTOTYPE_THREE
	= new RS_Bingham( "Newtonian_VelocityProfile", THREE ) ;

RS_Bingham const*
	RS_Bingham::PROTOTYPE_FOUR
	= new RS_Bingham( "Bingham_PipeVelocityProfile", FOUR ) ;

RS_Bingham const*
	RS_Bingham::PROTOTYPE_Five
	= new RS_Bingham( "Transient_Newtonian_VelocityProfile", Five ) ;

//----------------------------------------------------------------------
RS_Bingham:: RS_Bingham(
                               std::string const& a_name, Func an_expr )
//----------------------------------------------------------------------
   : PEL_Expression( a_name )
   , EXPR( an_expr )
{
   PEL_LABEL( "RS_Bingham" ) ;
}

//----------------------------------------------------------------------
RS_Bingham*
RS_Bingham:: create_replica(
          PEL_Object* a_owner, PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Bingham:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   RS_Bingham* result =
               new RS_Bingham( a_owner, name(), argument_list, EXPR ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
RS_Bingham:: RS_Bingham(
                 PEL_Object* a_owner, std::string const& a_name,
		 PEL_Sequence const* argument_list, Func an_expr  )
//----------------------------------------------------------------------
   : PEL_Expression( a_owner, a_name, argument_list )
   , EXPR( an_expr )
{
   PEL_LABEL( "RS_Bingham:: RS_Bingham" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
RS_Bingham:: ~RS_Bingham( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Bingham:: ~RS_Bingham" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_Data::Type
RS_Bingham:: data_type( void ) const
//----------------------------------------------------------------------
{
   return( Double ) ;
}

double ff(double ys, double Bn, double Q) {
  return( pow(ys,3) - 3. * ys * ( 1 + 2.*Q/Bn) + 2.) ;
}

double dff(double ys, double Bn, double Q){
  return( 3.*pow(ys,2) - 3. * ( 1 + 2.*Q/Bn) );
}

double ff_pipe(double ys, double Bn, double Q) {
   return( pow(ys,4) - 4. * ys * ( 1 + 3.*Q/Bn) + 3.) ;
}

double dff_pipe(double ys, double Bn, double Q){
   return( 4.*pow(ys,3) - 4. * ( 1 + 3.*Q/Bn) );
}


//----------------------------------------------------------------------
double
RS_Bingham:: to_double( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Bingham:: to_double" ) ;
   PEL_CHECK_PRE( to_double_PRE( ct ) ) ;

   double ys=0.;
   int iter=1, maxiter=100;
   double dx = 1;

   double Left=-1., Right=1., wall=1.;

   double Bn = 0.0, Bn1=0.0;
   double Q = 1.0;
   switch( EXPR )
   {
	case ONE:
      case TWO:
    	  Bn = arg(0)->to_double(ct) ;

	   switch( EXPR )
            {
		   case	ONE:
			   Q = arg(1)->to_double(ct) ;
		   break;
		   case TWO:
		     Left  = arg(3)->to_double(ct) ;
		     Right = arg(4)->to_double(ct) ;
		     Q = arg(5)->to_double(ct) ;

		     wall = (Right-Left)/2.;
		   break;
		   case THREE:
		   case FOUR:
		   case Five:
		   break;
            }

            Bn1 = Bn * pow(wall,2);

            ys = .5;
            if(Bn==0.0){
            	ys=0.0;
            }
            else {
            	while( fabs(dx) > 1.E-8 && iter < maxiter){
            		dx = ff(ys,Bn1,Q)/dff(ys,Bn1,Q);
            		ys -= dx;
            		iter++;
            	}            
                if( iter >= maxiter)
                  PEL_Error::object()->raise_plain( "Newton method failed: Too many iterations!" ) ;

                if( ys <= 0.0 || ys > 1.0 )
                  {
                    PEL::out() << "NEED ys>0 && ys<1 BUT ys = " << ys << std::endl;
                    PEL_Error::object()->raise_plain( "Newton method failed" ) ;
                  }
            }
            break;
	case THREE :
	   Bn = 0.0;
	break;
	case FOUR:
            Bn = arg(0)->to_double(ct) ;
            Left  = arg(3)->to_double(ct) ;
            Right = arg(4)->to_double(ct) ;
            Right = arg(5)->to_double(ct) ;

            wall = (Right-Left)/2.;
            Bn1 = Bn * pow(wall,3);

            ys = .5;
		while( fabs(dx) > 1.E-8 && iter < maxiter){
		   dx = ff_pipe(ys,Bn1,Q)/dff_pipe(ys,Bn1,Q);
		   ys -= dx;
		   iter++;
		}
		if( iter >= maxiter)
		   PEL_Error::object()->raise_plain( "Newton method failed: Too many iterations!" ) ;

		if( ys <= 0.0 || ys > 1.0 )
		{
		   PEL::out() << "NEED ys>0 && ys<1 BUT ys = " << ys << std::endl;
		   PEL_Error::object()->raise_plain( "Newton method failed" ) ;
		}
	break;
        case Five :
	   Bn = 0.0;
	break;  
   }


   double result = 0. ;
   switch( EXPR )
   {
      case ONE :
        result = ys ;
        break ;
      case TWO :
	case FOUR:
	{
        double const ycoord = arg(1)->to_double(ct) ;
        double const Uscale = arg(2)->to_double(ct) ;
	  double mid = (Right+Left)/2.;

	  ys = ys*wall;

        double tmp1;
	  if (Bn == 0.)
	  {
	     tmp1 = 1.;
           ys = 0.;
	  }
	  else
	  {
           tmp1 = Bn/(2.*ys);
	  }

        double tmp2 = (wall-ys)*(wall-ys); // mid + wall = Right
        double tmp3 = tmp1 * tmp2 ;
	  double tmp4;
	  if(ycoord >=mid)
             tmp4 = tmp1 * (tmp2 - pow((ycoord-mid-ys),2) );
        else
    	       tmp4 = tmp1 * (tmp2 - pow((mid-ycoord-ys),2) );

//	        double integral = tmp1*tmp2*ys+tmp1*tmp2*(wall-ys)-tmp1/3.*pow((wall-ys),3);


	  if( ycoord >= Left && ycoord <=Right) {
	     if( ycoord >= (mid-ys)  && ycoord <= (mid+ys) )
           {
             result = tmp3 ;
           } else {
             result = tmp4 ;
           }
	     result=result*Uscale;//integral;
	  } else {
	     result = 0.0;
	  }
	}

        break ;
      case THREE:
	{
	   ys = 0.;
	   double const ycoord = arg(0)->to_double(ct) ;
	   double const Uscale = arg(1)->to_double(ct) ;
	   Left  = arg(2)->to_double(ct) ;
	   Right = arg(3)->to_double(ct) ;

	   double diameter = Right-Left;
	   if( ycoord >= Left && ycoord <=Right)
	   {
	      result = (4.*(ycoord-Left) * (diameter -ycoord + Left))/(diameter*diameter) ;
		result = Uscale*result ;
	   }
	   else	result = 0.0;
	}
	   break ;
      case Five:
	{
	   ys = 0.;
	   double const ycoord = arg(0)->to_double(ct) ;
	   double const Uscale = arg(1)->to_double(ct) ;
	   Left  = arg(2)->to_double(ct) ;
	   Right = arg(3)->to_double(ct) ;
           double T = arg(4)->to_double(ct);
           double Tinterval =arg(5)->to_double(ct);
           double FT;
	   double diameter = Right-Left;
           if ( T <= (Tinterval) )
            {
            FT= (1.0/Tinterval)*T;
              }
            else FT=1.0;
          
	   if( ycoord >= Left && ycoord <=Right) 
	   {
	      result = (4.*(ycoord-Left) * (diameter -ycoord + Left))/(diameter*diameter) ;
		result = Uscale*result*FT ;
	   }
	   else	result = 0.0;
	}
	   break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const&
RS_Bingham:: usage( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Bingham:: usage" ) ;

   static std::string result ;

   switch( EXPR )
   {
      case ONE :
         result = "Bingham_YieldSurface($DS_Bn, $DS_Q)" ;
         break ;
      case TWO :
	 result = "Bingham_VelocityProfile($DS_Bn, $DS_Y, $DS_Uscale, $DS_Left, $DS_Right, $DS_Q)" ;
         break ;
      case THREE :
         result = "Newtonian_VelocityProfile($DS_Y, $DS_Uscale, $DS_Left, $DS_Right)" ;
	 break ;
      case FOUR :
	 result = "Bingham_PipeVelocityProfile($DS_Bn, $DS_Y, $DS_Uscale, $DS_Left, $DS_Right, $DS_Q)" ;
	 break ;
      case Five :
         result = "Transient_Newtonian_VelocityProfile($DS_Y, $DS_Uscale, $DS_Left, $DS_Right,$DS_T,$DS_Tinterval)" ;
	 break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
RS_Bingham:: valid_arguments(
                              PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "RS_Bingham:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = false ;
   switch( EXPR )
   {
      case ONE :
         result =  ( some_arguments->count() == 2 ) &&
			 ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::Double )&&
			 ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double );
         break ;
      case TWO :
         result = ( some_arguments->count() == 6 ) &&
           ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::Double )&&
           ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double )&&
           ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double )&&
           ( extract_arg( some_arguments, 3 )->data_type() == PEL_Data::Double )&&
           ( extract_arg( some_arguments, 4 )->data_type() == PEL_Data::Double )&&
           ( extract_arg( some_arguments, 5 )->data_type() == PEL_Data::Double ) ;
         break ;
	case THREE :
	   result = ( some_arguments->count() == 4 ) &&
		   ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 3 )->data_type() == PEL_Data::Double ) ;
	   break ;
	case FOUR :
	   result = ( some_arguments->count() == 6 ) &&
		   ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 3 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 4 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 5 )->data_type() == PEL_Data::Double ) ;
	   break ;
      case Five :
	   result = ( some_arguments->count() == 6 ) &&
		   ( extract_arg( some_arguments, 0 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 1 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 2 )->data_type() == PEL_Data::Double )&&
                   ( extract_arg( some_arguments, 3 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 4 )->data_type() == PEL_Data::Double )&&
		   ( extract_arg( some_arguments, 5 )->data_type() == PEL_Data::Double ) ;
	   break ;
   }

   return( result ) ;
}
