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

#include <PDE_2D_Q0isoQ1_4nodes.hh>

#include <PEL.hh>

#include <GE_Point.hh>
#include <GE_ReferenceSquare.hh>

using std::string ;

//----------------------------------------------------------------------
PDE_2D_Q0isoQ1_4nodes* PDE_2D_Q0isoQ1_4nodes::uniqueInstance = new PDE_2D_Q0isoQ1_4nodes() ;
//----------------------------------------------------------------------



//----------------------------------------------------------------------
PDE_2D_Q0isoQ1_4nodes:: ~PDE_2D_Q0isoQ1_4nodes( void )
//----------------------------------------------------------------------
{
   uniqueInstance = 0 ;
}


//----------------------------------------------------------------------
PDE_2D_Q0isoQ1_4nodes:: PDE_2D_Q0isoQ1_4nodes( void )
//----------------------------------------------------------------------
   : PDE_ReferenceElement( "PDE_2D_Q0isoQ1_4nodes", GE_ReferenceSquare::object() )
{
   append_node( GE_Point::create( this, 0.25, 0.25 ) ) ;
   append_node( GE_Point::create( this, 0.75, 0.25 ) ) ;
   append_node( GE_Point::create( this, 0.75, 0.75 ) ) ;
   append_node( GE_Point::create( this, 0.25, 0.75 ) ) ;
}



//----------------------------------------------------------------------
double
PDE_2D_Q0isoQ1_4nodes:: N_local( size_t node,
                      GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_2D_Q0isoQ1_4nodes:: N_local" ) ;
   PEL_CHECK_PRE( N_local_PRE( node, pt_ref ) ) ;

   double bf = 0.0 ;

   double x = pt_ref->coordinate( 0 ) ;
   double y = pt_ref->coordinate( 1 ) ;

   switch( node )
   {
    case 0 :
               if( x <= 0.5 && y <= 0.5) bf = 1.0;
               break ;
    case 1 :
               if( x >= 0.5 && y <= 0.5) bf = 1.0;
               break ;
    case 2 :
               if( x >= 0.5 && y >= 0.5) bf = 1.0;
               break ;
    case 3 :
               if( x <= 0.5 && y >= 0.5) bf = 1.0;
               break ;
   }

   return( bf ) ;
}



//----------------------------------------------------------------------
double
PDE_2D_Q0isoQ1_4nodes:: dN_local( size_t node,
                       size_t a,
                       GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_2D_Q0isoQ1_4nodes:: dN_local" ) ;
   PEL_CHECK_PRE( dN_local_PRE( node, a, pt_ref ) ) ;

   double dbf = 0.0 ;

   return( dbf );
}



//----------------------------------------------------------------------
double
PDE_2D_Q0isoQ1_4nodes:: d2N_local( size_t node,
                        size_t a, size_t b,
                        GE_Point const* pt_ref ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_2D_Q0isoQ1_4nodes:: d2N_local" ) ;
   PEL_CHECK_PRE( d2N_local_PRE( node, a, b, pt_ref ) ) ;

   double d2bf = 0.0 ;

   return( d2bf ) ;
}
