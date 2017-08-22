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

#include <EXT_SPARSKITsolver.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Bool.hh>
#include <PEL_Error.hh>
#include <PEL_Exec.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Variable.hh>

#include <LA_MatrixIterator.hh>
#include <LA_Matrix.hh>
#include <LA_SeqVector.hh>
#include <iostream>

using std::cout ;
using std::endl ;
using std::string ;

extern "C" 
{
   void lusol_(int * nrow,double * my_rhs,double * xran,
               double * au,int * jau,int * ju) ;
   void ilu0_(int * nrow, double * a, int * ja, int * ia,
              double * au, int * jau, int * ju, int * iw, int * ierr) ;
   void milu0_(int * nrow, double * a, int * ja, int * ia,
               double * au, int * jau, int * ju, int * iw, int * ierr) ;
   
   void runrc_(int * nrow,double * my_rhs,double * x,int * ipar,double * fpar,
               double * vv,double * xran,double * a,int * ja,int * ia,
               double * au,int * jau,int * ju, void (*)() ) ;
   
   void gmres_( void ) ;
   void fgmres_( void ) ;
}
static size_t const FORTRAN = 1 ;

EXT_SPARSKITsolver const* EXT_SPARSKITsolver::PROTOTYPE = 
                                                    new EXT_SPARSKITsolver() ;

//----------------------------------------------------------------------------
EXT_SPARSKITsolver:: EXT_SPARSKITsolver( void ) 
//----------------------------------------------------------------------------
   : LA_Solver( "EXT_SPARSKITsolver" )
   , EXP( 0 )
{   
   PEL_Bool* val = PEL_Bool::create( 0, true ) ;
   PEL_Exec::add_variable_to_execution_context(
                         PEL_Variable::object( "BS_with_SPARSKIT" ), val ) ;
}

//----------------------------------------------------------------------------
EXT_SPARSKITsolver:: ~EXT_SPARSKITsolver( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_SPARSKITsolver:: ~EXT_SPARSKITsolver" ) ;
   
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
   if( matrix_is_set() )
   {
      unset_matrix() ;
   }
}

//----------------------------------------------------------------------------
EXT_SPARSKITsolver*
EXT_SPARSKITsolver:: create( PEL_Object* a_owner,
                             PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_SPARSKITsolver:: create" ) ;

   EXT_SPARSKITsolver* result = new EXT_SPARSKITsolver( a_owner, exp ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( !result->matrix_is_set() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
EXT_SPARSKITsolver*
EXT_SPARSKITsolver:: create_replica( PEL_Object* a_owner, 
                                    PEL_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_SPARSKITsolver:: create_replica" ) ;
   PEL_CHECK_PRE( create_replica_PRE(  a_owner, exp ) ) ;
   
   EXT_SPARSKITsolver* result = new EXT_SPARSKITsolver( a_owner, exp ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
EXT_SPARSKITsolver:: EXT_SPARSKITsolver( PEL_Object* a_owner,
                                         PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
   : LA_Solver( a_owner )
   , EXP( exp->create_clone( this ) )
   , a( 0 )
   , my_solver( 0 )
{
   PEL_LABEL( "EXT_SPARSKITsolver:: EXT_SPARSKITsolver" ) ;
   raise_fatal_error_if_not_sequential() ;
   initialize() ;
}

//----------------------------------------------------------------------------
EXT_SPARSKITsolver*
EXT_SPARSKITsolver:: create_clone( PEL_Object* a_owner ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_SPARSKITsolver:: create_clone" ) ;
   
   EXT_SPARSKITsolver* result = new EXT_SPARSKITsolver( a_owner, this ) ;

   PEL_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
EXT_SPARSKITsolver:: EXT_SPARSKITsolver( PEL_Object* a_owner,
                                         EXT_SPARSKITsolver const* other ) 
//----------------------------------------------------------------------------
   : LA_Solver( a_owner, other )
   , EXP( other->EXP->create_clone( this ) )
   , a( 0 )
   , my_solver( 0 )
{
   PEL_LABEL( "EXT_SPARSKITsolver:: EXT_SPARSKITsolver" ) ;
   raise_fatal_error_if_not_sequential() ;
   initialize() ;
}

//----------------------------------------------------------------------------
void
EXT_SPARSKITsolver:: initialize( void ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_SPARSKITsolver:: initialize" ) ;

   PEL_ModuleExplorer* sexp =
      EXP->create_subexplorer( 0, "SPARSKIT_IterativeSolver" ) ;
   build_ksp( sexp ) ;
   sexp->destroy() ; sexp = 0 ;
   
   sexp = EXP->create_subexplorer( 0, "SPARSKIT_Preconditioner" ) ;
   build_pc( sexp ) ;
   sexp->destroy() ; sexp = 0 ;
   
}

//----------------------------------------------------------------------------
void
EXT_SPARSKITsolver:: build_ksp( PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_SPARSKITsolver:: build_ksp" ) ;
   PEL_CHECK( exp != 0 ) ;
   
   algo = exp->string_data( "type" ) ;
   double toler = exp->double_data( "tolerance" ) ;
   int maxit = exp->int_data( "nb_iterations_max" ) ;
   int restart = exp->int_data( "restart" ) ;
   ipar[1-FORTRAN] = 0 ;
   ipar[2-FORTRAN] = 1 ;
   ipar[3-FORTRAN] = 2 ;
   ipar[5-FORTRAN] = restart ;
   ipar[6-FORTRAN] = maxit ;
   fpar[1-FORTRAN] = 0.0 ;
   fpar[2-FORTRAN] = toler ;
   set_iterative( true );
   if( algo=="GMRES" )
   {
      my_solver = gmres_ ;
   }
   else if( algo=="FGMRES" )
   {
      my_solver = fgmres_ ;
   }
   else
   {
      PEL_Error::object()->raise_plain(
         "Allowed Krylov methods are GMRES and FGMRES " ) ;
   }
}

//----------------------------------------------------------------------------
void
EXT_SPARSKITsolver:: build_pc( PEL_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_SPARSKITsolver:: build_pc" ) ;
   PEL_CHECK( exp != 0 ) ;
   
   precond = exp->string_data( "type" ) ;
   if( precond=="none" )
   {
   }
   else if( precond=="ilu0" )
   {
   }
   else if( precond=="milu0" )
   {
   }
   else
   {
      PEL_Error::object()->raise_plain(
         "Unknown (not in : none,ilu0,milu0) precond. "+precond ) ;
   }

}

//----------------------------------------------------------------------------
void
EXT_SPARSKITsolver:: unset_matrix_self( void )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_SPARSKITsolver:: unset_matrix_self" ) ;
   PEL_CHECK( unset_matrix_self_PRE() ) ;
   
   if( a!=0 )
   {
      delete [] a ; a=0 ;
      delete [] ia ; ia=0 ;
      delete [] ja ; ja=0 ;
      delete [] au ; au=0 ;
      delete [] jau ; jau=0 ;
      delete [] ju ; ju=0 ;
      delete [] vv ; vv=0 ;
   }
   
}

//----------------------------------------------------------------------------
void
EXT_SPARSKITsolver:: set_matrix_self(
                           LA_Matrix const* mat, bool &ok, bool same_pattern )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_SPARSKITsolver:: initialize_self" ) ;
   PEL_CHECK( set_matrix_self_PRE( mat, same_pattern ) ) ;

   size_t nnz = mat->nb_stored_items() ;
   nrow = (int) mat->nb_rows() ;
   a = new double [nnz] ;
   ja = new int [nnz] ;
   ia = new int [nrow+1] ;
   au = new double [nnz+10*nrow] ;
   jau = new int [nnz+10*nrow] ;
   ju = new int [nrow+1] ;
   
   LA_MatrixIterator* it = mat->create_stored_item_iterator( this ) ;
   size_t el = 0 ;
   for( int i=0 ; i<nrow ; i++ )
   {
      ia[i] = el + FORTRAN ;
      for( it->start_row_items(i) ;
           it->is_valid() ;
           it->go_next() )
      {
         size_t j = it->col() ;
         double Aij = it->item() ;
         a[el] = Aij ;
         ja[el] = j + FORTRAN ;
         el++ ;
      }
   }
   ia[nrow]=el+ FORTRAN ;
   
   PEL_ASSERT( el <= nnz ) ;
   

   int n = nrow ;
   int m = ipar[5-FORTRAN] ;
   int lwk = PEL::max( (n+3)*(m+2) + (m+1)*m/2 , 2*n*(m+1) + (m+1)*m/2 + 3*m + 2 );
   ipar[4-FORTRAN] = lwk ;
   vv = new double [ lwk ] ;

   build( ok ) ;
   
   PEL_CHECK_POST( set_matrix_self_POST( mat, ok ) ) ;
}

//----------------------------------------------------------------------------
void
EXT_SPARSKITsolver:: build( bool& success )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_SPARSKITsolver:: build" ) ;
   PEL_CHECK( au!=0 ) ;
   int * iw = new int [ 3*nrow ] ;
   int ierr = 0 ;
   if( precond=="none" )
   {
      int j0 = nrow + 2 ;
      for( int i=0 ; i<nrow ; i++ )
      {
         au[i] = 1.0 ;
         ju[i] = j0 ;
         jau[i] = j0 ;
      }
      jau[nrow] = j0 ;
      ju[nrow] = j0 ;
   }
   else if( precond=="ilu0" )
   {
      ilu0_(&nrow,a, ja, ia, au, jau, ju, iw, &ierr) ;
   }
   else if( precond=="milu0" )
   {
      milu0_(&nrow,a, ja, ia, au, jau, ju, iw, &ierr) ;
   }
   else
   {
      PEL_Error::object()->raise_plain( "Unknown "+precond ) ;
   }
   
   delete [] iw ;
   
   success = ierr==0 ;

}

//----------------------------------------------------------------------------
void
EXT_SPARSKITsolver:: solve_self( LA_Vector const* b, LA_Vector* x,
                                 size_t &nb_iter, bool &ok )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_SPARSKITsolver:: solve_self" ) ;
   PEL_CHECK( solve_self_PRE( b, x ) ) ;
   
   PEL_CHECK( au!=0 ) ;
   
   PEL_CHECK( dynamic_cast<LA_SeqVector const* >( b ) != 0 ) ;
   LA_SeqVector const* brhs = static_cast<LA_SeqVector const* >( b ) ;
   
   PEL_CHECK( dynamic_cast<LA_SeqVector* >( x ) != 0 ) ;
   LA_SeqVector const* bx = static_cast<LA_SeqVector const* >( x ) ;
   
   double* my_rhs = new double [nrow] ;
   for( int i=0 ; i<nrow ; i++ )
   {
      my_rhs[i] = brhs->item(i) ;
   }
   double * my_x = new double [ nrow ] ;
   if( zero_initial_guess() )
   {
      lusol_(&nrow,my_rhs,my_x,au,jau,ju) ;
   }
   else
   {
      for( int i=0 ; i<nrow ; i++ ) my_x[i] = bx->item(i) ;
   }
   
   runrc_(&nrow,my_rhs,my_x,ipar,fpar,vv,my_x,a,ja,ia,au,jau,ju, my_solver ) ;
   
   for( int i=0 ; i<nrow ; i++ )
   {
      x->set_item( i, my_x[i] ) ;
   }
   delete [] my_x ;
   delete [] my_rhs ;
   
   ok = ipar[1-FORTRAN]==0 ;
   nb_iter = ipar[7-FORTRAN] ;
   
   PEL_CHECK_POST( solve_self_POST( b, x, nb_iter, ok ) ) ;
}

//----------------------------------------------------------------------------
void
EXT_SPARSKITsolver:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "EXT_SPARSKITsolver:: print" ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   std::string s( indent_width, ' ') ;
   os << s << "linear solver : \"SPARSKIT_" << algo << "\"" << std::endl ;
   os << s << "linear solver error bound : "
      << fpar[2-FORTRAN] << std::endl ;
   os << s << "linear solver maximal number of iterations : "
      << ipar[6-FORTRAN] << std::endl ;
   os << s << "linear solver restart iteration : "
      << ipar[5-FORTRAN] << std::endl ;   
   
   // Preconditioner :
   os << s << "preconditioner : \"SPARSKIT_" << precond << "\"" << std::endl ;
}
