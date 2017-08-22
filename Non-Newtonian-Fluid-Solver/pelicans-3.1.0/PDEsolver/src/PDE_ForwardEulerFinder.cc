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

#include <PDE_ForwardEulerFinder.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_VectorIterator.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_Point.hh>
#include <GE_Vector.hh>

#include <LA_GaussLU_DS.hh>
#include <LA_DenseMatrix.hh>
#include <LA_SeqVector.hh>

#include <PDE_BoundFE.hh>
#include <PDE_CellFE.hh>
#include <PDE_DiscreteField.hh>
#include <PDE_GridFE.hh>
#include <PDE_FaceFE.hh>

//---------------------------------------------------------------------------
PDE_ForwardEulerFinder*
PDE_ForwardEulerFinder:: create( PDE_LocalFEcell* fem,
                                 PEL_ModuleExplorer const* exp,
                                 PDE_GridFE const* grid )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ForwardEulerFinder:: create" ) ;
   PEL_CHECK_PRE( fem!=0 ) ;
   PEL_CHECK_PRE( exp!=0 ) ;
   PEL_CHECK_PRE( grid!=0 ) ;
   PEL_CHECK_PRE( fem->nb_space_dimensions()==grid->nb_space_dimensions() ) ;

   PDE_ForwardEulerFinder* result =
                         new PDE_ForwardEulerFinder( fem, exp, grid ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == fem ) ;
//   PEL_CHECK_POST( result->localFEcell() == fem ) ;
   PEL_CHECK_POST( !result->search_has_been_performed() ) ;
   return( result ) ;   
}

//---------------------------------------------------------------------------
PDE_ForwardEulerFinder:: PDE_ForwardEulerFinder(
                                       PDE_LocalFEcell* fem,
                                       PEL_ModuleExplorer const* exp,
                                       PDE_GridFE const* grid )
//---------------------------------------------------------------------------
   : PDE_CFootFinder( fem, exp ),
     MESHES( 0 ),
     ACCEPTABLE( 0 ),
     FOOT( 0 ),
     P( 0 ),
     Q( 0 ),
     tr_P( 0 ),
     tr_Q( 0 ),
     P_refC( 0 ),
     Q_refC( 0 ),
     P_refB( 0 ),
     Q_refB( 0 ),
     PM( 0 ),
     QM( 0 ),
     aP( 0 ),
     PT_REF( 0 )   
{
   PEL_LABEL( "PDE_ForwardEulerFinder:: PDE_ForwardEulerFinder(fem,exp,grid)" ) ;
   
   DIM = grid->nb_space_dimensions() ;

   P =  GE_Point::create( this, DIM ) ;
   Q = GE_Point::create( this, DIM ) ;
   tr_P = GE_Point::create( this, DIM ) ;
   tr_Q = GE_Point::create( this, DIM ) ;
   Q_refC = GE_Point::create( this, DIM ) ;
   P_refC = GE_Point::create( this, DIM ) ;
   Q_refB = GE_Point::create( this, DIM-1 ) ;
   P_refB = GE_Point::create( this, DIM-1 ) ;
   PM = GE_Vector::create( this, DIM ) ;
   QM = GE_Vector::create( this, DIM ) ;
   aP = GE_Vector::create( this, DIM ) ;

   ACCEPTABLE.re_initialize( grid->nb_cells() ) ;
   MESHES = PEL_List::create( this ) ;

   FOOT = GE_Point::create( this, DIM ) ;
   FOOT_REF = GE_Point::create( this, DIM ) ;

   CELL_FT = GE_Point::create( this, DIM ) ;
   CELL_FT_REF = GE_Point::create( this, DIM ) ;

   BOUND_FT = GE_Point::create( this, DIM ) ;
   BOUND_FT_REF = GE_Point::create( this, DIM-1 ) ;
   PT_REF = GE_Point::create( this, DIM ) ;

   CELL_ITER_MAX = exp->int_data( "nb_iter_max_when_searching_in_cell" ) ;

   if( CELL_ITER_MAX < 2 )
   {
      PEL_Error::object()->raise_bad_data_value( 
                                  exp,
                                  "nb_iter_max_when_searching_in_cell",
                                  "greater or equal to 2" ) ;
   }

   BOUND_ITER_MAX = exp->int_data( "nb_iter_max_when_searching_in_bound" ) ;
   if( BOUND_ITER_MAX < 2 )
   {
      PEL_Error::object()->raise_bad_data_value( 
                                  exp,
                                  "nb_iter_max_when_searching_in_bound",
                                  "greater or equal to 2" ) ;
   }

   REF_DIST_MAX = exp->double_data( "distance_max_characteristic_foot" ) ;

   SOLVER = LA_GaussLU_DS::create( this, 1.E-30 ) ;
   A = LA_DenseMatrix::create( this, DIM,  DIM ) ;
   b  = LA_SeqVector::create( this, DIM ) ;
   dX = LA_SeqVector::create( this, DIM ) ;

   PEL_CHECK_INV( invariant() ) ;
}

//---------------------------------------------------------------------------
PDE_ForwardEulerFinder:: ~PDE_ForwardEulerFinder( void )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ForwardEulerFinder:: ~PDE_ForwardEulerFinder" ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//---------------------------------------------------------------------------
void
PDE_ForwardEulerFinder:: search_foot( 
                                              PDE_DiscreteField const* aa,
                                              size_t level,
                                              double time_step,
                                              GE_Point const* head,
                                              PDE_CellFE const* head_cell )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ForwardEulerFinder:: search_foot" ) ;
   PEL_CHECK( search_foot_PRE( aa, level, time_step,
                                              head, head_cell ) ) ;
   PEL_CHECK_INV( invariant() ) ;
   
   initialize_trial_cells( head_cell ) ;

   FOUND = false ;
   do
   {
      PDE_CellFE const* trial_cell = pop_one_trial_cell() ;

      search_in_cell( aa, level, time_step, head, trial_cell ) ;

      if( FOUND_IN_CELL )
      {
         FOUND = true ;
         set_found_cell( trial_cell ) ;
         FOOT->set( CELL_FT ) ; 
         FOOT_REF->copy( CELL_FT_REF ) ; 
      }
      else
      {
         PEL_VectorIterator* it = 
             PEL_VectorIterator::create( 0, trial_cell->faces() ) ;
         for( it->start() ; it->is_valid() && !FOUND ; it->go_next() )
         {
            PDE_FaceFE const* side = static_cast<PDE_FaceFE*>( it->item() ) ;
               
            if( side->has_adjacent_bound() )
            {
               PDE_BoundFE const* trial_bound = side->adjacent_bound() ;

               search_in_bound( aa, level, time_step, head, trial_bound ) ;

               if( FOUND_IN_BOUND )
               {
                  FOUND = true ;
                  set_found_bound( trial_bound ) ;
                  FOOT->set( BOUND_FT ) ;
                  FOOT_REF->copy( BOUND_FT_REF ) ;
               }
            }
            else
            {
               if( side->is_periodic() ) raise_periodicity_error() ;
               extend_trial_cells( side->adjacent_cell_other_than( trial_cell ), 
                                   head_cell ) ;
            }
         }

         it->destroy() ; it = 0;
      }

   } while( !FOUND && more_trial_cells() > 0 ) ;

   PEL_CHECK_INV( invariant() ) ;
   PEL_CHECK( search_foot_POST( aa, level, time_step,
                                               head, head_cell ) ) ;
}

//--------------------------------------------------------------------------
GE_Point const*
PDE_ForwardEulerFinder:: foot( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ForwardEulerFinder:: foot" ) ;
   PEL_CHECK_PRE( foot_PRE() ) ;

   GE_Point const* result = FOOT ;

   PEL_CHECK_POST( foot_POST( result ) ) ;
   return( result ) ;
}

//--------------------------------------------------------------------------
GE_Point const*
PDE_ForwardEulerFinder:: foot_ref( void ) const
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ForwardEulerFinder:: foot_ref" ) ;
   PEL_CHECK_PRE( foot_ref_PRE() ) ;

   GE_Point const* result = FOOT_REF ;

   PEL_CHECK_POST( foot_ref_POST( result ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PDE_ForwardEulerFinder:: search_in_cell( 
                                   PDE_DiscreteField const* aa,
                                   size_t level, 
                                   double time_step,
                                   GE_Point const* head,
                                   PDE_CellFE const* trial_cell )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ForwardEulerFinder:: search_in_cell" ) ;
   
   FOUND_IN_CELL = false ;
   CELL_FT->set( GE_Point::origin(DIM) ) ;
   CELL_FT_REF->set( GE_Point::origin(DIM) ) ;
  
   double const dX_tol = 1.E-6 ;
   double h ;

   GE_Mpolyhedron const* poly = trial_cell->polyhedron() ;
   GE_ReferencePolyhedron const* ref_poly =  poly->reference_polyhedron() ;

   P_refC->set( ref_poly->center() ) ;

   size_t iter = 0 ;
   bool converged = false ;
   do
   {
      iter++ ;

      poly->apply_mapping( P_refC, P ) ;

      transport( aa, level, time_step, trial_cell, P_refC, P, tr_P ) ;

      for( size_t d=0 ; d<DIM ; d++ )
      {
         ref_poly->build_neighbor( P_refC, d, Q_refC, h ) ;

         poly->apply_mapping( Q_refC, Q ) ;
         transport( aa, level, time_step, trial_cell, Q_refC, Q, tr_Q ) ;
         for( size_t d2=0 ; d2<DIM ; d2++ )
         {
            A->set_item( d2, d, ( tr_Q->coordinate(d2) -
                                  tr_P->coordinate(d2) ) / h ) ;
         }
         if( PEL::abs(A->item(d,d))<1.E-12 ) 
         {
            A->set_item( d, d, 1.E-12 ) ;
         }
         b->set_item(d, -( tr_P->coordinate(d) - head->coordinate(d) ) ) ;
      }
      
      A->synchronize() ;
      b->synchronize() ;
      SOLVER->set_matrix( A ) ;
      SOLVER->solve( b, dX ) ;
      SOLVER->unset_matrix() ;

      converged = ( dX->max_norm() < dX_tol ) ;

      for( size_t d=0 ; d < DIM ; ++d )
      {
         P_refC->set_coordinate( d, P_refC->coordinate(d) + dX->item(d) ) ;
      }

      if( !ref_poly->contains( P_refC ) )
      {
         ref_poly->project( P_refC ) ;
         if( converged ) FOUND_IN_CELL = false ;
      }
      else
      {
         if( converged ) FOUND_IN_CELL = true ;
      }

   } while( !converged && iter<CELL_ITER_MAX ) ;

   CELL_FT_REF->set( P_refC ) ;
   poly->apply_mapping( CELL_FT_REF, CELL_FT ) ;
}

//----------------------------------------------------------------------
void
PDE_ForwardEulerFinder:: search_in_bound( 
                                   PDE_DiscreteField const* aa,
                                   size_t level, 
                                   double time_step,
                                   GE_Point const* head,
                                   PDE_BoundFE const* trial_bound )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ForwardEulerFinder:: search_in_bound" ) ;
   
   FOUND_IN_BOUND = false ;
   BOUND_FT->set( GE_Point::origin( DIM ) ) ;
   BOUND_FT_REF->set( GE_Point::origin( DIM-1 ) ) ;

   if( trial_bound->isInFlow( level, aa ) )
   {
      PEL_ASSERT( DIM == 2 ) ;
      PEL_ASSERT( trial_bound->polyhedron()->dimension() == 1 ) ;

      double const ds_tol = 1.E-6 ;
      double f, fe, ds, h ;

      GE_Mpolyhedron const* poly = trial_bound->polyhedron() ;
      GE_ReferencePolyhedron const* ref_poly = poly->reference_polyhedron() ;

      double s = ref_poly->center()->coordinate( 0 ) ;
      P_refB->set_coordinate( 0, s ) ;

      size_t iter = 0 ;
      bool converged = false ;
      do
      {
         iter++ ;

         ref_poly->build_neighbor( P_refB, 0, Q_refB, h ) ;

         poly->apply_mapping( P_refB, P ) ;
         poly->apply_mapping( Q_refB, Q ) ;

         PM->re_initialize( head, P ) ;
         QM->re_initialize( head, Q ) ;

         f  = PM->component(0)  * trial_bound->value( aa, level, P_refB, 1 ) -
              PM->component(1)  * trial_bound->value( aa, level, P_refB, 0 ) ;

         fe = QM->component(0) * trial_bound->value( aa, level, Q_refB, 1 ) -
              QM->component(1) * trial_bound->value( aa, level, Q_refB, 0 ) ;

         ds = - f * h/(fe-f) ;

         s += ds ;
         P_refB->set_coordinate( 0, s ) ;

         converged = PEL::abs(ds) < ds_tol ;
 
         if( !ref_poly->contains( P_refB ) )
         {
            ref_poly->project( P_refB ) ;
            s = P_refB->coordinate( 0 ) ;
         }
         else if( converged )
         {
            poly->apply_mapping( P_refB, P ) ;     
            PM->re_initialize( head, P ) ;
            for( size_t d=0 ; d<DIM ; ++d )
            {
               aP->set_component( d, trial_bound->value( aa, level, P_refB, d ) ) ;
            }
            if( aP->dot_product( PM )>0.0  &&  PM->norm()/aP->norm()<=time_step )
            {
               FOUND_IN_BOUND = true ;
               BOUND_FT->set( P ) ;
               BOUND_FT_REF->set( P_refB ) ;
            }
         }
      } while( !converged && iter<BOUND_ITER_MAX ) ;
   }
}

//----------------------------------------------------------------------
void
PDE_ForwardEulerFinder:: transport( PDE_DiscreteField const* aa,
                                    size_t level,
                                    double time_step,
                                    PDE_CellFE const* cell_of_X, 
                                    GE_Point const* X_ref,
                                    GE_Point const* X,
                                    GE_Point* tr_X ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ForwardEulerFinder:: transport" ) ;
   
   for( size_t d=0 ; d<tr_X->nb_coordinates() ; d++ )
   {
      double v = cell_of_X->value( aa, level, X_ref, d ) ;
      tr_X->set_coordinate( d, X->coordinate(d)+time_step*v ) ;
   }
}

//------------------------------------------------------------------------
void
PDE_ForwardEulerFinder:: initialize_trial_cells( PDE_CellFE const* head_cell )
//------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ForwardEulerFinder:: initialize_trial_cells" ) ;
   
   MESHES->clear() ;
   ACCEPTABLE.set( true ) ;

   MESHES->append( const_cast<PDE_CellFE*>( head_cell  ) ) ;
   ACCEPTABLE( head_cell->id_number() ) = false ;
}

//------------------------------------------------------------------------
PDE_CellFE const*
PDE_ForwardEulerFinder:: pop_one_trial_cell( void )
//------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ForwardEulerFinder:: pop_one_trial_cell" ) ;
   
   PDE_CellFE const* result = static_cast<PDE_CellFE const*>( MESHES->at(0) ) ;
   MESHES->remove_at( 0 ) ;

   return( result ) ;
}

//--------------------------------------------------------------------------
void
PDE_ForwardEulerFinder:: extend_trial_cells( PDE_CellFE const* a_cell,
                                             PDE_CellFE const* head_cell )
//--------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ForwardEulerFinder:: extend_trial_cells" ) ;
   
   if( ACCEPTABLE( a_cell->id_number() ) )
   {
      head_cell->polyhedron()->apply_inverse_mapping( 
                                      a_cell->polyhedron()->center(),
                                      PT_REF ) ;
      GE_Point const* ref_center =
                head_cell->polyhedron()->reference_polyhedron()->center() ;
      bool ok = ( PT_REF->distance(ref_center) <= REF_DIST_MAX ) ;

      if( ok )
      {
         if( a_cell->polyhedron()->contains( CELL_FT ) )
         {
            MESHES->prepend( const_cast<PDE_CellFE*>(a_cell) ) ;
         }
         else
         {
            MESHES->append( const_cast<PDE_CellFE*>(a_cell) ) ;
         }
      }

      ACCEPTABLE( a_cell->id_number() ) = false ;
   }
}

//-------------------------------------------------------------------------
bool
PDE_ForwardEulerFinder:: more_trial_cells( void ) const
//-------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_ForwardEulerFinder:: more_trial_cells" ) ;
   
   return( MESHES->count() > 0 ) ;
}

//-------------------------------------------------------------------------
bool
PDE_ForwardEulerFinder:: invariant( void ) const
//-------------------------------------------------------------------------
{
   PEL_ASSERT( PDE_CFootFinder::invariant() ) ;
   return( true ) ;
}
