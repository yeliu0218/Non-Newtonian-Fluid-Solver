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

#include <PDE_AdaptationRequestFromIndicator.hh>

#include <PEL.hh>
#include <PEL_List.hh>
#include <PEL_ListIterator.hh>
#include <PEL_ModuleExplorer.hh>

#include <GE_Mpolyhedron.hh>

#include <PDE_AdaptationIndicator.hh>
#include <PDE_BasisFunctionCell.hh>
#include <PDE_CellFE.hh>
#include <PDE_ReferenceElement.hh>

#include <iostream>

using std::cout ;
using std::endl ;

//---------------------------------------------------------------------------
PDE_AdaptationRequestFromIndicator*
PDE_AdaptationRequestFromIndicator:: create( 
                                     PEL_Object* a_owner, 
                                     PEL_List const* indicators,
                                     size_t highest_refinement_level,
                                     size_t verbose_level  )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequestFromIndicator:: create" ) ;
   
   PDE_AdaptationRequestFromIndicator* result = 
     new PDE_AdaptationRequestFromIndicator( a_owner, 
                                             indicators, 
                                             highest_refinement_level,
                                             verbose_level ) ;
   
   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_AdaptationRequestFromIndicator:: PDE_AdaptationRequestFromIndicator( 
                                         PEL_Object* a_owner,
                                         PEL_List const* indicators,
                                         size_t highest_refinement_level,
                                         size_t verbose_level )
//---------------------------------------------------------------------------
   : PDE_AdaptationRequest( a_owner, highest_refinement_level, verbose_level )
   , INDICS_IT( indicators->create_iterator( this ) )
{
}

//---------------------------------------------------------------------------
PDE_AdaptationRequestFromIndicator:: ~PDE_AdaptationRequestFromIndicator( void  )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
bool
PDE_AdaptationRequestFromIndicator:: to_be_refined( 
                                       PDE_CellFE const* cell,
                                       PDE_BasisFunctionCell const* bf,
                                       PDE_ReferenceElement const* elm,
                                       size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequestFromIndicator:: to_be_refined" ) ;
   PEL_CHECK_PRE( to_be_refined_PRE( cell, bf, elm, local_node ) ) ;

   bool result = false ;
   INDICS_IT->start() ;
   for( ; INDICS_IT->is_valid() ; INDICS_IT->go_next() )
   {
      PDE_AdaptationIndicator const* indic = 
          static_cast< PDE_AdaptationIndicator* >( INDICS_IT->item() ) ;
          
      //??? on pourrait faire tous les calculs ensemble
      double bf_indicator = basis_function_indicator( indic, bf, cell ) ;

      result = result || indic->to_be_refined( bf_indicator, 
                                               cell->polyhedron(), 
                                               elm, 
                                               local_node ) ;
   }
   
   return( result ) ;
}

//---------------------------------------------------------------------------
bool
PDE_AdaptationRequestFromIndicator:: to_be_unrefined(  
                                         PDE_CellFE const* cell,
                                         PDE_BasisFunctionCell const* bf,
                                         PDE_ReferenceElement const* elm,
                                         size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequestFromIndicator:: to_be_unrefined" ) ;
   PEL_CHECK_PRE( to_be_unrefined_PRE( cell, bf, elm, local_node ) ) ;

   bool result = true ;
   INDICS_IT->start() ;
   for( ; INDICS_IT->is_valid() ; INDICS_IT->go_next() )
   {
      PDE_AdaptationIndicator const* indic = 
          static_cast< PDE_AdaptationIndicator* >( INDICS_IT->item() ) ;
          
      //??? on pourrait faire tous les calculs ensemble
      double bf_indicator = basis_function_indicator( indic, bf, cell ) ;

      result = result && indic->to_be_unrefined( bf_indicator, 
                                        cell->polyhedron(), 
                                        elm, 
                                        local_node ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------------
double
PDE_AdaptationRequestFromIndicator:: basis_function_indicator(
                                         PDE_AdaptationIndicator const* indic,
                                         PDE_BasisFunctionCell const* bf,
                                         PDE_CellFE const* cell ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequestFromIndicator:: basis_function_indicator" ) ;

   double tot_mes = 0.0 ;
   double result = 0.0 ;
   for( size_t i=0 ; i<bf->nb_cells() ; ++i )
   {
      get_cell_contribution( indic, bf->cell( i ), tot_mes, result ) ;
   }
   result /= tot_mes ;

   return( result ) ;
}

//----------------------------------------------------------------------------
void
PDE_AdaptationRequestFromIndicator:: get_cell_contribution( 
                                        PDE_AdaptationIndicator const* indic,
                                        PDE_CellFE const* cell,
                                        double& total_measure,
                                        double& total_indicator ) const
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequestFromIndicator:: get_cell_contribution" ) ;

   if( cell->is_active() )
   {
      double mes = cell->polyhedron()->measure() ;
      total_measure += mes ;

      size_t id = cell->id_number() ;
      total_indicator += indic->cell_indicator( id ) * mes ;
      return ;
   }
   else
   {
      for( size_t i=0 ; i<cell->nb_childs() ; ++i )
      {
         get_cell_contribution( indic, cell->child( i ), 
                                total_measure, total_indicator ) ;
      }
   }
}

