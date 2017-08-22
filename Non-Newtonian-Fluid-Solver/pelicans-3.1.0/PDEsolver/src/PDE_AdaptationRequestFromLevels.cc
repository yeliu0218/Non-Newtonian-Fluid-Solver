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

#include <PDE_AdaptationRequestFromLevels.hh>

#include <PEL.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Vector.hh>

#include <PDE_AdaptationIndicator.hh>
#include <PDE_BasisFunctionCell.hh>
#include <PDE_CellFE.hh>
#include <PDE_ReferenceElement.hh>

#include <iostream>

using std::cout ;
using std::endl ;

//---------------------------------------------------------------------------
PDE_AdaptationRequestFromLevels*
PDE_AdaptationRequestFromLevels:: create_for_unrefinement(
                                          PEL_Object* a_owner,
                                          size_t target_level,
                                          size_t verbose_level  )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequestFromLevels:: create_for_unrefinement" ) ;
   PEL_CHECK_PRE( target_level != 0 ) ;

   bool refine = false ;
   PDE_AdaptationRequestFromLevels* result =
       new PDE_AdaptationRequestFromLevels( a_owner, target_level, refine,
                                            verbose_level ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_AdaptationRequestFromLevels*
PDE_AdaptationRequestFromLevels:: create_for_refinement(
                                          PEL_Object* a_owner,
                                          size_t target_level,
                                          size_t verbose_level  )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequestFromLevels:: create_for_refinement" ) ;

   bool refine = true ;
   PDE_AdaptationRequestFromLevels* result =
       new PDE_AdaptationRequestFromLevels( a_owner, target_level, refine,
                                            verbose_level ) ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_AdaptationRequestFromLevels:: PDE_AdaptationRequestFromLevels(
                                         PEL_Object* a_owner,
                                         size_t target_level,
                                         bool refine,
                                         size_t verbose_level )
//---------------------------------------------------------------------------
   : PDE_AdaptationRequest( a_owner, PEL::bad_index(), verbose_level )
   , TLEVEL( target_level )
   , REF( refine )
{
}

//---------------------------------------------------------------------------
PDE_AdaptationRequestFromLevels:: ~PDE_AdaptationRequestFromLevels( void  )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
bool
PDE_AdaptationRequestFromLevels:: to_be_refined(
                                       PDE_CellFE const* cell,
                                       PDE_BasisFunctionCell const* bf,
                                       PDE_ReferenceElement const* elm,
                                       size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequestFromLevels:: to_be_refined" ) ;
   PEL_CHECK_PRE( to_be_refined_PRE( cell, bf, elm, local_node ) ) ;

   bool result = false ;
   if( REF )
   {
      PEL_ASSERT( bf->refinement_level() <= TLEVEL ) ;
      if( bf->refinement_level() == TLEVEL ) result = true ;
   }

   return( result ) ;
}

//---------------------------------------------------------------------------
bool
PDE_AdaptationRequestFromLevels:: to_be_unrefined(
                                         PDE_CellFE const* cell,
                                         PDE_BasisFunctionCell const* bf,
                                         PDE_ReferenceElement const* elm,
                                         size_t local_node ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_AdaptationRequestFromLevels:: to_be_unrefined" ) ;
   PEL_CHECK_PRE( to_be_unrefined_PRE( cell, bf, elm, local_node ) ) ;

   bool result = false ;
   if( !REF )
   {
      PEL_ASSERT( bf->refinement_level() <= TLEVEL ) ;
      if( bf->refinement_level() == TLEVEL-1 ) result = true ;
  }

   return( result ) ;
}

