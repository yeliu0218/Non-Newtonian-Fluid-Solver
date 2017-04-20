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

#include <FE_AdaptationStepHN.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>

#include <PDE_AdapterHN.hh>
#include <PDE_ResultSaver.hh>

#include <FE_TimeIterator.hh>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

using std::endl ;
using std::cout ;
using std::string ;
using std::setprecision ; using std::setw ;

FE_AdaptationStepHN const* 
FE_AdaptationStepHN:: PROTOTYPE = new FE_AdaptationStepHN() ;

//----------------------------------------------------------------------
FE_AdaptationStepHN:: FE_AdaptationStepHN( void )
//----------------------------------------------------------------------
   : FE_OneStepIteration( "FE_AdaptationStepHN" )
{
}

//----------------------------------------------------------------------
FE_AdaptationStepHN*
FE_AdaptationStepHN:: create_replica( PEL_Object* a_owner,
                                      PDE_DomainAndFields const* dom,
                                      FE_SetOfParameters const* prms,
                                      PEL_ModuleExplorer* exp ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationStepHN:: create_replica" ) ;
   PEL_CHECK( create_replica_PRE( a_owner, dom, prms, exp ) ) ;

   FE_AdaptationStepHN* result = 
                        new FE_AdaptationStepHN( a_owner, dom, prms, exp ) ;

   PEL_CHECK( create_replica_POST( result, a_owner, dom, prms, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
FE_AdaptationStepHN:: FE_AdaptationStepHN( PEL_Object* a_owner,
                                           PDE_DomainAndFields const* dom,
                                           FE_SetOfParameters const* prms,
                                           PEL_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : FE_OneStepIteration( a_owner, dom, exp )
   , TA( dom->adapter_HN() )
{
}

//----------------------------------------------------------------------
FE_AdaptationStepHN:: ~FE_AdaptationStepHN( void )
//----------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void
FE_AdaptationStepHN:: do_one_inner_iteration( FE_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationStepHN:: do_one_inner_iteration" ) ;
   PEL_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "FE_AdaptationStepHN:: do_one_inner_iteration" ) ;
   
   if( TA != 0 )
   {
      TA->adapt() ;
   }
   
   stop_total_timer() ;
}

//------------------------------------------------------------------------
void
FE_AdaptationStepHN:: do_additional_savings( FE_TimeIterator const* t_it,
                                             PDE_ResultSaver* rs )
//------------------------------------------------------------------------
{
   PEL_LABEL( "FE_AdaptationStepHN:: do_additional_savings" ) ;
   PEL_CHECK_PRE( do_additional_savings_PRE( t_it, rs ) ) ;

   start_total_timer( "FE_AdaptationStepHN:: do_additional_savings" ) ;

   if( TA != 0 )
   {
      if( t_it->is_started() ) 
      {
         if( verbose_level() >= 2 )
         {
            PEL::out() << indent() << "   save grid..." << endl ;
         }
         rs->save_grid() ;
      }
   }

   stop_total_timer() ;
}
