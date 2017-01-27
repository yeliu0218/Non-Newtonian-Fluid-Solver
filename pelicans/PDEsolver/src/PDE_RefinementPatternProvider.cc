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

#include <PDE_RefinementPatternProvider.hh>

#include <PEL_Error.hh>
#include <PEL_Root.hh>

#include <GE_Mpolyhedron.hh>
#include <GE_ReferencePolyhedron.hh>
#include <GE_ReferencePolyhedronRefiner.hh>

#include <PDE_CellFE.hh>
#include <PDE_ReferenceElement.hh>
#include <PDE_ReferenceElementRefiner.hh>

#include <iostream>
#include <sstream>
#include <string>

using std::map ;
using std::endl ;

//---------------------------------------------------------------------------
PDE_RefinementPatternProvider const*
PDE_RefinementPatternProvider:: object( std::string const& a_name )
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_RefinementPatternProvider:: object" ) ;

   static PDE_RefinementPatternProvider const* result = 0 ;
   if( result == 0 )
   {
      result = new PDE_RefinementPatternProvider( PEL_Root::object() ) ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == PEL_Root::object() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_RefinementPatternProvider:: PDE_RefinementPatternProvider( 
                                                         PEL_Object* a_owner )
//---------------------------------------------------------------------------
   : PEL_Object( a_owner )
{
   add_pattern( 
      PDE_ReferenceElement::object( "PDE_2D_Q0_1node" ),
      PDE_ReferenceElementRefiner::object( "PDE_2D_Q0_1node_RefinerA" ) ) ;
   
   add_pattern( 
         PDE_ReferenceElement::object( "PDE_2D_Q1_4nodes" ),
         PDE_ReferenceElementRefiner::object( "PDE_2D_Q1_4nodes_RefinerA" ) ) ;

   add_pattern( 
         PDE_ReferenceElement::object( "PDE_2D_Q1isoNonConfA_4nodes" ),
         PDE_ReferenceElementRefiner::object( "PDE_2D_Q1isoNonConfA_4nodes_RefinerA" ) ) ;

   add_pattern( 
      PDE_ReferenceElement::object( "PDE_2D_Q2_9nodes" ),
      PDE_ReferenceElementRefiner::object( "PDE_2D_Q2_9nodes_RefinerA" ) ) ;
   
   add_pattern( 
         PDE_ReferenceElement::object( "PDE_2D_P0_1node" ),
         PDE_ReferenceElementRefiner::object( "PDE_2D_P0_1node_RefinerA" ) ) ;
   
   add_pattern( 
      PDE_ReferenceElement::object( "PDE_2D_P1_3nodes" ),
      PDE_ReferenceElementRefiner::object( "PDE_2D_P1_3nodes_RefinerA" ) ) ;

   add_pattern( 
      PDE_ReferenceElement::object( "PDE_2D_P2_6nodes" ),
      PDE_ReferenceElementRefiner::object( "PDE_2D_P2_6nodes_RefinerA" ) ) ;

   add_pattern( 
      PDE_ReferenceElement::object( "PDE_3D_Q1_8nodes" ),
      PDE_ReferenceElementRefiner::object( "PDE_3D_Q1_8nodes_RefinerA" ) ) ;

   add_pattern( 
      PDE_ReferenceElement::object( "PDE_3D_Q2_27nodes" ),
      PDE_ReferenceElementRefiner::object( "PDE_3D_Q2_27nodes_RefinerA" ) ) ;
}

//---------------------------------------------------------------------------
PDE_RefinementPatternProvider:: ~PDE_RefinementPatternProvider( void )
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
GE_ReferencePolyhedronRefiner const*
PDE_RefinementPatternProvider:: cell_refiner( PDE_CellFE const* cell ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_RefinementPatternProvider:: cell_refiner" ) ;
   PEL_CHECK_PRE( cell != 0 ) ;

   GE_ReferencePolyhedron const* ref_poly = 
                          cell->polyhedron()->reference_polyhedron() ;
   IT_C = CELL_RFS.find( ref_poly ) ;

   if( IT_C == CELL_RFS.end() )
   {
      std::ostringstream mesg ;
      mesg << "*** PDE_RefinementPatternProvider:" << std::endl ;
      mesg << "    cells with \"" << ref_poly->name() << "\"" << std::endl ; 
      mesg << "    as reference polyhedron are not handled" ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   GE_ReferencePolyhedronRefiner const* result = IT_C->second ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == PEL_Root::object() ) ;
   return( result ) ;
}

//---------------------------------------------------------------------------
PDE_ReferenceElementRefiner const* 
PDE_RefinementPatternProvider:: reference_element_refiner(
                                      PDE_ReferenceElement const* elm ) const
//---------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_RefinementPatternProvider:: reference_element_refiner" ) ;
   PEL_CHECK_PRE( elm != 0 ) ;

   IT_E = ELEM_RFS.find( elm ) ;

   if( IT_E == ELEM_RFS.end() )
   {
      std::ostringstream mesg ;
      mesg << "*** PDE_RefinementPatternProvider:" << std::endl ;
      mesg << "    the reference element \"" << elm->name() 
           << "\"" << std::endl ; 
      mesg << "    is not handled" ;
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   PDE_ReferenceElementRefiner const* result = IT_E->second ;

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == PEL_Root::object() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
void
PDE_RefinementPatternProvider:: add_pattern( 
                                PDE_ReferenceElement const* elm,
                                PDE_ReferenceElementRefiner const* refiner )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PDE_RefinementPatternProvider:: add_pattern" ) ;
   PEL_CHECK( elm != 0 ) ;
   PEL_CHECK( refiner != 0 ) ;

   if( ELEM_RFS.find( elm ) != ELEM_RFS.end() )
   {
      std::ostringstream mesg ;
      mesg << "*** PDE_RefinementPatternProvider:" << std::endl ;
      mesg << "    multiple treatment of the" << std::endl ;
      mesg << "    reference element \"" << elm->name() << "\"" << std::endl ; 
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }

   ELEM_RFS[ elm ] = refiner ;

   GE_ReferencePolyhedron const* ref_poly = elm->reference_polyhedron() ;
   std::map< GE_ReferencePolyhedron const*,
             GE_ReferencePolyhedronRefiner const* >::const_iterator it = 
      CELL_RFS.find( ref_poly ) ;
   if( it == CELL_RFS.end() )
   {
      CELL_RFS[ ref_poly ] = refiner->cell_refiner() ;
   }
   else if( it->second != refiner->cell_refiner() )
   {
      std::ostringstream mesg ;
      mesg << "*** PDE_RefinementPatternProvider:" << std::endl ;
      mesg << "    multiple inconsistent definitions of cell refiners" 
           << std::endl ;
      mesg << "       " << it->second->type_name() << std::endl ;
      mesg << "       " << refiner->cell_refiner()->type_name() << std::endl ; 
      PEL_Error::object()->raise_plain( mesg.str() ) ;
   }
}
