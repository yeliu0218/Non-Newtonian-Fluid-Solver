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

#ifndef PDE_VALUES_ON_MESHING_HH
#define PDE_VALUES_ON_MESHING_HH

#include <PEL_Object.hh>

class PEL_Context ;
class PEL_DataWithContext ;
class PEL_DoubleVector ;
class PEL_ModuleExplorer ;
class PEL_Vector ;

class GE_Color ;
class GE_Mpolyhedron ;
class GE_Point ;
class GE_SetOfPoints ;

class PDE_ReferenceElement ;
class PDE_ResultReader ;

#include <boolVector.hh>
#include <doubleArray2D.hh>

#include <string>

class PEL_EXPORT PDE_ValuesOnMeshing : public PEL_Object
{
   public: //----------------------------------------------------------------

   //-- Instance delivery and initialization

      static PDE_ValuesOnMeshing* create( PEL_Object* a_owner,
                                          std::string const& a_field_name,
                                          size_t nb_field_components,
                                          PEL_ModuleExplorer const* exp,
                                          GE_SetOfPoints const* vs,
                                          PDE_ResultReader const* reader ) ;

   //-- Description

      std::string const& field_name( void ) const ;
      
      size_t nb_components( void ) const ;
      
   //-- Values
      
      void set_mesh( GE_Mpolyhedron const* mesh_poly,
                     GE_Color const* mesh_color ) ;
      
      void compute_value( GE_Point const* pt, doubleVector& result ) const ;

   protected: //-------------------------------------------------------------
      
   private: //---------------------------------------------------------------

      PDE_ValuesOnMeshing( void ) ;
     ~PDE_ValuesOnMeshing( void ) ;
      PDE_ValuesOnMeshing( PDE_ValuesOnMeshing const& other ) ;
      PDE_ValuesOnMeshing& operator=( PDE_ValuesOnMeshing const& other ) ;

      PDE_ValuesOnMeshing( PEL_Object* a_owner,
                           std::string const& a_field_name,
                           size_t nb_field_components,
                           PEL_ModuleExplorer const* exp,
                           GE_SetOfPoints const* vs,
                           PDE_ResultReader const* reader ) ;

      enum FieldInitializationType
      {
         uniformly_defined,
         vertex_defined,
         mesh_defined,
         reader_defined
      } ;

      PEL_DataWithContext const* formula( GE_Color const* color ) ;

   //-- Attributes

      std::string const NAME ;
      size_t const NB_COMPS ;
      FieldInitializationType TYPE ;
      GE_SetOfPoints const* const VM ;
      PDE_ResultReader const* const READER ;
      PEL_Vector* const COLORS ;
      PEL_Vector* const FORMULAS ;
      boolVector USED_COLORS ;
      PEL_DataWithContext const* DEFAULT_FORMULA ;
      PEL_Context* CTX ;
      PEL_DoubleVector* COORDS ;
      GE_Mpolyhedron const* POLY ;
      PEL_DataWithContext const* CURRENT_FORMULA ;
      PDE_ReferenceElement const* ELM ;
      doubleArray2D VVALS ;
} ;

#endif
