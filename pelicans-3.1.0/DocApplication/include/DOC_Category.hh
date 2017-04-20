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

#ifndef CATEGORY_HH
#define CATEGORY_HH
#include <PEL_Object.hh>

#include <string>

class PEL_List ;

// Category class.
// Each category can be associated a ranking value by adding its value
//  between parenthesis after the category name.
// Default rank is 1000.0.

class PEL_EXPORT DOC_Category : public PEL_Object
{
   public://-------------------------------------------------

   //-- Creation
      static DOC_Category const* create( std::string const& a_name ) ;

   //-- Status
      std::string const& name( void ) const ;
      
   //-- Comparison
      int three_way_comparison( PEL_Object const* other ) const ;
      
   //-- Input - Output(2.5)
      static void display_list( std::ostream& os, size_t indent_width ) ;
      
   protected://----------------------------------------------
   private://------------------------------------------------
      DOC_Category( PEL_Object* a_owner,
                    std::string const& a_name,
                    double a_rank ) ;
      ~DOC_Category( void ) ;
      
      DOC_Category( DOC_Category const& other ) ;
      DOC_Category( void ) ;
      DOC_Category& operator=( DOC_Category const& other ) ;
      int index( void ) const ;
      
      size_t my_index ;
      std::string my_name ;
      static size_t global_index ;
      static PEL_List * categories( void ) ;
      double rank ;
} ;
#endif
