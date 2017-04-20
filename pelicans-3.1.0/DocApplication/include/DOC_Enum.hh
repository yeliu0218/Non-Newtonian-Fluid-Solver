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

#ifndef ENUM_HH
#define ENUM_HH
#include <PEL_List.hh>

#include <DOC_ClassItem.hh>

// DOC_Class DOC_Enum.

class DOC_Category ;
class DOC_Type ;
class PEL_EXPORT DOC_Enum : public DOC_ClassItem
{
   public: //-------------------------------------------------

   //-- Status
      virtual std::string const& name( void ) const ;
      virtual bool declare( std::string const& a__name ) const ;
      
      
   //-- Input - Output
      virtual std::string prototype( DOC_Writer& sullizer ) const ;
      virtual void display( std::ostream& out ) const ;
      virtual std::string signature( void ) const ;
      
      // Creation
      static DOC_Enum* create( std::string const& a_name,
                           Protection maprotection,
                           DOC_Category const* macategory,
                           DOC_Text * mycomment,
                           PEL_List const* lst  ) ;
      
   protected: //-----------------------------------------------
   private: //-------------------------------------------------
      DOC_Enum( std::string const& a_name,
            Protection maprotection,
            DOC_Category const* macategory,
            DOC_Text * mycomment,
            PEL_List const* lst ) ;
      ~DOC_Enum( void ) ;

      DOC_Enum( void ) ;
      DOC_Enum( DOC_Enum const& other ) ;
      DOC_Enum& operator=( DOC_Enum const& other ) ;
      
      std::string myName ;
      PEL_List const* my_list ;
} ;


#endif
