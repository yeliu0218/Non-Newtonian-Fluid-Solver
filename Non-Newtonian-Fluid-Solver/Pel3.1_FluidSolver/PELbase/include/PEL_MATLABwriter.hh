/*
 *  Copyright :
 *    "Institut de Radioprotection et de Sret�Nucl�ire - IRSN" (1995-2008)
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

#ifndef PEL_MATLABwriter_HH
#define PEL_MATLABwriter_HH

#include <PEL_DataOnMeshingWriter.hh>

#include <PEL_MATLABio.hh>
class PEL_ContextSimple ;
class PEL_DoubleVector ;
class PEL_Data ;

/*
writers for the data postprocessor MATLAB

PUBLISHED
*/

class PEL_MATLABwriter : public PEL_DataOnMeshingWriter
{
   public: //-----------------------------------------------------------

   //-- Write

      virtual void write_cycle( PEL_ModuleExplorer const* exp ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~PEL_MATLABwriter( void ) ;
      PEL_MATLABwriter( PEL_MATLABwriter const& other ) ;
      PEL_MATLABwriter& operator=( PEL_MATLABwriter const& other ) ;

      PEL_MATLABwriter( PEL_Object* a_owner,
		     PEL_ModuleExplorer const* exp ) ;

   //-- Plug in

      PEL_MATLABwriter( void ) ;

      virtual PEL_MATLABwriter* create_replica(
                                       PEL_Object* a_owner,
				       PEL_ModuleExplorer const* exp ) const ;

   //-- Write

      int nb_variables( PEL_ModuleExplorer const* exp ) const ;

      void write_grid( PEL_ModuleExplorer const* exp ) const ;

      void write_field( PEL_ModuleExplorer const* exp ) const ;

      void write_integration_domain( PEL_ModuleExplorer const* exp ) const ;

      void write_one_variable( std::string const& name,
                               PEL_Data const* val ) const ;


      //-- Output file

      std::string output_file_name( size_t nb, std::string add_string ) ;

   //-- Class attributes

      static PEL_MATLABwriter const* PROTOTYPE ;

   //-- Attributes

      size_t ICYCLE ;
      std::string FILEBASENAME ;
      std::string FILEXTENSION ;
      std::string FILE_NAME ;
      std::string FILE_NAME_U ;
      std::string FILE_NAME_STRESS ;

      PEL_MATLABio::MATLAB_FORMAT FORMAT ;

      PEL_ContextSimple* CONTEXT ;
      PEL_DoubleVector* COORDS;

} ;

#endif
