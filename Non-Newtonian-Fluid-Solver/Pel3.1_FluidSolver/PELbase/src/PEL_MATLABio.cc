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

#include <PEL_MATLABio.hh>

#include <PEL.hh>
#include <PEL_assertions.hh>
#include <PEL_Double.hh>
#include <PEL_DoubleArray2D.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_Error.hh>
#include <PEL_Int.hh>
#include <PEL_IntArray2D.hh>
#include <PEL_IntVector.hh>
#include <PEL_Module.hh>

#include <doubleArray2D.hh>
#include <doubleVector.hh>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

struct PEL_MATLABio_ERROR
{
   static void n0( void ) ;
   static void n1( void ) ;
   static void n2( std::string const& file_name ) ;
   static void n3( std::string const& file_name ) ;
   static void n4( std::string const& file_name, size_t i_cycle ) ;
} ;

//----------------------------------------------------------------------------
void
PEL_MATLABio:: create_file( std::string const& file_name,
                         MATLAB_FORMAT file_format )
//----------------------------------------------------------------------------
{
    PEL_LABEL( "PEL_MATLABio:: create_file" ) ;
    PEL_CHECK_PRE( !file_name.empty() ) ;
    PEL_CHECK_PRE( file_format != Unspecified ) ;

    std::ofstream file ;
    open_file( file_name, file_format, false, file ) ;
    file.close() ;
}

//----------------------------------------------------------------------
void
PEL_MATLABio:: write_int_variable( std::string const& file_name,
                                MATLAB_FORMAT file_format,
                                std::string const& name,
                                int val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: write_int_variable" ) ;
   PEL_CHECK_PRE( !file_name.empty() ) ;
   PEL_CHECK_PRE( file_format != Unspecified ) ;
   PEL_CHECK_PRE( !name.empty() ) ;

   int const IDIM  = -1 ;
   int const IEND  = 1  ;
   double* XVAR = new double [ IEND ] ;
   *XVAR = val;

   write_variable( file_name, file_format,
                   name,
                   true, IDIM, IEND, XVAR ) ;

   delete [] XVAR ;
}

//----------------------------------------------------------------------
void
PEL_MATLABio:: write_double_variable( std::string const& file_name,
                                   MATLAB_FORMAT file_format,
                                   std::string const& name,
                                   double val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: write_double_variable" ) ;
   PEL_CHECK_PRE( !file_name.empty() ) ;
   PEL_CHECK_PRE( file_format != Unspecified ) ;
   PEL_CHECK_PRE( !name.empty() ) ;

   int const IDIM = -1 ;
   int const IEND = 1 ;
   double* XVAR = new double [ IEND ] ;
   *XVAR = val ;

   write_variable( file_name, file_format,
                   name,
                   false, IDIM, IEND, XVAR ) ;

   delete [] XVAR ;
}

//----------------------------------------------------------------------
void
PEL_MATLABio:: write_gridXY_variable(
                         std::string const& file_name,
                         MATLAB_FORMAT file_format, int dimension,
                         std::string const& name,
                         doubleVector const& valx, doubleVector const& valy, doubleVector const& valz )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: write_doubleXY_variable" ) ;
   PEL_CHECK_PRE( !file_name.empty() ) ;
   PEL_CHECK_PRE( file_format == Text ) ;
   PEL_CHECK_PRE( !name.empty() ) ;

   std::string file_nameXYZ;

   if( dimension == 1){
	   file_nameXYZ = file_name+name+"X.dat";
   }
   if( dimension == 2 ){
	   file_nameXYZ = file_name+name+"XY.dat";
   }

   if( dimension == 3 ){
	   file_nameXYZ = file_name+name+"XYZ.dat";
   }

   write_variableXYZ( file_nameXYZ, file_format, dimension,
		   "X", "Y", "Z", valx, valy, valz );
}

//----------------------------------------------------------------------
void
PEL_MATLABio:: write_doubleVector_variable( std::string const& file_name,
                                         MATLAB_FORMAT file_format,
                                         std::string const& name,
                                         doubleVector const& val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: write_doubleVector_variable" ) ;
   PEL_CHECK_PRE( !file_name.empty() ) ;
   PEL_CHECK_PRE( file_format != Unspecified ) ;
   PEL_CHECK_PRE( !name.empty() ) ;

   int const IDIM = 0 ;
   int const IEND = val.size() ;
   double* XVAR = new double [ IEND ] ;
   for( size_t i=0 ; i<val.size() ; i++ ) XVAR[i] = val(i) ;

   write_variable( file_name, file_format,
                   name,
                   false, IDIM, IEND, XVAR ) ;

   delete [] XVAR ;
}

//----------------------------------------------------------------------
void
PEL_MATLABio:: write_doubleVector_variableUU( std::string const& file_name,
                                         MATLAB_FORMAT file_format,
                                         std::string const& name,
                                         doubleVector const& val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: write_doubleVector_variable" ) ;
   PEL_CHECK_PRE( !file_name.empty() ) ;
   PEL_CHECK_PRE( file_format != Unspecified ) ;
   PEL_CHECK_PRE( !name.empty() ) ;

   int const IDIM = 0 ;
   int const IEND = val.size() ;
   double* XVAR = new double [ IEND ] ;
   for( size_t i=0 ; i<val.size() ; i++ ) XVAR[i] = val(i) ;

   write_variableUU( file_name, file_format,
                   name,
                   false, IDIM, IEND, XVAR ) ;

   delete [] XVAR ;
}

//----------------------------------------------------------------------
void
PEL_MATLABio:: write_doubleVector_variableSTRESS( std::string const& file_name,
                                         MATLAB_FORMAT file_format,
                                         std::string const& name,
                                         doubleVector const& val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: write_doubleVector_variable" ) ;
   PEL_CHECK_PRE( !file_name.empty() ) ;
   PEL_CHECK_PRE( file_format != Unspecified ) ;
   PEL_CHECK_PRE( !name.empty() ) ;

   int const IDIM = 0 ;
   int const IEND = val.size() ;
   double* XVAR = new double [ IEND ] ;
   for( size_t i=0 ; i<val.size() ; i++ ) XVAR[i] = val(i) ;

   write_variableSTRESS( file_name, file_format,
                   name,
                   false, IDIM, IEND, XVAR ) ;

   delete [] XVAR ;
}

//----------------------------------------------------------------------
void
PEL_MATLABio:: write_intVector_variable( std::string const& file_name,
                                      MATLAB_FORMAT file_format,
                                      std::string const& name,
                                      intVector const& val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: write_intVector_variable" ) ;
   PEL_CHECK_PRE( !file_name.empty() ) ;
   PEL_CHECK_PRE( file_format != Unspecified ) ;
   PEL_CHECK_PRE( !name.empty() ) ;

   int const IDIM = 0 ;
   int const IEND = val.size() ;
   double* XVAR = new double [ IEND ] ;
   for( size_t i=0 ; i<val.size() ; i++ ) XVAR[i] = val(i) ;

   write_variable( file_name, file_format,
                   name,
                   true, IDIM, IEND, XVAR ) ;

   delete [] XVAR ;
}

//----------------------------------------------------------------------
void
PEL_MATLABio:: write_doubleArray2D_variable( std::string const& file_name,
                                          MATLAB_FORMAT file_format,
                                          std::string const& name,
                                          doubleArray2D const& val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: write_doubleArray2D_variable" ) ;
   PEL_CHECK_PRE( !file_name.empty() ) ;
   PEL_CHECK_PRE( file_format != Unspecified ) ;
   PEL_CHECK_PRE( !name.empty() ) ;

   int const IDIM = val.index_bound(0) ;
   int const IEND = val.index_bound(0)*val.index_bound(1) ;
   double* XVAR = new double [ IEND ] ;
   int idx = 0 ;
   for( size_t j=0 ; j<val.index_bound(1) ; j++ )
      for( size_t i=0 ; i<val.index_bound(0); i++ )
      {
	 XVAR[idx++] = val(i,j) ;
      }

   write_variable( file_name, file_format,
                   name,
                   false, IDIM, IEND, XVAR ) ;

   delete [] XVAR ;
}

//----------------------------------------------------------------------
void
PEL_MATLABio:: write_intArray2D_variable( std::string const& file_name,
                                       MATLAB_FORMAT file_format,
                                       std::string const& name,
                                       intArray2D const& val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: write_intArray2D_variable" ) ;
   PEL_CHECK_PRE( !file_name.empty() ) ;
   PEL_CHECK_PRE( file_format != Unspecified ) ;
   PEL_CHECK_PRE( !name.empty() ) ;

   int const IDIM = val.index_bound(0) ;
   int const IEND = val.index_bound(0)*val.index_bound(1) ;
   double* XVAR = new double [ IEND ] ;
   int idx = 0 ;
   for( size_t j=0 ; j<val.index_bound(1) ; j++ )
      for( size_t i=0 ; i<val.index_bound(0); i++ )
      {
	 XVAR[idx++] = val(i,j) ;
      }

   write_variable( file_name, file_format,
                   name,
                   true, IDIM, IEND, XVAR ) ;

   delete [] XVAR ;
}

//----------------------------------------------------------------------------
PEL_MATLABio::MATLAB_FORMAT
PEL_MATLABio:: file_format( std::string const& file_name )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: file_format" ) ;
   PEL_CHECK_PRE( !file_name.empty() ) ;

   std::ifstream file ;
   MATLAB_FORMAT result = Unspecified ;

   open_file( file_name, result, file ) ;
   file.close() ;

   return( result ) ;
}


//----------------------------------------------------------------------
void
PEL_MATLABio:: convert_to_local( MATLAB_FORMAT file_format, int* val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: convert_to_local(int)" ) ;

}

//----------------------------------------------------------------------
void
PEL_MATLABio:: convert_to_local( MATLAB_FORMAT file_format, float* val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: convert_to_local(float)" ) ;

}

//----------------------------------------------------------------------
void
PEL_MATLABio:: check_file( std::ofstream const& file )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: check_file(output)" ) ;

   if( !file )
   {
      PEL_MATLABio_ERROR:: n0() ;
   }

   PEL_CHECK_POST( file ) ;
}

//----------------------------------------------------------------------
void
PEL_MATLABio:: open_file( std::string const& file_name,
                       MATLAB_FORMAT file_format,
                       bool append,
                       std::ofstream& file )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: open_file(output)" ) ;
   PEL_CHECK( !file_name.empty() ) ;
   PEL_CHECK( file_format != Unspecified ) ;

   if( file.is_open() )
   {
      PEL_MATLABio_ERROR:: n3( file_name ) ;
   }
   std::ios_base::openmode flags = std::ios_base::out ;
   if( file_format != Text ) flags |= std::ios_base::binary ;
   if( append ) flags |= std::ios_base::app ;
   if(file_name == "saveXY.dat")
	   flags &= std::ios_base::trunc;

   file.open( file_name.c_str(), flags ) ;
   if( !file.is_open() )
   {
      PEL_MATLABio_ERROR:: n2( file_name ) ;
   }

   PEL_CHECK_POST( file.is_open() ) ;
}


//----------------------------------------------------------------------
void
PEL_MATLABio:: write_variableXYZ(std::string const& file_name,
        MATLAB_FORMAT file_format, int dimension,
        std::string const& nameX, std::string const& nameY, std::string const& nameZ,
        doubleVector const& valx, doubleVector const& valy, doubleVector const& valz)
//----------------------------------------------------------------------
{
  PEL_LABEL( "PEL_MATLABwriter:: PEL_MATLABio:: write_variable" ) ;
  PEL_CHECK_PRE( !file_name.empty() ) ;
  PEL_CHECK_PRE( file_format == Text ) ;
  PEL_CHECK_PRE( valx.size() != 0 ) ;

  int const endX = valx.size() ;
  double* XVAR = new double [ endX ] ;
  for( size_t i=0 ; i<valx.size() ; i++ ) XVAR[i] = valx(i) ;

  int const endY = valy.size() ;
  double* YVAR = new double [ endY ] ;
  for( size_t j=0 ; j<valy.size() ; j++ ) YVAR[j] = valy(j) ;

  int const endZ = valz.size() ;
  double* ZVAR = new double [ endY ] ;
  for( size_t k=0 ; k<valz.size() ; k++ ) ZVAR[k] = valz(k) ;

  std::ofstream file ;
  open_file( file_name, file_format, true, file ) ;

  size_t NB_PER_ROW = 1;//( is_integer ? 10 : 5 ) ;

  write_ASCII_int( file, endX );
  file << "   " << nameX.c_str()  ;
  if(dimension >= 2){
	file << "   " << nameY.c_str();
	//<< "   " ;
	//write_ASCII_int( file, endY ) ;
  }
  if(dimension >= 3){
	file << "   " << nameZ.c_str();
	//<< "   "  ;
	//write_ASCII_int( file, endZ ) ;
  }
  file << std::endl ;
  size_t i ;
  file << " " ;
  for( i=0 ; i<(size_t) endX ; ++i )
  {
	  file << " " ;
      write_ASCII_double( file, XVAR[i] ) ;
      if(dimension >= 2 && endX == endY )
    	  write_ASCII_double( file, YVAR[i] ) ;
      if(dimension >= 3 && endX == endZ )
    	  write_ASCII_double( file, ZVAR[i] ) ;

      if( (i+1)%NB_PER_ROW == 0 )
      {
        file << std::endl ;
        if( i+1 != (size_t) endX )
        {
          file << " " ;
        }
      }
  }
  if( ! ( i % NB_PER_ROW ) == 0 ) file << std::endl ;

  file.close() ;
}

//----------------------------------------------------------------------
void
PEL_MATLABio:: write_variable( std::string const& file_name,
                            MATLAB_FORMAT file_format,
                            std::string const& name,
                            bool is_integer,
                            int dim, int end,
                            double* XVAR )
//----------------------------------------------------------------------

{
  PEL_LABEL( "PEL_MATLABwriter:: PEL_MATLABio:: write_variable" ) ;
  PEL_CHECK_PRE( !file_name.empty() ) ;
  PEL_CHECK_PRE( file_format == Text ) ;
  PEL_CHECK_PRE( XVAR != 0 ) ;

  std::ofstream file ;
  open_file( file_name, file_format, true, file ) ;

  static char CNAME[5] ;
  for( size_t i=0 ; i<4 ; i++ ) CNAME[i] = ' ' ;
  CNAME[4] = '\0' ;
  for( size_t i=0 ; i< PEL::min( name.length(), (size_t)4 ) ; i++ )
     CNAME[i] = name[i] ;

  size_t NB_PER_ROW = 1;//( is_integer ? 10 : 5 ) ;
  if( dim == -1 )
  {
    file << "'" << std::setw( 4 ) << CNAME << "' 0 " << std::endl ;
    file << "  " ;
    if( is_integer )
    {
      write_ASCII_int( file, (int) *XVAR ) ;
    }
    else
    {
      write_ASCII_double( file, *XVAR ) ;
    }
    file << std::endl ;
  }
  else if( dim == 0 )
  {
    file << "'" << std::setw( 4 ) << CNAME << "' 1 " << "    " ;
    write_ASCII_int( file, end ) ;
    file << std::endl ;
    size_t i ;
    file << " " ;
    for( i=0 ; i<(size_t) end ; ++i )
    {
      file << " " ;
      if( is_integer )
      {
        write_ASCII_int( file, (int) XVAR[i] ) ;
      }
      else
      {
        write_ASCII_double( file, XVAR[i] ) ;
      }
      if( (i+1)%NB_PER_ROW == 0 )
      {
        file << std::endl ;
        if( i+1 != (size_t) end )
        {
          file << " " ;
        }
      }
    }
    if( ! ( i % NB_PER_ROW ) == 0 ) file << std::endl ;
  }
  else
  {
    size_t nbcol = (size_t) end/dim ;
    file << "'" << std::setw( 4 ) << CNAME << "' 2 "
    << "    " ;
    write_ASCII_int( file, dim ) ;
    file << "    " ;
    write_ASCII_int( file, nbcol ) ;
    file << std::endl ;
    file << " " ;
    size_t idx = 0 ;
    for( size_t i=0 ; i<(size_t) dim ; ++i )
    {
      size_t j = 0 ;

      for(  ; j<nbcol ; ++j )
      {
        file << " " ;
        if( is_integer )
        {
          write_ASCII_int( file, (int) XVAR[idx++] ) ;
        }
        else
        {
          write_ASCII_double( file, XVAR[idx++] ) ;
        }
        if( (j+1)%NB_PER_ROW == 0 )
        {
          file << std::endl ;
          if( j+1 != (size_t) dim )
          {
            file << " " ;
          }
        }
      }
      if( ! ( j % NB_PER_ROW  ) == 0 ) file << std::endl ;
    }
  }

  file.close() ;
}

//----------------------------------------------------------------------
void
PEL_MATLABio:: write_variableUU( std::string const& file_name,
                            MATLAB_FORMAT file_format,
                            std::string const& name,
                            bool is_integer,
                            int dim, int end,
                            double* XVAR )
//----------------------------------------------------------------------

{
  PEL_LABEL( "PEL_MATLABwriter:: PEL_MATLABio:: write_variable" ) ;
  PEL_CHECK_PRE( !file_name.empty() ) ;
  PEL_CHECK_PRE( file_format == Text ) ;
  PEL_CHECK_PRE( XVAR != 0 ) ;

  std::ofstream file ;
  open_file( file_name, file_format, true, file ) ;

  static char CNAME[5] ;
  for( size_t i=0 ; i<4 ; i++ ) CNAME[i] = ' ' ;
  CNAME[4] = '\0' ;
  for( size_t i=0 ; i< PEL::min( name.length(), (size_t)4 ) ; i++ )
     CNAME[i] = name[i] ;
     
  size_t NB_PER_ROW = 2;//( is_integer ? 10 : 5 ) ;

  if( dim == -1 )
  {
    file << "'" << std::setw( 4 ) << CNAME << "' 0 " << std::endl ;
    file << "  " ;
    if( is_integer )
    {
      write_ASCII_int( file, (int) *XVAR ) ;
    }
    else
    {
      write_ASCII_double( file, *XVAR ) ;
    }
    file << std::endl ;
  }
  else if( dim == 0 )
  {
    file << "'" << std::setw( 4 ) << CNAME << "' 1 " << "    " ;
    int tmp=end/2; 
    write_ASCII_int( file, tmp ) ;
    file << std::endl ;
    size_t i ;
    file << " " ;
    for( i=0 ; i<(size_t) end ; ++i )
    {
      file << " " ;
      if( is_integer )
      {
        write_ASCII_int( file, (int) XVAR[i] ) ;
      }
      else
      {
        write_ASCII_double( file, XVAR[i] ) ;
      }
      if( (i+1)%NB_PER_ROW == 0 )
      {
        file << std::endl ;
        if( i+1 != (size_t) end )
        {
          file << " " ;
        }
      }
    }
    if( ! ( i % NB_PER_ROW ) == 0 ) file << std::endl ;
  }
  else
  {
    size_t nbcol = (size_t) end/dim ;
    file << "'" << std::setw( 4 ) << CNAME << "' 2 "
    << "    " ;
    write_ASCII_int( file, dim ) ;
    file << "    " ;
    write_ASCII_int( file, nbcol ) ;
    file << std::endl ;
    file << " " ;
    size_t idx = 0 ;
    for( size_t i=0 ; i<(size_t) dim ; ++i )
    {
      size_t j = 0 ;

      for(  ; j<nbcol ; ++j )
      {
        file << " " ;
        if( is_integer )
        {
          write_ASCII_int( file, (int) XVAR[idx++] ) ;
        }
        else
        {
          write_ASCII_double( file, XVAR[idx++] ) ;
        }
        if( (j+1)%NB_PER_ROW == 0 )
        {
          file << std::endl ;
          if( j+1 != (size_t) dim )
          {
            file << " " ;
          }
        }
      }
      if( ! ( j % NB_PER_ROW  ) == 0 ) file << std::endl ;
    }
  }

  file.close() ;
}


//----------------------------------------------------------------------
void
PEL_MATLABio:: write_variableSTRESS( std::string const& file_name,
                            MATLAB_FORMAT file_format,
                            std::string const& name,
                            bool is_integer,
                            int dim, int end,
                            double* XVAR )
//----------------------------------------------------------------------

{
  PEL_LABEL( "PEL_MATLABwriter:: PEL_MATLABio:: write_variable" ) ;
  PEL_CHECK_PRE( !file_name.empty() ) ;
  PEL_CHECK_PRE( file_format == Text ) ;
  PEL_CHECK_PRE( XVAR != 0 ) ;

  std::ofstream file ;
  open_file( file_name, file_format, true, file ) ;

  static char CNAME[5] ;
  for( size_t i=0 ; i<4 ; i++ ) CNAME[i] = ' ' ;
  CNAME[4] = '\0' ;
  for( size_t i=0 ; i< PEL::min( name.length(), (size_t)4 ) ; i++ )
     CNAME[i] = name[i] ;

  size_t NB_PER_ROW = 3;//( is_integer ? 10 : 5 ) ;

  if( dim == -1 )
  {
    file << "'" << std::setw( 4 ) << CNAME << "' 0 " << std::endl ;
    file << "  " ;
    if( is_integer )
    {
      write_ASCII_int( file, (int) *XVAR ) ;
    }
    else
    {
      write_ASCII_double( file, *XVAR ) ;
    }
    file << std::endl ;
  }
  else if( dim == 0 )
  {
    file << "'" << std::setw( 4 ) << CNAME << "' 1 " << "    " ;
    int tmp=end/2; 
    write_ASCII_int( file, tmp ) ;
    file << std::endl ;
    size_t i ;
    file << " " ;
    for( i=0 ; i<(size_t) end ; ++i )
    {
      file << " " ;
      if( is_integer )
      {
        write_ASCII_int( file, (int) XVAR[i] ) ;
      }
      else
      {
        write_ASCII_double( file, XVAR[i] ) ;
      }
      if( (i+1)%NB_PER_ROW == 0 )
      {
        file << std::endl ;
        if( i+1 != (size_t) end )
        {
          file << " " ;
        }
      }
    }
    if( ! ( i % NB_PER_ROW ) == 0 ) file << std::endl ;
  }
  else
  {
    size_t nbcol = (size_t) end/dim ;
    file << "'" << std::setw( 4 ) << CNAME << "' 2 "
    << "    " ;
    write_ASCII_int( file, dim ) ;
    file << "    " ;
    write_ASCII_int( file, nbcol ) ;
    file << std::endl ;
    file << " " ;
    size_t idx = 0 ;
    for( size_t i=0 ; i<(size_t) dim ; ++i )
    {
      size_t j = 0 ;

      for(  ; j<nbcol ; ++j )
      {
        file << " " ;
        if( is_integer )
        {
          write_ASCII_int( file, (int) XVAR[idx++] ) ;
        }
        else
        {
          write_ASCII_double( file, XVAR[idx++] ) ;
        }
        if( (j+1)%NB_PER_ROW == 0 )
        {
          file << std::endl ;
          if( j+1 != (size_t) dim )
          {
            file << " " ;
          }
        }
      }
      if( ! ( j % NB_PER_ROW  ) == 0 ) file << std::endl ;
    }
  }

  file.close() ;
}

//----------------------------------------------------------------------
void
PEL_MATLABio:: write_ASCII_int( std::ofstream &file, int value )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: write_ASCII_int" ) ;
   PEL_CHECK( file ) ;
   PEL_CHECK( file.is_open() ) ;

   static int const INT_W = 7 ;
   static size_t const INT_M = (size_t) PEL::pow( 10., (double) INT_W ) ;
   if( PEL::abs( value )>=INT_M )
   {
      std::string msg =
         "*** PEL_MATLABio error :\n"
         "    too big int for \"text\" format\n";
      PEL_Error::object()->raise_plain( msg ) ;
   }
   file << std::setw( INT_W ) << value ;
}

//----------------------------------------------------------------------
void
PEL_MATLABio:: write_ASCII_double( std::ofstream &file, double value )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: write_ASCII_double" ) ;
   PEL_CHECK( file ) ;
   PEL_CHECK( file.is_open() ) ;

   static int const FLT_W = 15 ;
   static int const FLT_PREC = 7 ;

   file.precision( FLT_PREC ) ;
   file << std::setw( FLT_W )
        << std::setiosflags( std::ios::scientific )
        << std::setiosflags( std::ios::uppercase )
        << value ;
}

//----------------------------------------------------------------------
void
PEL_MATLABio:: check_file( std::ifstream const& file )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: check_file(input)" ) ;

   if( !file )
   {
      PEL_MATLABio_ERROR:: n1() ;
   }

   PEL_CHECK_POST( file ) ;
}

//----------------------------------------------------------------------
void
PEL_MATLABio:: open_file( std::string const& file_name,
                       MATLAB_FORMAT& file_format,
                       std::ifstream& file )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: open_file(input)" ) ;
   PEL_CHECK( !file_name.empty() ) ;
   PEL_CHECK( file_format != Unspecified ) ;

   using namespace std ; // for strncmp

   if( file.is_open() )
   {
      PEL_MATLABio_ERROR:: n3( file_name ) ;
   }

   // Try to open the file :
   //file.open( file_name.c_str(), std::ios_base::in | std::ios_base::binary ) ;
   file.open( file_name.c_str(), std::ios_base::in ) ;
   if( !file.is_open() )
   {
      PEL_MATLABio_ERROR:: n2( file_name ) ;
   }

   bool ok = false ;

   int len ;

   // Try ASCII FORMAT :
   if( !ok )
   {
      file_format = Text ;
      file.open( file_name.c_str(), std::ios_base::in ) ;
      std::string format ;
      check_file( file ) ;
      file >> format >> len ;
      ok = file && ( len == 4 ) && ( format == "'GENE'" ) ;
      if( !ok ) file.close() ;
   }

   // Test :
   if( !ok )
   {
      PEL_MATLABio_ERROR:: n1() ;
   }

   if( !file.is_open() )
   {
      PEL_MATLABio_ERROR:: n2( file_name ) ;
   }

   PEL_CHECK_POST( file.is_open() ) ;
   PEL_CHECK_POST( file ) ;
   PEL_CHECK_POST( file_format == Text ) ;
}

//----------------------------------------------------------------------
void
PEL_MATLABio:: close_file( std::ifstream& file )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_MATLABio:: close_file" ) ;
   file.close() ;
}


//internal--------------------------------------------------------------
void
PEL_MATLABio_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
         "*** PEL_MATLABio error :\n"
         "    error writing MATLAB file" ) ;
}

//internal--------------------------------------------------------------
void
PEL_MATLABio_ERROR:: n1( void )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
         "*** PEL_MATLABio error :\n"
         "    error reading MATLAB file" ) ;
}

//internal--------------------------------------------------------------
void
PEL_MATLABio_ERROR:: n2( std::string const& file_name )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
         "*** PEL_MATLABio error :\n"
         "    Unable to open : \""+file_name+"\"" ) ;
}

//internal--------------------------------------------------------------
void
PEL_MATLABio_ERROR:: n3( std::string const& file_name )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
         "*** PEL_MATLABio error :\n"
         "    File already open : \""+file_name+"\"" ) ;
}

//internal--------------------------------------------------------------
void
PEL_MATLABio_ERROR:: n4( std::string const& file_name, size_t i_cycle )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** PEL_MATLABio error:" << std::endl
       << "    Reading MATLAB file: \""+file_name+"\"" << std::endl
       << "    Unable to restore cycle number: " << i_cycle ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
