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

#include <PEL_TICio.hh>

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

struct PEL_TICio_ERROR
{
   static void n0( void ) ;
   static void n1( void ) ;
   static void n2( std::string const& file_name ) ;
   static void n3( std::string const& file_name ) ;
   static void n4( std::string const& file_name, size_t i_cycle ) ;
} ;

//----------------------------------------------------------------------------
void
PEL_TICio:: create_file( std::string const& file_name,
                         TIC_FORMAT file_format )
//----------------------------------------------------------------------------
{
    PEL_LABEL( "PEL_TICio:: create_file" ) ;
    PEL_CHECK_PRE( !file_name.empty() ) ;
    PEL_CHECK_PRE( file_format != Unspecified ) ;

    std::ofstream file ;
    open_file( file_name, file_format, false, file ) ;
    file.close() ;
}

//----------------------------------------------------------------------------
void
PEL_TICio:: write_new_cycle( std::string const& file_name,
                             TIC_FORMAT file_format,
                             size_t nb_vars, size_t i_cycle )
//----------------------------------------------------------------------------
{
    PEL_LABEL( "PEL_TICio:: follows_fortran_convention" ) ;
    PEL_CHECK_PRE( !file_name.empty() ) ;
    PEL_CHECK_PRE( file_format != Unspecified ) ;
    PEL_CHECK_PRE( i_cycle > 0 ) ;
    PEL_CHECK_PRE( nb_vars > 0 ) ;

    std::ofstream file ;
    open_file( file_name, file_format, true, file ) ;

    int const NBVAR = (int) nb_vars ;
    int const ICYCLE = (int) i_cycle ;
    
    if( file )
    {
       if( i_cycle == 1 )
       {
          if( file_format == Text )
          {
             file << "'GENE'        4" << std::endl ;
          }
          else if( file_format == CBinary )
          {
             int k = 4 ;
             char label[] = "GENE" ;
             file.write( label, 4 ) ;
             file.write( (const char*)&k, sizeof(int) ) ;
          }
          else if( file_format == Binary )
          {
             int k = 4 ;
             char label[] = "GENE" ;
             write_length( file, file_format, 4+sizeof(int) ) ;
             file.write( label, 4 ) ;
             file.write( (const char*)&k, sizeof(int) ) ;
             write_length( file, file_format, 4+sizeof(int) ) ;
          }
          else
          {
             PEL_TICio_ERROR:: n0() ;
          }  
       }
       if( file_format == Text )
       {
          file << " " << ICYCLE << " " << NBVAR << std::endl ;
       }
       else if( file_format == CBinary ||
                file_format == Binary ||
                file_format == Binary_no_local )
       {
          write_length( file, file_format, 2*sizeof(int) ) ;
          int c = ICYCLE ;
          convert_to_local( file_format, &c ) ;
          file.write( (const char * ) &c, sizeof(int) ) ;
          int v = NBVAR ;
          convert_to_local( file_format, &v ) ;
          file.write( (const char * ) &v, sizeof(int) ) ;
          write_length( file, file_format, 2*sizeof(int) ) ;
       }
       else
       {
          PEL_TICio_ERROR:: n0() ;
       }
       file.close() ;
   }
}

//----------------------------------------------------------------------
void
PEL_TICio:: write_int_variable( std::string const& file_name,
                                TIC_FORMAT file_format,
                                std::string const& name,
                                int val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: write_int_variable" ) ;
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
PEL_TICio:: write_double_variable( std::string const& file_name,
                                   TIC_FORMAT file_format,
                                   std::string const& name,
                                   double val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: write_double_variable" ) ;
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
PEL_TICio:: write_doubleVector_variable( std::string const& file_name,
                                         TIC_FORMAT file_format,
                                         std::string const& name,
                                         doubleVector const& val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: write_doubleVector_variable" ) ;
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
PEL_TICio:: write_intVector_variable( std::string const& file_name,
                                      TIC_FORMAT file_format,
                                      std::string const& name,
                                      intVector const& val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: write_intVector_variable" ) ;
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
PEL_TICio:: write_doubleArray2D_variable( std::string const& file_name,
                                          TIC_FORMAT file_format,
                                          std::string const& name,
                                          doubleArray2D const& val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: write_doubleArray2D_variable" ) ;
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
PEL_TICio:: write_intArray2D_variable( std::string const& file_name,
                                       TIC_FORMAT file_format,
                                       std::string const& name,
                                       intArray2D const& val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: write_intArray2D_variable" ) ;
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
PEL_TICio::TIC_FORMAT
PEL_TICio:: file_format( std::string const& file_name )
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: file_format" ) ;
   PEL_CHECK_PRE( !file_name.empty() ) ;

   std::ifstream file ;
   TIC_FORMAT result = Unspecified ;
   
   open_file( file_name, result, file ) ;
   file.close() ;

   return( result ) ;
}

//----------------------------------------------------------------------------
PEL_Module*
PEL_TICio:: create_from_gene_file( PEL_Object* a_owner,
                                   std::string const& mod_name,
                                   std::string const& file_name,
                                   size_t i_cycle ) 
//----------------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: create_from_gene_file" ) ;
   PEL_CHECK_PRE( !mod_name.empty() ) ;
   PEL_CHECK_PRE( !file_name.empty() ) ;
   
   PEL_Module* result = PEL_Module::create( a_owner, mod_name ) ;
   
   std::ifstream file ;
   TIC_FORMAT file_format = Unspecified ;

   open_file( file_name, file_format, file ) ;

   bool ok_cycle = false ;
   size_t nb_vars = PEL::bad_index() ;
   size_t i_cyc = PEL::bad_index()  ;
   for(;;)
   {
      read_new_cycle( file, file_format, nb_vars, i_cyc ) ;
      if( nb_vars==0 ) break ;
      PEL_Module* cycle = 0 ;
      if( i_cycle == PEL::bad_index() || i_cycle == i_cyc )
      {
         std::ostringstream os ;
         os << "cycle_" << (int) i_cyc ;
         
         cycle = PEL_Module::create( result, os.str() ) ;
         result->add_module( cycle ) ;
      }
         
      for( size_t i=0 ; i<nb_vars ; ++i )
      {
         size_t type ;
         size_t dim ;
         size_t end ;
         bool is_integer ;
         std::string name ;
         read_new_variable( file, file_format, name, type, dim, end, is_integer ) ;
            
         float* fptr = 0 ;
         int*   iptr = 0 ;
         if( !is_integer )
         {
            fptr = new float[end] ;
            read_float( file, file_format, dim, end, fptr ) ;

            if( cycle != 0 )
            {
               if( type == 0 )
               {
                  cycle->add_entry( name,
                                    PEL_Double::create( cycle, fptr[0] ) ) ;
               }
               else if( type == 1 )
               {
                  doubleVector v(end) ;
                  for( size_t j=0 ; j<end ; ++j )
                  {
                     v(j)=fptr[j] ;
                  }
                  cycle->add_entry( name,
                                    PEL_DoubleVector::create( cycle, v ) ) ;
               }
               else
               {
                  size_t nbcols = end/dim ;
                  doubleArray2D arr( dim, nbcols ) ;
                  size_t cpt=0 ;
                  for( size_t j=0 ; j<nbcols ; ++j )
                  {
                     for( size_t l=0; l<dim ; ++l )
                     {
                        arr(l,j)=fptr[cpt++] ;
                     }
                  }
                  cycle->add_entry( name,
                                    PEL_DoubleArray2D::create( cycle, arr ) ) ;
               }
            }
            delete fptr ; fptr = 0 ;
         }
         else
         {
               
            iptr = new int[end] ;
            read_int( file, file_format, dim, end, iptr ) ;

            if( cycle != 0 )
            {
               if( type == 0 )
               {
                  cycle->add_entry( name, PEL_Int::create( cycle, iptr[0] ) ) ;
               }
               else if( type == 1 )
               {
                  intVector v(end) ;
                  for( size_t j=0 ; j<end ; j++ )
                  {
                     v(j)=iptr[j] ;
                  }
                  cycle->add_entry( name,
                                    PEL_IntVector::create( cycle, v ) ) ;
                  
               }
               else
               {
                  size_t nbcols = end/dim ;
                  intArray2D arr( dim, nbcols ) ;
                  size_t cpt=0 ;
                  for( size_t j=0; j<nbcols ; ++j )
                  {
                     for( size_t l=0 ; l<dim ; ++l )
                     {
                        arr(l,j)=iptr[cpt++] ;
                     }
                  }
                  cycle->add_entry( name,
                                    PEL_IntArray2D::create( cycle, arr ) ) ;
               }
            }
            delete iptr ; iptr = 0 ;
         }
      }
      if( i_cycle == i_cyc )
      {
         ok_cycle = true ;
         break ;
      }
   }
   close_file( file ) ;

   if( i_cycle != PEL::bad_index() && !ok_cycle )
   {
      PEL_TICio_ERROR::n4( file_name, i_cycle ) ;
   }

   PEL_CHECK_POST( result != 0 ) ;
   PEL_CHECK_POST( result->owner() == a_owner ) ;
   PEL_CHECK_POST( result->name() == mod_name ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
PEL_TICio:: convert_to_local( TIC_FORMAT file_format, int* val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: convert_to_local(int)" ) ;

   if( file_format == Binary_no_local )
   {
      static int const sz  = sizeof(int)/2 ;
      static int const max = sizeof(int)-1 ;
   
      char x ;
      char *n = reinterpret_cast<char*>(val) ;
      for( int i=0 ; i<sz ; i++)
      {
         x        = n[i] ;
         n[i]     = n[max-i] ;
         n[max-i] = x ;
      }
   }
}

//----------------------------------------------------------------------
void
PEL_TICio:: convert_to_local( TIC_FORMAT file_format, float* val )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: convert_to_local(float)" ) ;
   
   if( file_format == Binary_no_local )
   {
      static int const sz  = sizeof(float)/2 ;
      static int const max = sizeof(float)-1 ;
   
      char *n = reinterpret_cast<char*>(val), x ;
      for( int i=0 ; i<sz ; i++)
      {
         x        = n[i] ;
         n[i]     = n[max-i] ;
         n[max-i] = x ;
      }
   }
}

//----------------------------------------------------------------------
void
PEL_TICio:: check_file( std::ofstream const& file )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: check_file(output)" ) ;
 
   if( !file )
   {
      PEL_TICio_ERROR:: n0() ;
   }

   PEL_CHECK_POST( file ) ;
}

//----------------------------------------------------------------------
void
PEL_TICio:: open_file( std::string const& file_name,
                       TIC_FORMAT file_format,
                       bool append,
                       std::ofstream& file )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: open_file(output)" ) ;
   PEL_CHECK( !file_name.empty() ) ;
   PEL_CHECK( file_format != Unspecified ) ;
   
   if( file.is_open() )
   {
      PEL_TICio_ERROR:: n3( file_name ) ;
   }
   std::ios_base::openmode flags = std::ios_base::out ;
   if( file_format != Text ) flags |= std::ios_base::binary ;
   if( append ) flags |= std::ios_base::app ;
   
   file.open( file_name.c_str(), flags ) ;
   if( !file.is_open() )
   {
      PEL_TICio_ERROR:: n2( file_name ) ;
   }

   PEL_CHECK_POST( file.is_open() ) ;
}

//----------------------------------------------------------------------------
void
PEL_TICio:: check_variable_name( std::string const& name, bool is_integer )
//----------------------------------------------------------------------------
{
    PEL_LABEL( "PEL_TICio:: check_variable_name" ) ;
   
    bool ok = name.length() > 0 ;
   
    if( ok )
    {
       char type_id = name[0] ;
      
       if( ( type_id>='I' && type_id<='N' ) ||
           ( type_id>='i' && type_id<='n' ) )
       {
          ok = is_integer ;
       }
       else
       {
          ok = !is_integer ;
       }
    }

    if( !ok )
    {
       PEL_Error::object()->raise_plain(
          "*** PEL_TICio error :\n"
          "    the TIC variable \""+name+"\"\n"
          "    of "+( is_integer ? "integer" : "float" )+" value(s)\n"
          "    does not follow the TIC convention\n"
          "      integer : begin with I,J,K,L,M,N,i,j,k,l,m or n\n"
          "      float : otherwise" ) ;
    }
}

//----------------------------------------------------------------------
void
PEL_TICio:: write_variable( std::string const& file_name,
                            TIC_FORMAT file_format,
                            std::string const& name,
                            bool is_integer,
                            int dim, int end,
                            double* XVAR ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICwriter:: PEL_TICio:: write_variable" ) ;
   PEL_CHECK_PRE( !file_name.empty() ) ;
   PEL_CHECK_PRE( file_format != Unspecified ) ;
   PEL_CHECK_PRE( XVAR != 0 ) ;

   check_variable_name( name, is_integer ) ;
   
   std::ofstream file ;
   open_file( file_name, file_format, true, file ) ;

   static char CNAME[5] ;
   for( size_t i=0 ; i<4 ; i++ ) CNAME[i] = ' ' ;
   CNAME[4] = '\0' ;
   for( size_t i=0 ; i< PEL::min( name.length(), (size_t)4 ) ; i++ )
      CNAME[i] = name[i] ;
   
   if( file_format == Text )
   {
      size_t NB_PER_ROW = ( is_integer ? 10 : 5 ) ;
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
   }
   else
   {
      if( dim == -1 )
      {
         int k = 0 ;
         float x = (float) *XVAR ;
         int ix = (int) *XVAR ;
         convert_to_local( file_format, &k ) ;
         convert_to_local( file_format, &x ) ;
         convert_to_local( file_format, &ix ) ;
         
         write_length( file, file_format, 4+sizeof(int) ) ;
         file.write( CNAME, 4 ) ;
         file.write( (const char*)&k, sizeof(int) ) ;
         write_length( file, file_format, 4+sizeof(int) ) ;
         
         int const length = ( is_integer ? sizeof(int) : sizeof(float) ) ;
         const char* v = ( is_integer ? (const char*) &ix : (const char*) &x ) ;
         write_length( file, file_format, length ) ;
         file.write( v, length ) ;
         write_length( file, file_format, length ) ;  
      }
      else if( dim == 0 )
      {
         int k = 1 ;
         int e = end ;
         convert_to_local( file_format, &k ) ;
         convert_to_local( file_format, &e ) ;
         
         write_length( file, file_format, 4+2*sizeof(int) ) ;
         file.write( CNAME, 4 ) ;
         file.write( (const char*)&k, sizeof(int) ) ;
         file.write( (const char*)&e, sizeof(int) ) ;
         write_length( file, file_format, 4+2*sizeof(int) ) ;

         int const length = ( is_integer ? sizeof(int) : sizeof(float) ) ;
         write_length( file, file_format, end*length ) ;
         for( size_t i=0 ; i<(size_t) end ; i++ )
         {
            float x = (float) XVAR[i] ;
            int ix = (int) XVAR[i] ;
            convert_to_local( file_format, &x ) ;
            convert_to_local( file_format, &ix ) ;
            const char* v =
               ( is_integer ? (const char*) &ix : (const char*) &x ) ;
            file.write( v, length ) ;
         }
         write_length( file, file_format, end*length ) ;
      }
      else 
      {
         int nbcol = end/dim ;
         int k = 2 ;
         int d = dim ;
         int c = nbcol ;
         convert_to_local( file_format, &k ) ;
         convert_to_local( file_format, &d ) ;
         convert_to_local( file_format, &c ) ;
         
         write_length( file, file_format, 4+3*sizeof(int) ) ;
         file.write( CNAME, 4 ) ;
         file.write( (const char*)&k, sizeof(int) ) ;
         file.write( (const char*)&d, sizeof(int) ) ;
         file.write( (const char*)&c, sizeof(int) ) ;
         write_length( file, file_format, 4+3*sizeof(int) ) ;

         int const length = ( is_integer ? sizeof(int) : sizeof(float) ) ;
         write_length( file, file_format, end*length ) ;
         size_t idx = 0 ;
         for( size_t i=0 ; i<(size_t) dim ; i++ )
         {
            for( size_t j=0 ; j<(size_t) nbcol ; j++ )
            {
               float x = (float) XVAR[idx] ;
               int ix = (int) XVAR[idx++] ;
               convert_to_local( file_format, &x ) ;
               convert_to_local( file_format, &ix ) ;
               const char* v =
                  ( is_integer ? (const char*) &ix : (const char*) &x ) ;
               file.write( v, length ) ;
            }
         }
         write_length( file, file_format, end*length ) ;
      }
   }
   
   file.close() ; 
}

//----------------------------------------------------------------------
void
PEL_TICio:: write_ASCII_int( std::ofstream &file, int value )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: write_ASCII_int" ) ;
   PEL_CHECK( file ) ;
   PEL_CHECK( file.is_open() ) ;
   
   static int const INT_W = 5 ;
   static size_t const INT_M = (size_t) PEL::pow( 10., (double) INT_W ) ;
   if( PEL::abs( value )>=INT_M )
   {
      std::string msg =
         "*** PEL_TICio error :\n"
         "    too big int for \"text\" format\n"
         "    (use \"binary\" or \"Cbinary\")" ;
      PEL_Error::object()->raise_plain( msg ) ;
   }
   file << std::setw( INT_W ) << value ;
}

//----------------------------------------------------------------------
void
PEL_TICio:: write_ASCII_double( std::ofstream &file, double value ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: write_ASCII_double" ) ;
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
PEL_TICio:: write_length( std::ofstream &file, TIC_FORMAT file_format,
                          int length ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: write_length" ) ;
   PEL_CHECK( file ) ;
   PEL_CHECK( file.is_open() ) ;
   PEL_CHECK( file_format != Unspecified ) ;
   PEL_CHECK( length > 0 ) ;

   if( file_format == Binary || file_format == Binary_no_local )
   {
      int l = length ;
      convert_to_local( file_format, &l ) ;
      check_file( file ) ;
      file.write( (const char*) &l, sizeof( int ) ) ;
   }
}

//----------------------------------------------------------------------
void
PEL_TICio:: check_file( std::ifstream const& file )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: check_file(input)" ) ;
 
   if( !file )
   {
      PEL_TICio_ERROR:: n1() ;
   }

   PEL_CHECK_POST( file ) ;
}

//----------------------------------------------------------------------
void
PEL_TICio:: open_file( std::string const& file_name,
                       TIC_FORMAT& file_format,
                       std::ifstream& file )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: open_file(input)" ) ;
   PEL_CHECK( !file_name.empty() ) ;
   PEL_CHECK( file_format == Unspecified ) ;

   if( file.is_open() )
   {
      PEL_TICio_ERROR:: n3( file_name ) ;
   }

   // Try to open the file :
   file.open( file_name.c_str(), std::ios_base::in | std::ios_base::binary ) ;
   if( !file.is_open() )
   {
      PEL_TICio_ERROR:: n2( file_name ) ;
   }

   bool ok = false ;

   char buff[5] ;
   int len ;
   
   // Try C FORMAT :
   {
      file_format = CBinary ;
      file.read( buff, 4 ) ;
      std::string sbuff( buff, 4 ) ;
      
      ok = (file &&  file.gcount()==4 && sbuff== "GENE" ) ;
      if( ok ) file.read( (char*)&len, sizeof( int ) ) ;
      ok = ok && (file &&  file.gcount()==sizeof( int ) ) && ( len==4 ) ;
      if( !ok ) file.close() ;
   }

   // Try FORTRAN FORMAT :
   if( !ok )
   {
      file_format = Binary ;
      file.open( file_name.c_str(), std::ios_base::in | std::ios_base::binary ) ;

      file.read( (char*)&len,sizeof( int ) ) ;
      
      if( len != 4+sizeof( int ) )
      {
         file_format = Binary_no_local ;
         convert_to_local( file_format, &len ) ;
      }

      ok = file && file.gcount()==sizeof( int ) && ( len==4+sizeof( int ) ) ;
      
      if( ok ) file.read( buff, 4 ) ;
      std::string sbuff( buff, 4 ) ;
      ok = ok && ( file && file.gcount()==4 && sbuff== "GENE" ) ;
      
      if( ok )
      {
         file.read( (char*)&len, sizeof( int ) ) ;
         convert_to_local( file_format, &len ) ;
      }
      ok = ok && ( file && file.gcount()==sizeof( int ) ) && ( len==4 ) ;
      if( ok ) read_length( file, file_format, 4+sizeof( int ) ) ;
      if( !ok ) file.close() ;
   }

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
      PEL_TICio_ERROR:: n1() ;
   }

   if( !file.is_open() )
   {
      PEL_TICio_ERROR:: n2( file_name ) ;
   }

   PEL_CHECK_POST( file.is_open() ) ;
   PEL_CHECK_POST( file ) ;
   PEL_CHECK_POST( file_format == Text ||
                   file_format == CBinary ||
                   file_format == Binary ||
                   file_format == Binary_no_local ) ;
}

//----------------------------------------------------------------------
void
PEL_TICio:: close_file( std::ifstream& file )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: close_file" ) ;
   file.close() ;
}

//----------------------------------------------------------------------
void
PEL_TICio:: read_new_cycle( std::ifstream& file,
                            TIC_FORMAT file_format,
                            size_t& nb_vars, size_t& i_cycle )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: read_new_cycle" ) ;
   PEL_CHECK( file.is_open() ) ;
   PEL_CHECK( file_format != Unspecified ) ;
   
   int NBVARS = 0 ;
   int ICYCLE = 0 ;
   
   if( file_format == Text )
   {
      check_file( file ) ;
      file >> ICYCLE ;
      if( file ) // Not end of file !
      {
         check_file( file ) ;
         file >> NBVARS ;
      }
   }
   else if( file_format == CBinary )
   {
      check_file( file ) ;
      file.read( (char*)&ICYCLE, sizeof(int) ) ;
      if( file ) // Not end of file !
      {
         file.read( (char*)&NBVARS, sizeof(int) ) ;
      }
   }
   else if( file_format == CBinary ||
            file_format == Binary ||
            file_format == Binary_no_local )
   {
      if( file_format == CBinary )
      {
         check_file( file ) ;
         file.read( (char*)&ICYCLE, sizeof(int) ) ;
      }
      else
      {
         int length = read_length( file, file_format, -1 ) ; // Check end of file.
         if( file ) // Not end of file !
         {
            if( length != 2*sizeof(int) )
            {
               PEL_TICio_ERROR:: n1() ;
            }
            check_file( file ) ;
            file.read( (char*)&ICYCLE, sizeof(int) ) ;
            convert_to_local( file_format, &ICYCLE ) ;
         }
      }
      if( file ) // Not end of file !
      {
         check_file( file ) ;
         file.read( (char*)&NBVARS, sizeof(int) ) ;
         convert_to_local( file_format, &NBVARS ) ;
         read_length( file, file_format, 2*sizeof(int) ) ;
      }
   }
   else
   {
      PEL_TICio_ERROR:: n1() ;
   }

   nb_vars = (size_t) NBVARS ;
   i_cycle = (size_t) ICYCLE ;
}

//----------------------------------------------------------------------
void
PEL_TICio:: read_new_variable( std::ifstream& file,  TIC_FORMAT file_format,
                               std::string& name,
                               size_t& type,
                               size_t& dim, size_t& end,
                               bool& is_integer )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: read_new_variable" ) ;
   PEL_CHECK( file.is_open() ) ;
   PEL_CHECK( file_format != Unspecified ) ;

   char CNAME[5] ;
   int IDIM ;
   intVector LEN(2) ;

   if( file_format == Text )
   {
      char c =' ' ;
      while( c!='\'' )
      {
         check_file( file ) ;
         file.get(c) ;
      }
      for( size_t i=0 ; i<5 ; i++ )
      {
         check_file( file ) ;
         file.get(CNAME[i]) ;
      }
      check_file( file ) ;
      file >> IDIM ;
      if( IDIM<0 || IDIM>2 )
      {
         PEL_TICio_ERROR:: n1() ;
      }
      for( size_t i=0 ; i<(size_t) IDIM ; i++ )
      {
         check_file( file ) ;
         file >> LEN(i) ;
      }
   }
   else if( file_format == CBinary ||
            file_format == Binary ||
            file_format == Binary_no_local )
   {
      size_t lenght = read_length( file, file_format, -1 ) ; // Unknown size
      check_file( file ) ;
      file.read( CNAME, 4 ) ;
      check_file( file ) ;
      file.read( (char*)&IDIM, sizeof(int) ) ;
      convert_to_local( file_format, &IDIM ) ;
      
      if( IDIM<0 || IDIM>2 )
      {
         PEL_TICio_ERROR:: n1() ;
      }
      for( size_t i=0 ; i<(size_t) IDIM ; i++ )
      {
         check_file( file ) ;
         file.read( (char*)&LEN(i), sizeof(int) ) ;
         convert_to_local( file_format, &LEN(i) ) ;
      }
      read_length( file, file_format, lenght ) ;
   }
   else
   {
      PEL_TICio_ERROR:: n1() ;
   }
   type = (size_t) IDIM ;
   dim = type ;
   if( dim == 0 )
   {
      end = 1 ;
   }
   else if( dim == 1 )
   {
      end = (size_t) LEN(0) ;
   }
   else
   {
      dim=(size_t) LEN(0) ;
      end=(size_t) LEN(0)*LEN(1) ;
   }
   is_integer = ( CNAME[0]>='i' && CNAME[0]<='n' ) ||
                ( CNAME[0]>='I' && CNAME[0]<='N' ) ;
   
   CNAME[4]=0 ;
   for( size_t c=3 ; c>0 ; c-- )
   {
      if( CNAME[c]==' ' ) CNAME[c] = 0 ;
   }
   name = std::string( CNAME ) ;
}

//----------------------------------------------------------------------
void
PEL_TICio:: read_float( std::ifstream& file,
                        TIC_FORMAT file_format,
                        size_t dim, size_t end,
                        float* XVAR )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: read_float" ) ;
   PEL_CHECK( file.is_open() ) ;
   PEL_CHECK( file_format != Unspecified ) ;
   PEL_CHECK( XVAR != 0 ) ;
   
   double x ;

   if( dim == 0 )
   {
      if( file_format == Text )
      {
         check_file( file ) ;
         file >> x ;
         *XVAR = (float) x ;
      }
      else if( file_format == CBinary ||
               file_format == Binary ||
               file_format == Binary_no_local )
      {
         read_length( file, file_format, sizeof( float )  ) ;
         check_file( file ) ;
         file.read( (char*)XVAR, sizeof( float) ) ;
         read_length( file, file_format, sizeof( float )  ) ;
         convert_to_local( file_format, XVAR ) ;
      }
      else
      {
         PEL_TICio_ERROR:: n1() ;
      }
   }
   else if( dim == 1 )
   {
      if( file_format == Text )
      {
         for( size_t i=0 ; i<end ; ++i )
         {
            check_file( file ) ;
            file >> x ;
            XVAR[i] = (float) x ;
         }
      }
      else if( file_format == CBinary ||
               file_format == Binary ||
               file_format == Binary_no_local )
      {
         read_length( file, file_format, end*sizeof( float ) ) ;
         file.read( (char*)XVAR, end*sizeof( float) ) ;
         read_length( file, file_format, end*sizeof( float )  ) ;
         for( size_t i=0 ; i<end ; ++i )
         {
            convert_to_local( file_format, &XVAR[i] ) ;
         }
      }
      else
      {
         PEL_TICio_ERROR:: n1() ;
      }
   }
   else 
   {
      size_t const nb_lines = dim ;
      size_t const nb_cols = end/dim ;

      if( file_format == Text )
      {
         for( size_t j=0 ; j<nb_cols ; j++ )
         {
            size_t idx = j*nb_lines ;
            for( size_t i=0 ; i<nb_lines ; i++ )
            {
               check_file( file ) ;
               file >> x ;
               XVAR[idx+i] = (float) x ;
            }
         }
      }
      else if( file_format == CBinary ||
               file_format == Binary ||
               file_format == Binary_no_local )
      {
         read_length( file, file_format, end*sizeof( float ) ) ;
         file.read( (char*)XVAR, end*sizeof( float) ) ;
         read_length( file, file_format, end*sizeof( float )  ) ;
      }
      else
      {
         PEL_TICio_ERROR:: n1() ;
      }
      
      for( size_t i=0 ; i<end ; ++i )
      {
         convert_to_local( file_format, &XVAR[i] ) ;
      }
   }
}

//----------------------------------------------------------------------
void
PEL_TICio:: read_int( std::ifstream& file,
                      TIC_FORMAT file_format,
                      size_t dim, size_t end,
                      int* IVAR )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: read_float" ) ;
   PEL_CHECK( file.is_open() ) ;
   PEL_CHECK( file_format != Unspecified ) ;
   PEL_CHECK( IVAR != 0 ) ;
   
   if( dim == 0 )
   {
      if( file_format == Text )
      {
         check_file( file ) ;
         file >> *IVAR ;
      }
      else if( file_format == CBinary ||
               file_format == Binary ||
               file_format == Binary_no_local )
      {
         read_length( file, file_format, sizeof( int )  ) ;
         check_file( file ) ;
         file.read( (char*)IVAR, sizeof( int ) ) ;
         read_length( file, file_format, sizeof( int )  ) ;
         convert_to_local( file_format, IVAR ) ;
      }
      else
      {
         PEL_TICio_ERROR:: n1() ;
      }
   }
   else if( dim == 1 )
   {
      if( file_format == Text )
      {
         for( size_t i=0 ; i<end ; ++i )
         {
            check_file( file ) ;
            file >> IVAR[i] ;
         }
      }
      else if( file_format == CBinary ||
               file_format == Binary ||
               file_format == Binary_no_local )
      {
         read_length( file, file_format, end*sizeof( int ) ) ;
         file.read( (char*)IVAR, end*sizeof( int ) ) ;
         read_length( file, file_format, end*sizeof( int )  ) ;
         for( size_t i=0 ; i<end ; ++i )
         {
            convert_to_local( file_format, &IVAR[i] ) ;
         }
      }
      else
      {
         PEL_TICio_ERROR:: n1() ;
      }
   }
   else 
   {
      size_t const nb_lines = dim ;
      size_t const nb_cols = end/dim ;
      size_t idx = 0 ;
      if( file_format == Text )
      {
         for( size_t i=0 ; i<nb_lines ; i++ )
         {
            for( size_t j=0 ; j<nb_cols ; j++ )
            {
               check_file( file ) ;
               file >> IVAR[idx++] ;
            }
         }
      }
      else if( file_format == CBinary ||
               file_format == Binary ||
               file_format == Binary_no_local )
      {
         read_length( file, file_format, end*sizeof( int ) ) ;
         file.read( (char*)IVAR, end*sizeof( int ) ) ;
         read_length( file, file_format, end*sizeof( int )  ) ;
      }
      else
      {
         PEL_TICio_ERROR:: n1() ;
      }
      for( size_t i=0 ; i<end ; ++i )
      {
         convert_to_local( file_format, &IVAR[i] ) ;
      }
   }
}

//----------------------------------------------------------------------
int
PEL_TICio:: read_length( std::ifstream &file, TIC_FORMAT file_format,
                         int expected_length )  
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_TICio:: read_length" ) ;
   PEL_CHECK( file.is_open() ) ;
   PEL_CHECK( expected_length == -1 || expected_length > 0 ) ;

   int result = expected_length ;

   if( file_format == Binary || file_format == Binary_no_local )
   {
      check_file( file ) ;
      file.read( ( char*) &result, sizeof( int ) ) ;
      if( file_format == Binary_no_local )
         convert_to_local( file_format, &result ) ;
      if( expected_length != -1 && result != expected_length )
      {
         PEL_TICio_ERROR:: n1() ;
      }
   }

   PEL_CHECK_POST( IMPLIES( expected_length != -1,
                            result == expected_length )  ) ;
   return( result ) ;
}

//internal--------------------------------------------------------------
void 
PEL_TICio_ERROR:: n0( void )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
         "*** PEL_TICio error :\n"
         "    error writing TIC file" ) ;
}

//internal--------------------------------------------------------------
void 
PEL_TICio_ERROR:: n1( void )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
         "*** PEL_TICio error :\n"
         "    error reading TIC file" ) ;
}

//internal--------------------------------------------------------------
void 
PEL_TICio_ERROR:: n2( std::string const& file_name )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
         "*** PEL_TICio error :\n"
         "    Unable to open : \""+file_name+"\"" ) ;
}

//internal--------------------------------------------------------------
void 
PEL_TICio_ERROR:: n3( std::string const& file_name )
//internal--------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
         "*** PEL_TICio error :\n"
         "    File already open : \""+file_name+"\"" ) ;
}

//internal--------------------------------------------------------------
void 
PEL_TICio_ERROR:: n4( std::string const& file_name, size_t i_cycle )
//internal--------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << "*** PEL_TICio error:" << std::endl
       << "    Reading TIC file: \""+file_name+"\"" << std::endl
       << "    Unable to restore cycle number: " << i_cycle ;
   PEL_Error::object()->raise_plain( msg.str() ) ;
}
