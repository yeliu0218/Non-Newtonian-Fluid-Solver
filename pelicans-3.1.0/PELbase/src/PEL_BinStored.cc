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

#include <PEL_BinStored.hh>

#include <PEL_assertions.hh>
#include <PEL_Double.hh>
#include <PEL_BoolVector.hh>
#include <PEL_DoubleVector.hh>
#include <PEL_DoubleArray2D.hh>
#include <PEL_DoubleArray3D.hh>
#include <PEL_IntVector.hh>
#include <PEL_IntArray2D.hh>
#include <PEL_IntArray3D.hh>
#include <PEL_Int.hh>
#include <PEL_Error.hh>
#include <PEL_List.hh>
#include <PEL_Module.hh>
#include <PEL_Vector.hh>
#include <PEL_Root.hh>
#include <PEL_Sequence.hh>
#include <PEL_String.hh>
#include <PEL_System.hh>

#include <boolVector.hh>
#include <doubleVector.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>
#include <intArray2D.hh>
#include <intArray3D.hh>
#include <intVector.hh>

#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>

//internal---------------------------------------------------------------
class PEL_ThisFileDirExp : public PEL_Expression
{
   public:
      PEL_ThisFileDirExp( PEL_Object* a_owner ) ;
     ~PEL_ThisFileDirExp( void ) ; 
      virtual Type data_type( void ) const ;
      virtual std::string const& to_string( PEL_Context const* ct = 0 ) const ;
      virtual PEL_ThisFileDirExp* create_replica( 
                                   PEL_Object * a_owner,
                                   PEL_Sequence const* argument_list ) const ;
      
      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( PEL_Sequence const* some_arguments ) const ;
} ;
//internal---------------------------------------------------------------

//----------------------------------------------------------------------
PEL_BinStored const* PEL_BinStored::static_OpComponent = new PEL_BinStored() ;
const size_t PEL_BinStored::bad_record = (size_t)~0 ;
size_t PEL_BinStored::last_record = bad_record ;
std::string PEL_BinStored::last_file = "" ;
struct magic 
{     char magic_char[256]  ;
      int magic_int ;
      size_t magic_size_t ;
      double magic_double ;
      bool magic_bool ;
} ;
static struct magic magic_stamp = { "PELICANS binary output format rev 1.1",
                                    -2,
                                    (size_t)~0,
                                    1./3.0,
                                    true } ;

size_t PEL_BinStored::preamble_length = sizeof( magic ) ;
//----------------------------------------------------------------------

//----------------------------------------------------------------------
PEL_BinStored:: PEL_BinStored( void ) 
//----------------------------------------------------------------------
      : PEL_TransferExp( "binary" ),
        my_data( 0 )
{
   PEL_LABEL( "PEL_BinStored:: PEL_BinStored" ) ;
}

//----------------------------------------------------------------------
PEL_BinStored:: PEL_BinStored( PEL_Object* a_owner,
                             PEL_Sequence const* argument_list ) 
//----------------------------------------------------------------------
      : PEL_TransferExp( a_owner, "binary", argument_list ),
        my_data( 0 )
{
   PEL_LABEL( "PEL_BinStored:: PEL_BinStored" ) ;
   PEL_CHECK_PRE( argument_list!=0 ) ;
   PEL_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
PEL_BinStored:: ~PEL_BinStored( void ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BinStored:: ~PEL_BinStored" ) ;
   PEL_CHECK_INV( invariant() ) ;
   if( my_data!=0 )  my_data->destroy() ;
}

//----------------------------------------------------------------------
PEL_BinStored*
PEL_BinStored:: create_replica( PEL_Object* a_owner,
                             PEL_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BinStored:: create_replica" ) ;
   PEL_CHECK_PRE( create_replica_PRE( a_owner, argument_list ) ) ;

   PEL_BinStored* result = new PEL_BinStored( a_owner, argument_list ) ;

   PEL_CHECK_POST( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
PEL_BinStored:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result = "binary(SS,SS,IS)" ;
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_BinStored:: valid_arguments( PEL_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BinStored:: valid_arguments" ) ;
   PEL_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = some_arguments->count()==3 &&
      extract_arg(some_arguments,0)->data_type() == String &&
      extract_arg(some_arguments,1)->data_type() == String &&
      extract_arg(some_arguments,2)->data_type() == Int ;
   return result ;
}

//----------------------------------------------------------------------
bool
PEL_BinStored:: is_type_supported( PEL_Data::Type a_type )
//----------------------------------------------------------------------
{
   bool result =
      a_type==Double ||
      a_type==DoubleVector ||
      a_type==DoubleArray2D ||
      a_type==DoubleArray3D ||
      a_type==IntArray2D ||
      a_type==IntArray3D ||
      a_type==IntVector ||
      a_type==BoolVector ;
   PEL_CHECK_POST( EQUIVALENT( result, ( a_type==Double ||
                                         a_type==DoubleVector ||
                                         a_type==DoubleArray2D ||
                                         a_type==DoubleArray3D ||
                                         a_type==IntArray2D ||
                                         a_type==IntArray3D ||
                                         a_type==IntVector ||
                                         a_type==BoolVector ) ) ) ;
   return result ;
}

//----------------------------------------------------------------------
PEL_Data::Type
PEL_BinStored:: data_type( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BinStored:: kind" ) ;
   std::string tmp = arg(0)->to_string() ;
   Type result = Undefined ;
   
   if( tmp=="Double" )
   {
      result = Double ;
   }
   else if( tmp=="DoubleVector" )
   {
      result = DoubleVector ;
   }
   else if( tmp=="DoubleArray2D" )
   {
      result = DoubleArray2D ;
   }
   else if( tmp=="DoubleArray3D" )
   {
      result = DoubleArray3D ;
   }
   else if( tmp=="IntArray2D" )
   {
      result = IntArray2D ;
   }
   else if( tmp=="IntArray3D" )
   {
      result = IntArray3D ;
   }
   else if( tmp=="IntVector" )
   {
      result = IntVector ;
   }
   else if( tmp=="BoolVector" )
   {
      result = BoolVector ;
   }
   else
   {
      PEL_Error::object()->raise_plain( "Unsupported data type " + tmp ) ;
   }
   
   return result ;
}

//----------------------------------------------------------------------
PEL_Data const*
PEL_BinStored:: data( PEL_Context const* ct ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BinStored:: data" ) ;
   PEL_CHECK( data_PRE( ct ) ) ;
   
   if( my_data==0 )
   {
      std::string const& in_file = arg(1)->to_string( ct ) ;
      if( !is_valid_binary_file( in_file ) )
      {
         PEL_Error::object()->raise_plain(
            " File " + in_file
            + " is not a valid binary file on this system " ) ;
      }
      
      size_t record_number = arg(2)->to_int( ct ) ;
      my_data = restore_from_binary( 0, in_file, record_number ) ;
      
      PEL_ASSERT( my_data->data_type() == data_type() ) ;
   }

   PEL_CHECK_POST( data_POST( my_data, ct ) ) ;
   return my_data ;
}

//----------------------------------------------------------------------
bool 
PEL_BinStored:: is_valid_binary_file( std::string const& file_name ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BinStored:: is_valid_binary_file" ) ;
   
   std::ifstream stream( file_name.c_str(), std::ios::binary | std::ios::in ) ;
   // Check for no-empty file
   bool result = stream && !stream.eof() ;
   // Check for magic stamp
   if( result )
   {
      struct magic tmp ;
      
      stream.read( (char*)&tmp, sizeof(magic) ) ;
      result = result &&
         ( stream.gcount()==sizeof(magic) ) ;
      std::string str = tmp.magic_char ;
      result = result &&
         str == magic_stamp.magic_char &&
         tmp.magic_int == magic_stamp.magic_int &&
         tmp.magic_size_t == magic_stamp.magic_size_t &&
         tmp.magic_double == magic_stamp.magic_double &&
         tmp.magic_bool == magic_stamp.magic_bool ;
      stream.close() ;
   }
   return result ;
}

//----------------------------------------------------------------------
size_t 
PEL_BinStored:: last_record_number( std::string const& file_name ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BinStored:: last_record_number" ) ;
   PEL_CHECK( is_valid_binary_file( file_name ) ) ;
   
   size_t result = bad_record ;
   // CACHING
   if( file_name==last_file )
   {
      result = last_record ;
   }
   else
   {
      std::ifstream stream( file_name.c_str(), std::ios::binary | std::ios::in ) ;
      // Ignore preamble
      stream.ignore( preamble_length ) ;

      BinaryRecord record ;
      size_t record_size = sizeof( BinaryRecord ) ;
      while( !stream.eof() )
      {
         stream.read( (char*)&record, record_size ) ;
      
         if( (size_t)stream.gcount()!=record_size )
         {
            break ;
         }
         result = record.number ;
         stream.ignore( record.length ) ;
      }
      stream.close() ;
   }
   return result ;
}

//----------------------------------------------------------------------
PEL_Data const*
PEL_BinStored:: restore_from_binary( PEL_Object * a_owner,
                                     std::string const& file_name,
                                     size_t record_number ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BinStored:: restore_from_binary" ) ;
   PEL_CHECK_PRE( record_number<=last_record_number( file_name ) ) ;

   static std::string oldname ;
   static size_t old_record = 0-1 ;
   static std::ifstream stream ;
   // Garder une reference sur un fichier externe est certe plus rapide mais peut generer
   // des erreurs innattendues si le fichier disparait en cours de route.
   static bool keep_reference = false ;
   
   static bool first = true ;
   static off_t last_modif ;
   static struct stat status ;
   if( stat(file_name.c_str(),&status)!=0 )
   {
      PEL_Error::object()->raise_plain(
         "Unable to inquire for binary file "+file_name ) ;
   }
   off_t last = status.st_size ;
   PEL_Data const* result = 0 ;
   if( !keep_reference ||
       ( first || record_number<=old_record || oldname!=file_name || last!=last_modif ) )
   {
      if( !first && keep_reference ) stream.close() ;
      stream.open( file_name.c_str(), std::ios::binary | std::ios::in ) ;
      if( !stream )
      {
         PEL_Error::object()->raise_plain(
            "Unable to open binary file "+file_name ) ;
      }
      PEL_CHECK( stream.good() ) ;
      stream.ignore( preamble_length ) ;
      PEL_CHECK( stream.good() ) ;
      first = false ;
   }
   oldname = file_name ;
   old_record = record_number ;
   last_modif = last ;
   
   BinaryRecord record ;
   size_t record_size = sizeof( BinaryRecord ) ;
   bool found = false ;
   while( !stream.eof() )
   {
      stream.read( (char*)&record, record_size ) ;
      if( (size_t)stream.gcount()!=record_size )
      {
	 PEL_Error::object()->raise_plain( 
	    "Failed to restore corrupted binary file "+file_name ) ;
      }
      if( record.number==record_number )
      {
         found = true ;
         break ;
      }
      PEL_CHECK( record.length>0 ) ;
      stream.ignore( record.length ) ;
      PEL_CHECK( stream.good() ) ;
   }
   if( !found )
   {
      PEL_Error::object()->raise_plain( "Unable to find record" ) ;
   }
   if( record.type==Double )
   {
      PEL_ASSERT( record.length==sizeof( double ) ) ;
      double val ;
      stream.read( (char*)&val, sizeof( double ) ) ;
      result = PEL_Double::create( a_owner, val ) ;
   }
   else if( record.type==DoubleVector )
   {
      size_t nb = record.length/sizeof( double ) ;
      double *dval = new double[ nb ] ;
      stream.read( (char*)dval, record.length ) ;
      doubleVector dv( nb ) ;
      for( size_t i=0 ; i<nb ; i++ )
      {
         dv(i) = dval[i] ;
      }
      result = PEL_DoubleVector::create( a_owner, dv ) ;
      delete [] dval ;
   }
   else if( record.type==DoubleArray2D )
   {
      size_t nb = record.length/sizeof( double ) ;
      double *dval = new double[ nb ] ;
      stream.read( (char*)dval, record.length ) ;
      size_t dim0 = record.foo[0] ;
      size_t dim1 = record.foo[1] ; 
      PEL_ASSERT( nb==dim0*dim1 ) ;
      doubleArray2D da( dim0, dim1 ) ;
      for( size_t i=0 ; i<dim0 ; i++ )
      {
         for( size_t j=0 ; j<dim1 ; j++ )
         {
            da(i,j) = dval[i*dim1+j] ;
         }
      }
      result = PEL_DoubleArray2D::create( a_owner, da ) ;
      delete [] dval ;
   }
   else if( record.type==DoubleArray3D )
   {

      size_t dim0 = record.foo[0] ;
      size_t dim1 = record.foo[1] ;
      size_t dim2 = record.foo[2] ;
      size_t nb = record.length/sizeof( double ) ;
      PEL_ASSERT( nb==dim0*dim1*dim2 ) ;
      
      double *dval = new double[ nb ] ;
      stream.read( (char*)dval, record.length ) ;
      doubleArray3D da( dim0,dim1,dim2 ) ;
      size_t cpt = 0 ;
      for( size_t i=0 ; i<dim0 ; i++ )
      {
         for( size_t j=0 ; j<dim1 ; j++ )
         {
            for( size_t k=0 ; k<dim2 ; k++ )
            {
               da(i,j,k) = dval[cpt++] ;
            }
         }
      }
      result = PEL_DoubleArray3D::create( a_owner, da ) ;
      delete [] dval ;
   }
   else if( record.type==IntArray2D )
   {
      size_t nb = record.length/sizeof( int ) ;
      int *ival = new int[ nb ] ;
      stream.read( (char*)ival, record.length ) ;
      size_t dim0 = record.foo[0] ;
      size_t dim1 = record.foo[1] ;
      PEL_ASSERT( dim0*dim1==nb ) ;
      intArray2D ia( dim0, dim1 ) ;
      for( size_t i=0 ; i<dim0 ; i++ )
      {
         for( size_t j=0 ; j<dim1 ; j++ )
         {
            ia(i,j) = ival[i*dim1+j] ;
         }
      }
      result = PEL_IntArray2D::create( a_owner, ia ) ;
      delete [] ival ;
   }
   else if( record.type==IntArray3D )
   {

      size_t dim0 = record.foo[0] ;
      size_t dim1 = record.foo[1] ;
      size_t dim2 = record.foo[2] ;
      size_t nb = record.length/sizeof( int ) ;
      PEL_ASSERT( nb==dim0*dim1*dim2 ) ;
      
      int *ival = new int[ nb ] ;
      stream.read( (char*)ival, record.length ) ;
      intArray3D ia( dim0,dim1,dim2 ) ;
      size_t cpt = 0 ;
      for( size_t i=0 ; i<dim0 ; i++ )
      {
         for( size_t j=0 ; j<dim1 ; j++ )
         {
            for( size_t k=0 ; k<dim2 ; k++ )
            {
               ia(i,j,k) = ival[cpt++] ;
            }
         }
      }
      result = PEL_IntArray3D::create( a_owner, ia ) ;
      delete [] ival ;
   }
   else if( record.type==IntVector )
   {
      size_t nb = record.length/sizeof( int ) ;
      int *ival = new int[ nb ] ;
      stream.read( (char*)ival, record.length ) ;
      intVector iv( nb ) ;
      for( size_t i=0 ; i<nb ; i++ )
      {
         iv(i) = ival[i] ;
      }
      result = PEL_IntVector::create( a_owner, iv ) ;
      delete [] ival ;
   }
   else if( record.type==BoolVector )
   {
      size_t nb = record.length/sizeof( bool ) ;
      bool *bval = new bool[ nb ] ;
      stream.read( (char*)bval, record.length ) ;
      boolVector bv( nb ) ;
      for( size_t i=0 ; i<nb ; i++ )
      {
         bv(i) = bval[i] ;
      }
      result = PEL_BoolVector::create( a_owner, bv ) ;
      delete [] bval ;
   }
   if( !keep_reference ) stream.close() ;
   PEL_CHECK_POST( result!=0 ) ;
   return result ;
}

//----------------------------------------------------------------------
void
PEL_BinStored:: init_binary_file( std::string const& file_name )
//----------------------------------------------------------------------
{
   std::ofstream stream( file_name.c_str(), std::ios::binary | std::ios::out ) ;
   if( !stream )
   {
      PEL_Error::object()->raise_plain( "Unable to open "+file_name ) ;
   }
   stream.write( (char*)&magic_stamp, sizeof(magic) ) ;
   stream.close() ;
   PEL_CHECK_POST( is_valid_binary_file( file_name ) ) ;
}

//----------------------------------------------------------------------
PEL_BinStored const*
PEL_BinStored:: create_reference( PEL_Object* a_owner,
                                  PEL_Data const* a_data,
                                  std::string const& file_name,
                                  bool local_reference_file_name ) 
//----------------------------------------------------------------------
{
   PEL_LABEL( "PEL_BinStored:: create_reference" ) ;
   PEL_CHECK_PRE( a_data!=0 ) ;
   PEL_CHECK_PRE( a_data->value_can_be_evaluated(0) ) ;
   PEL_CHECK_PRE( is_type_supported( a_data->data_type() ) ) ;
   PEL_CHECK_PRE( !file_name.empty() ) ;
   PEL_CHECK_PRE( is_valid_binary_file( file_name ) ) ;

   size_t record_number = last_record_number( file_name ) ;
   if( record_number==bad_record )
   {
      record_number=0 ;
   }
   else
   {
      record_number++ ;
   }
   
   std::ofstream stream( file_name.c_str(), std::ios::binary
                         | std::ios::out | std::ios::app ) ;
   
   if( !stream )
   {
      PEL_Error::object()->raise_plain( "Unable to open "+file_name ) ;
   }
   BinaryRecord record ;
   record.type = a_data->data_type() ;
   record.number = record_number ;
   size_t nb = 0 ;
   if( record.type==Double )
   {
      nb = 1 ;
      record.length=sizeof( double ) ;
   }
   else if( record.type==DoubleVector )
   {
      nb = a_data->to_double_vector().size() ;
      record.length=nb*sizeof( double ) ;
   }
   else if( record.type==DoubleArray2D )
   {
      doubleArray2D const& ref = a_data->to_double_array2D() ;
      nb = ref.index_bound(0)*ref.index_bound(1) ;
      record.length=nb*sizeof( double ) ;
      record.foo[0] = ref.index_bound(0) ;
      record.foo[1] = ref.index_bound(1) ;
   }
   else if( record.type==DoubleArray3D )
   {
      doubleArray3D const& refd3D = a_data->to_double_array3D() ;
      nb = refd3D.index_bound(0)*
         refd3D.index_bound(1)*refd3D.index_bound(2) ;
      record.length=nb*sizeof( double ) ;
      record.foo[0] = refd3D.index_bound(0) ;
      record.foo[1] = refd3D.index_bound(1) ;
      record.foo[2] = refd3D.index_bound(2) ;
   }
   else if( record.type==IntArray2D )
   {
      intArray2D const& refi2D = a_data->to_int_array2D() ;
      nb = refi2D.index_bound(0)*
         refi2D.index_bound(1) ;
      record.length=nb*sizeof( int ) ;
      record.foo[0] = refi2D.index_bound(0) ;
      record.foo[1] = refi2D.index_bound(1) ;
   }
   else if( record.type==IntArray3D )
   {
      intArray3D const& refi3D = a_data->to_int_array3D() ;
      nb = refi3D.index_bound(0)*
         refi3D.index_bound(1)*refi3D.index_bound(2) ;
      record.length=nb*sizeof( int ) ;
      record.foo[0] = refi3D.index_bound(0) ;
      record.foo[1] = refi3D.index_bound(1) ;
      record.foo[2] = refi3D.index_bound(2) ;
   }
   else if( record.type==IntVector )
   {
      nb = a_data->to_int_vector().size() ;
      record.length=nb*sizeof( int ) ;
   }
   else if( record.type==BoolVector )
   {
      nb = a_data->to_bool_vector().size() ;
      record.length=nb*sizeof( bool ) ;
   }
   else
   {
      PEL_Error::object()->raise_plain( "Internal error" ) ;
   }
   
   stream.write( (char*)&record, sizeof( BinaryRecord ) ) ;
   if( record.type==Double )
   {
      double val = a_data->to_double() ;
      stream.write( (char*)&val, record.length ) ;
   }
   else if( record.type==DoubleVector )
   {
      double *dval = new double[ nb ] ;
      doubleVector const& dv = a_data->to_double_vector() ;
      for( size_t i=0 ; i<nb ; i++ )
      {
         dval[i] = dv(i) ;
      }
      stream.write( (char*)dval, record.length ) ;
      delete [] dval ;
   }
   else if( record.type==DoubleArray2D )
   {
      double *dval = new double[ nb ] ;
      doubleArray2D const& da = a_data->to_double_array2D() ;
      size_t nb0 = da.index_bound(0) ;
      size_t nb1 = da.index_bound(1) ;
      
      for( size_t i=0 ; i<nb0 ; i++ )
      {
         for( size_t j=0 ; j<nb1 ; j++ )
         {
            dval[i*nb1+j] = da(i,j) ;
         }
      }
      stream.write( (char*)dval, record.length ) ;
      delete [] dval ;
   }
   else if( record.type==DoubleArray3D )
   {
      double *dval = new double[ nb ] ;
      doubleArray3D const& da = a_data->to_double_array3D() ;
      size_t nb0 = da.index_bound(0) ;
      size_t nb1 = da.index_bound(1) ;
      size_t nb2 = da.index_bound(2) ;
      size_t cpt = 0 ;
      for( size_t i=0 ; i<nb0 ; i++ )
      {
         for( size_t j=0 ; j<nb1 ; j++ )
         {
            for( size_t k=0 ; k<nb2 ; k++ )
            {
               dval[cpt++] = da(i,j,k) ;
            }
         }
      }
      stream.write( (char*)dval, record.length ) ;
      delete [] dval ;
   }
   else if( record.type==IntArray3D )
   {
      int *ival = new int[ nb ] ;
      intArray3D const& ia = a_data->to_int_array3D() ;
      size_t nb0 = ia.index_bound(0) ;
      size_t nb1 = ia.index_bound(1) ;
      size_t nb2 = ia.index_bound(2) ;
      size_t cpt = 0 ;
      for( size_t i=0 ; i<nb0 ; i++ )
      {
         for( size_t j=0 ; j<nb1 ; j++ )
         {
            for( size_t k=0 ; k<nb2 ; k++ )
            {
               ival[cpt++] = ia(i,j,k) ;
            }
         }
      }
      stream.write( (char*)ival, record.length ) ;
      delete [] ival ;
   }
   else if( record.type==IntArray2D )
   {
      int *ival = new int[ nb ] ;
      intArray2D const& ia = a_data->to_int_array2D() ;
      size_t nb0 = ia.index_bound(0) ;
      size_t nb1 = ia.index_bound(1) ;
      
      for( size_t i=0 ; i<nb0 ; i++ )
      {
         for( size_t j=0 ; j<nb1 ; j++ )
         {
            ival[i*nb1+j] = ia(i,j) ;
         }
      }
      stream.write( (char*)ival, record.length ) ;
      delete [] ival ;
   }
   else if( record.type==IntVector )
   {
      int *ival = new int[ nb ] ;
      intVector const& iv = a_data->to_int_vector() ;
      for( size_t i=0 ; i<nb ; i++ )
      {
         ival[i] = iv(i) ;
      }
      stream.write( (char*)ival, record.length ) ;
      delete [] ival ;
   }
   else if( record.type==BoolVector )
   {
      bool *bval = new bool[ nb ] ;
      boolVector const& bv = a_data->to_bool_vector() ;
      for( size_t i=0 ; i<nb ; i++ )
      {
         bval[i] = bv(i) ;
      }
      stream.write( (char*)bval, record.length ) ;
      delete [] bval ;
   }
   
   PEL_List* argument_list = PEL_List::create( 0 ) ;
   argument_list->append(
      PEL_String::create( argument_list, type_name( a_data->data_type() ) ) ) ;
   if( local_reference_file_name )
   {
      PEL_List* args = PEL_List::create( 0 ) ;
      args->append( new PEL_ThisFileDirExp( args ) ) ;
      args->append( PEL_String::create( args, PEL_System::basename( file_name ) ) ) ;
      PEL_Data* join = PEL_Expression::create( argument_list, "join", args ) ;
      args->set_owner( join ) ;
      argument_list->append( join ) ;
   }
   else
   {
      argument_list->append(
         PEL_String::create( argument_list, file_name  ) ) ;
   }
   argument_list->append(
      PEL_Int::create( argument_list, record_number  ) ) ;
   
   PEL_BinStored* result = new PEL_BinStored( a_owner, argument_list ) ;
   argument_list->set_owner( result ) ;
   if( !stream.good() )
   {
      PEL_Error::object()->raise_plain( "Unable to write "+file_name ) ;
   }   
   stream.close() ;
   
   // CACHING
   last_file = file_name ;
   last_record = record_number ;
   
   PEL_CHECK_POST( result!=0 ) ;
   PEL_CHECK_POST( result->owner()==a_owner ) ;
   
   return result ;
}

//internal---------------------------------------------------------------
PEL_ThisFileDirExp:: PEL_ThisFileDirExp( PEL_Object* a_owner ) 
//internal---------------------------------------------------------------
      : PEL_Expression( a_owner, "this_file_dir", PEL_List::create( a_owner ) )
{
}

//internal---------------------------------------------------------------
PEL_ThisFileDirExp:: ~PEL_ThisFileDirExp( void ) 
//internal---------------------------------------------------------------
{
}

//internal---------------------------------------------------------------
PEL_ThisFileDirExp*
PEL_ThisFileDirExp:: create_replica(
           PEL_Object* a_owner, PEL_Sequence const* argument_list ) const
//internal---------------------------------------------------------------
{
   PEL_Error::object()->raise_internal( "Can not be called" ) ;
   return( 0 ) ;
}

//internal---------------------------------------------------------------
PEL_Data::Type
PEL_ThisFileDirExp:: data_type( void ) const
//internal---------------------------------------------------------------
{
   return( String ) ;
}

//internal---------------------------------------------------------------
std::string const&
PEL_ThisFileDirExp:: to_string( PEL_Context const* ct ) const
//internal---------------------------------------------------------------
{
   PEL_Error::object()->raise_plain(
      "*** PEL_BinStored error:\n"
      "    Expression \"binary\" can not be evaluated.\n"
      "    Use PEL_BinStored::create_reference with\n"
      "    local_reference_file_name = false" ) ;
   static std::string result ;
   return result ;
}

//internal---------------------------------------------------------------
bool
PEL_ThisFileDirExp:: valid_arguments(
                              PEL_Sequence const* some_arguments ) const
//internal---------------------------------------------------------------
{
   return true ;
}

//internal---------------------------------------------------------------
std::string const&
PEL_ThisFileDirExp:: usage( void ) const
//internal---------------------------------------------------------------
{
   PEL_Error::object()->raise_internal( "Can not be called" ) ;
   static std::string result ;
   return result ;
}
