#include <PDE_FacesOfCellFE.hh>

#include <PEL_VectorIterator.hh>
#include <PEL_assertions.hh>

#include <PDE_FaceFE.hh>

#include <iostream>

using std::cout ; using std::endl ;

//----------------------------------------------------------------------
PDE_FacesOfCellFE:: PDE_FacesOfCellFE( PEL_Object* a_owner,
                                       PEL_Vector const* faces )
//----------------------------------------------------------------------
   : PEL_Object( a_owner )
   , ALL_FACES()
   , IT()
   , OK( false )
{
   PEL_LABEL( "PDE_FacesOfCellFE:: PDE_FacesOfCellFE" ) ;
   
   PEL_VectorIterator* it = PEL_VectorIterator::create( 0, faces ) ;
   it->start() ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      PEL_CHECK( dynamic_cast< PDE_FaceFE* >( it->item() ) != 0 ) ;
      PDE_FaceFE* face = static_cast< PDE_FaceFE* >( it->item() ) ;
      extend_faces( ALL_FACES, face ) ;
   }
   it->destroy() ;
}

//----------------------------------------------------------------------
PDE_FacesOfCellFE:: ~PDE_FacesOfCellFE( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
PDE_FacesOfCellFE:: start( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FacesOfCellFE:: start" ) ;
   
   IT = ALL_FACES.begin() ;
   OK = ( IT != ALL_FACES.end() ) ;
}

//----------------------------------------------------------------------
bool
PDE_FacesOfCellFE:: is_valid( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FacesOfCellFE:: is_valid" ) ;
   
   return( OK ) ;
}

//----------------------------------------------------------------------
void
PDE_FacesOfCellFE:: go_next( void )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FacesOfCellFE:: go_next" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   
   ++IT ;
   OK = ( IT != ALL_FACES.end() ) ; 
}

//----------------------------------------------------------------------
PDE_FaceFE*
PDE_FacesOfCellFE:: item( void ) const
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FacesOfCellFE:: item" ) ;
   PEL_CHECK_PRE( is_valid() ) ;
   
   PDE_FaceFE* result = *IT ;
   
   PEL_CHECK_POST( result != 0 ) ;
   return( result ) ; 
}

//----------------------------------------------------------------------
void
PDE_FacesOfCellFE:: extend_faces( std::vector< PDE_FaceFE* >& vec_of_faces,
                                  PDE_FaceFE* a_face )
//----------------------------------------------------------------------
{
   PEL_LABEL( "PDE_FacesOfCellFE:: extend_faces" ) ;
   
   vec_of_faces.push_back( a_face ) ;
   if( a_face->nb_childs() != 0 )
   {
      for( size_t i=0 ; i<a_face->nb_childs() ; ++i )
      {
         extend_faces( vec_of_faces, a_face->child( i ) ) ;
      }
   }
}

