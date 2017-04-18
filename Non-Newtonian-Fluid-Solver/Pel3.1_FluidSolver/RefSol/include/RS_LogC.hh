#ifndef RS_LOG_C_HH
#define RS_LOG_C_HH

#include <PEL_Expression.hh>
#include <doubleVector.hh>

class LA_Matrix ;
class LA_SeqVector ;
class LA_Vector ;

/*

PUBLISHED
*/

class PEL_EXPORT RS_LogC : public PEL_Expression
{
   public: //---------------------------------------------------------------
   
   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual doubleVector const& to_double_vector( 
                                          PEL_Context const* ct = 0 ) const ;

   protected: //-----------------------------------------------------------
      
   private: //-------------------------------------------------------------
    
      RS_LogC( void ) ;
     ~RS_LogC( void ) ;
      RS_LogC( RS_LogC const& other ) ;
      RS_LogC& operator=( RS_LogC const& other ) ;

      enum Func {  Chan, Pipe } ;
            
      RS_LogC( PEL_Object* a_owner,
                   std::string const& a_name,
                   PEL_Sequence const* argument_list,
                   Func an_expr ) ;

   //-- Plug in

      RS_LogC( std::string const& a_name, Func an_expr ) ;

      virtual RS_LogC* create_replica( 
                                   PEL_Object * a_owner,
                                   PEL_Sequence const* argument_list ) const ;
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( PEL_Sequence const* some_arguments ) const ;

   //-- Class attributes
      static RS_LogC const* PROTOTYPE_Chan ;
      static RS_LogC const* PROTOTYPE_Pipe ;
     
   //-- Attributes

      Func EXPR ;

      mutable doubleVector DV_result_3 ;
} ;

#endif
