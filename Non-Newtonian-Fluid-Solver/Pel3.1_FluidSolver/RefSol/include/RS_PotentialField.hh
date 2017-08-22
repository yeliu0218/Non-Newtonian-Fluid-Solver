#ifndef RS_PotentialField_HH
#define RS_PotentialField_HH

#include <PEL_Expression.hh>

/*
    turns electrical field ON and OFF
    via potential field phi. E=-grad phi.

PUBLISHED
*/

class RS_PotentialField : public PEL_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual double to_double( PEL_Context const* ct = 0 ) const ;

   protected: //-----------------------------------------------------
      
   private: //-------------------------------------------------------

      RS_PotentialField( void ) ;
     ~RS_PotentialField( void ) ;
      RS_PotentialField( RS_PotentialField const& other ) ;
      RS_PotentialField& operator=(
                                 RS_PotentialField const& other ) ;
            
      enum Func { ONE, TWO } ;

      RS_PotentialField( PEL_Object* a_owner,
				std::string const& a_name,
				PEL_Sequence const* argument_list,
				Func an_expr ) ;

   //-- Plug in

      RS_PotentialField( std::string const& a_name, Func an_expr ) ;

      virtual RS_PotentialField* create_replica( 
                                   PEL_Object * a_owner,
                                   PEL_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( PEL_Sequence const* some_arguments ) const ;
      
   //-- Class attributes

      static RS_PotentialField const* PROTOTYPE_ONE ;
      static RS_PotentialField const* PROTOTYPE_TWO ;  

   //-- Attributes

      Func const EXPR ;
} ;

#endif
