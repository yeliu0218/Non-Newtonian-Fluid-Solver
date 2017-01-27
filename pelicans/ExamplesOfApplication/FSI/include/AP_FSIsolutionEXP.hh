#ifndef AP_FSI_SOLUTION_EXP_HH
#define AP_FSI_SOLUTION_EXP_HH

#include <PEL_Expression.hh>

#include <doubleVector.hh>
#include <doubleArray2D.hh>

/*
Expressions related to the "Flow in a flexible duct" FSI benchmark
   
PUBLISHED
*/

class AP_FSIsolutionEXP : public PEL_Expression
{
   public: //-------------------------------------------------------
      
   //-- Type
      
      virtual PEL_Data::Type data_type( void ) const ;
      
   //-- Value

      virtual doubleVector const& to_double_vector( 
                                PEL_Context const* ct = 0 ) const ;

      virtual doubleArray2D const& to_double_array2D( 
                                PEL_Context const* ct = 0 ) const ;

   protected: //-----------------------------------------------------
      
   private: //-------------------------------------------------------

      AP_FSIsolutionEXP( void ) ;
     ~AP_FSIsolutionEXP( void ) ;
      AP_FSIsolutionEXP( AP_FSIsolutionEXP const& other ) ;
      AP_FSIsolutionEXP& operator=( AP_FSIsolutionEXP const& other ) ;

      enum Func { UU, PP, MU, DD, SF, SS, UD, SD } ;
            
      AP_FSIsolutionEXP( PEL_Object* a_owner,
                         std::string const& a_name,
                         PEL_Sequence const* argument_list,
                         Func an_expr ) ;
      
   //-- Plug in

      AP_FSIsolutionEXP( std::string const& a_name, Func an_expr ) ;

      virtual AP_FSIsolutionEXP* create_replica( 
                            PEL_Object * a_owner,
                            PEL_Sequence const* argument_list ) const ;

   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( PEL_Sequence const* some_arguments ) const ;

   //-- Class attributes

      static AP_FSIsolutionEXP const* PROTOTYPE_UU ;
      static AP_FSIsolutionEXP const* PROTOTYPE_PP ;
      static AP_FSIsolutionEXP const* PROTOTYPE_MU ;
      static AP_FSIsolutionEXP const* PROTOTYPE_DD ;
      static AP_FSIsolutionEXP const* PROTOTYPE_SF ;
      static AP_FSIsolutionEXP const* PROTOTYPE_SS ;
      static AP_FSIsolutionEXP const* PROTOTYPE_UD ;
      static AP_FSIsolutionEXP const* PROTOTYPE_SD ;

   //-- Attributes

      Func EXPR ;
      mutable doubleVector DV_result_1 ;
      mutable doubleVector DV_result_2 ;
      mutable doubleArray2D doubleArray2D_result ;
} ;

#endif
