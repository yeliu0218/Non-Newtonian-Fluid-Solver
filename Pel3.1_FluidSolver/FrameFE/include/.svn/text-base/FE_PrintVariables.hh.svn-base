#ifndef FE_PrintVariables_HH
#define FE_PrintVariables_HH

#include <PEL_Object.hh>

#include <doubleArray2D.hh>
#include <doubleVector.hh>

class PEL_Data ;
class PEL_Double ;
class PEL_ModuleExplorer ;

/*
Print variables: Reymolds, Bingham, etc.
 * MODULE FE_PrintVariables 
 *       names=<"Reynolds" "Bn1 (displacing)" "Bn2 (displaced)" > 
 *       values=vector($DS_Re, $DS_Bn1, $DS_Bn2) 
 *    END MODULE FE_PrintVariables
PUBLISHED
*/

class FE_PrintVariables : public PEL_Object
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static FE_PrintVariables* create( PEL_Object* a_owner,
                                      PEL_ModuleExplorer const* exp ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~FE_PrintVariables( void ) ;
      FE_PrintVariables( void ) ;
      FE_PrintVariables( FE_PrintVariables const& other ) ;
      FE_PrintVariables& operator=( FE_PrintVariables const& other ) ;

      FE_PrintVariables( PEL_Object* a_owner,
                       PEL_ModuleExplorer const* exp ) ;


} ;

#endif
