#ifndef AP_HH
#define AP_HH

class doubleArray2D ;
class doubleArray4D ;

class PDE_LocalEquation ;
class PDE_LocalFEcell ;

class AP
{
   public: //-----------------------------------------------------------------
     
//?????? pourquoi pas des const dPio2dGL etc....
      static void add_C_gradsym_row_gradsym_col( PDE_LocalEquation* leq,
                                                 PDE_LocalFEcell const* fe,
                                                 doubleArray4D& dPio2dGL,
                                                 double coef ) ;

      static void add_S_gradsym_row( PDE_LocalEquation* leq,
                                     PDE_LocalFEcell const* fe,
                                     doubleArray2D& S,
                                     double coef ) ;

      static void add_C_gradGL_row_gradGL_col( PDE_LocalEquation* leq,
                                               PDE_LocalFEcell const* fe,
                                               doubleArray2D const& grad_disp,
                                               doubleArray4D const& dPio2dGL,
                                               double coef ) ;

      static void add_S_graddGL_row_col( PDE_LocalEquation* leq,
                                         PDE_LocalFEcell const* fe,
                                         doubleArray2D const& S,
                                         double coef ) ;

      static void add_S_gradGL_row( PDE_LocalEquation* leq,
                                    PDE_LocalFEcell const* fe,
                                    doubleArray2D const& grad_disp,
                                    doubleArray2D const& S,
                                    double coef ) ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

      AP( void ) ;
     ~AP( void ) ;
      AP( AP const& other ) ;
      AP& operator=( AP const& other ) ;

} ;

#endif
