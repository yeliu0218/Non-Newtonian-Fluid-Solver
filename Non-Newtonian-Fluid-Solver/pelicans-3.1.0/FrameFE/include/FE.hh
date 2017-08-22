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

#ifndef FE_HH
#define FE_HH

#include <PEL_Object.hh>

class doubleArray2D ;
class doubleVector ;

class PDE_CursorFEside ;
class PDE_DiscreteField ;
class PDE_LocalEquation ;
class PDE_LocalFE ;
class PDE_LocalFEbound ;
class PDE_LocalFEcell ;

/*
PUBLISHED
*/

class PEL_EXPORT FE
{
   public: //--------------------------------------------------------------

   //-- Geometrical characteristics

      enum FE_Geometry
      {
         cartesian,
         axisymmetrical,
         unspecified
      } ;

      static FE_Geometry geometry( void ) ;
      static void set_geometry( FE_Geometry geom ) ;

      static double bound_measure( PDE_LocalFEbound const* fe ) ;
      static double side_measure( PDE_CursorFEside const* fe ) ;
      static double side_measure( PDE_CursorFEside const* fe, size_t i_adj ) ;
      static double cell_measure( PDE_LocalFEcell const* fe ) ;

   //-- Terms of variational approximations : contributions to local systems

      static void add_row_col_S( PDE_LocalEquation* leq,
                                 PDE_LocalFE const* fe, 
                                 double coef ) ;

      static void add_row_col_NS( PDE_LocalEquation* leq,
                                  PDE_LocalFE const* fe, 
                                  double coef ) ;

      static void add_lumped_row_col( PDE_LocalEquation* leq,
                                      PDE_LocalFE const* fe, 
                                      double coef ) ;

      static void add_grad_row_grad_col_S( PDE_LocalEquation* leq,
                                           PDE_LocalFEcell const* fe, 
                                           double coef ) ;

      static void add_grad_row_grad_col_NS( PDE_LocalEquation* leq,
                                            PDE_LocalFEcell const* fe, 
                                            double coef ) ;

      static void add_grad_row_D_col_S( PDE_LocalEquation* leq,
                                        PDE_LocalFEcell const* fe,
                                        double coef ) ;

      static void add_row( PDE_LocalEquation* leq,
                           PDE_LocalFE const* fe,
                           double coef ) ;

      static void add_row( PDE_LocalEquation* leq,
                           PDE_LocalFE const* fe,
                           doubleVector const& coef_1,
                           double coef_2 = 1.0 ) ;

      static void add_grad_row( PDE_LocalEquation* leq,
                                PDE_LocalFE const* fe,
                                doubleVector const& coef ) ;

      static void add_grad_row( PDE_LocalEquation* leq,
                                PDE_LocalFE const* fe,
                                doubleArray2D const& coef ) ;

      static void add_div_row( PDE_LocalEquation* leq,
                               PDE_LocalFE const* fe,
                               double coef ) ;

      static void add_row_vvgrad_col( PDE_LocalEquation* leq,
                                      PDE_LocalFEcell const* fe, 
                                      doubleVector const& aa,
                                      double coef ) ;

      static void add_grad_row_vv_otimes_col( PDE_LocalEquation* leq,
                                              PDE_LocalFEcell const* fe,
                                              doubleVector const& aa,
                                              double coef ) ;

      static void add_vvgrad_row_col( PDE_LocalEquation* leq,
                                      PDE_LocalFEcell const* fe,
                                      doubleVector const& aa,
                                      double coef ) ;

      static void add_vvgrad_row_vvgrad_col( PDE_LocalEquation* leq,
                                             PDE_LocalFEcell const* fe, 
                                             doubleVector const& aa,
                                             double coef ) ;

      static void add_vvgrad_row_lapl_col( PDE_LocalEquation* leq,
                                           PDE_LocalFEcell const* fe,
                                           doubleVector const& aa,
                                           double coef ) ;

      static void add_vvgrad_row_div_D_col( PDE_LocalEquation* leq,
                                            PDE_LocalFEcell const* fe,
                                            doubleVector const& aa,
                                            double coef,
                                            doubleVector const& dcoef ) ;

      static void add_vvgrad_row( PDE_LocalEquation* leq,
                                  PDE_LocalFEcell const* fe,
                                  doubleVector const& aa,
                                  double coef ) ;

      static void add_vvgrad_row( PDE_LocalEquation* leq,
                                  PDE_LocalFEcell const* fe,
                                  doubleVector const& aa,
                                  doubleVector const& coef_1,
                                  double coef_2 = 1.0 ) ;

      static void add_row_div_col( PDE_LocalEquation* leq,
                                   PDE_LocalFEcell const* fe,
                                   double coef ) ;

      static void add_row_grad_col( PDE_LocalEquation* leq,
                                    PDE_LocalFEcell const* fe,
                                    double coef ) ;

      static void add_div_row_div_col( PDE_LocalEquation* leq,
                                       PDE_LocalFEcell const* fe,
                                       double coef ) ;

      static void add_row_vv_col( PDE_LocalEquation* leq,
                                  PDE_LocalFEcell const* fe,
                                  doubleVector const& aa,
                                  double coef ) ;

      static void add_lapl_row_vvgrad_col( PDE_LocalEquation* leq,
                                           PDE_LocalFEcell const* fe,
                                           doubleVector const& aa,
                                           double coef ) ;

      static void add_row_lapl_col( PDE_LocalEquation* leq,
                                    PDE_LocalFEcell const* fe,
                                    double coef ) ;

      static void add_lapl_row_lapl_col( PDE_LocalEquation* leq,
                                         PDE_LocalFEcell const* fe,
                                         double coef ) ;

      static void add_lapl_row_col( PDE_LocalEquation* leq,
                                    PDE_LocalFEcell const* fe,
                                    double coef ) ;

      static void add_lapl_row( PDE_LocalEquation* leq,
                                PDE_LocalFE const* fe,
                                double coef ) ;

      static void add_lapl_row( PDE_LocalEquation* leq,
                                PDE_LocalFE const* fe,
                                doubleVector const& coef_1,
                                double coef_2 = 1.0 ) ;

   //-- Differential operators applied to fields : computation of point values

      // `target' += value at the current integration point of `fe' of
      // `coef'*[grad`u'+(grad`u')^t]:[grad`u'+(grad`u')^t] 
      static void add_D_u_D_u_at_IP( double& target, 
                                     PDE_DiscreteField const* u,
                                     size_t u_level,
                                     PDE_LocalFE const* fe,
                                     double coef ) ;
      
      // `target' += value at the current calculation point of `fe' of
      // `coef'*[grad`u'+(grad`u')^t]:[grad`u'+(grad`u')^t] 
      static void add_D_u_D_u_at_pt( double& target,
                                     PDE_DiscreteField const* u,
                                     size_t u_level,
                                     PDE_LocalFE const* fe,
                                     double coef ) ;

      // `target' += value at the current integration point of `fe' of
      // `coef'* div`u'            
      static void add_div_u_at_IP( double& target,
                                   PDE_DiscreteField const* u,
                                   size_t u_level,
                                   PDE_LocalFE const* fe,
                                   double coef ) ;

      // `target' += value at the current calculation point of `fe' of
      // `coef'*[div`u']              
      static void add_div_u_at_pt( double& target,
                                   PDE_DiscreteField const* u,
                                   size_t u_level,
                                   PDE_LocalFE const* fe,
                                   double coef ) ;

   protected: //-----------------------------------------------------------

   private: //-------------------------------------------------------------

   //-- Static attribute :

      static FE_Geometry GEOM ;

} ;

#endif
