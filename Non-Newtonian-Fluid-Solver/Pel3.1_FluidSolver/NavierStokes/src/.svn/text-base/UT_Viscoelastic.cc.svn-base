#include <UT_Viscoelastic.hh>

#include <FE_Parameter.hh>
#include <FE_TimeIterator.hh>

#include <LA_SeqVector.hh>
#include <LA_SymmetricMatrix.hh>

#include <PEL.hh>
#include <PEL_Error.hh>
#include <PEL_ModuleExplorer.hh>
#include <PEL_Communicator.hh>
#include <PEL_Exec.hh>

#include <PDE_DiscreteField.hh>
#include <PDE_LocalFEcell.hh>
#include <PDE_LocalEquation.hh>

#include <FE.hh>

#include <doubleVector.hh>
#include <math.h>

/** @brief Elastic model (Phan-Thien and Thanner model)
 * 
 * @param[out] c_minus_I factor multiplied by (conformation matrix minus unit matrix) 
 * @param[in] eigenvals conformation matrix is given bei the eigenvalues (diagonal matrix)
 * @param[in] relax_time relaxation time
 * @param[in] eps material constant related to the extensional viscosity      
 */
/*------------------------------------------------------------*/
void UT_Viscoelastic::elastic_model_PPT(LA_SeqVector* c_minus_I, 
                                        LA_SeqVector* const eigenvals,
	                                    double const relax_time, double const eps)
/*------------------------------------------------------------*/
{
    PEL_LABEL( "UT_Viscoelastic:: elastic_model_PPT" ) ;

    c_minus_I->nullify() ;
    size_t dim = eigenvals->nb_rows();

    double trace_c = 0.0;
    for(size_t i=0; i<dim; i++)
    	trace_c +=  eigenvals->item(i);

    for(size_t i=0; i<dim; i++)
    	c_minus_I->set_item(i, eigenvals->item(i) - 1.);

    //PPT model
  	double factor = -PEL::exp(eps*(trace_c - 3.))/relax_time;

  	c_minus_I->scale(factor) ;
}

/** @brief Compute the norm of the difference between the
      tensors DD and GAMMADOT
      @param cFE finite elements
      @param DD DEVSS
      @param GAMMADOT strain rate tensor field as an independent variable */
/*------------------------------------------------------------*/
double UT_Viscoelastic::compute_norm_DminusGAMMADOT( PDE_LocalFEcell* cFE,
        PDE_DiscreteField const* DD,
        PDE_DiscreteField const* GAMMADOT)
/*------------------------------------------------------------*/
{
    PEL_LABEL( "UT_Viscoelastic:: compute_norm_DminusGAMMADOT" ) ;


    double dd_at_node=0.,gammadot_at_node=0.,norm=0.,diff,norm_collective;
    size_t DcompNumber=DD->nb_components();
    size_t node=0;

    for ( cFE->start() ; cFE->is_valid() ; cFE->go_next() )
    {
        for ( size_t nodeIdx = 0; nodeIdx < cFE->nb_local_nodes( DD ); ++nodeIdx )
        {
            // Get the global node
            node=cFE->global_node( DD, nodeIdx);

            for (size_t i=0; i < DcompNumber; ++i)
            {
                dd_at_node=DD->DOF_value(0,node,i);
                gammadot_at_node=GAMMADOT->DOF_value(0,node,i);
                diff = dd_at_node-gammadot_at_node;
                if (fabs(dd_at_node)>1e-2) diff/=fabs(dd_at_node);
                norm += pow( diff, 2. );
            }
        }
    }

    PEL_Communicator const* pelCOMM=PEL_Exec::communicator();
    norm_collective = pelCOMM->sum( norm );

    norm_collective=sqrt(norm_collective);

    return( norm_collective );

}
