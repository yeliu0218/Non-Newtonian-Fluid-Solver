#ifndef IS_FACE_TO_STREAM_CELL_NODES_HH
#define IS_FACE_TO_STREAM_CELL_NODES_HH

#include <PEL_Object.hh>

#include <PDE_DiscreteField.hh>

class FE_TimeIterator ;

class PDE_CursorFEside ;
class PDE_DomainAndFields ;
class PDE_LocalFEbound ;

class PEL_Vector ;

class MY_OneBoundStreamCellNodes ;
class MY_OneSideStreamCellNodes ;

/*
Servers that provide the connectivity from a `PDE_CursorFEside::' or 
a `PDE_LocalFEbound::' to the nodes of the related cells following 
the scheme :

A.) from `PDE_CursorFEside::' :

                               N(0->1)
                                -->  
             |         |    0    |    1    |         |
             |----x----|----x----|----x----|----x----|
             |         |         |         |         |
                  ^         ^    ^    ^         ^
 Up-upstream node_|         |    |    |         |_ Down-downstream node
             Upstream node _|    |    |_ Downstream node 
                                 |_ Current face

The notion of upstream direction is a purely geometrical one, and must be 
understood as the opposite direction to the normal of the `PDE_CursorFEside::'
that is given by `PDE_CursorFEside::normal' and goes from adjacent cell 
number 0 to adjacent cell number 1.

Specific cases:

   A.1.) No up-upstream node :
                                N(0->1)
                                -->  
                     /|    0    |    1    |         |
                     /|----x----|----x----|----x----|
                     /|         |         |         |
                      ^    ^    ^    ^         ^
            Boundary _|    |    |    |         |_ Down-downstream node
            Upstream node _|    |    |_ Downstream node 
                                |_ Current face

   A.2.) No down-downstream node :

                                N(0->1)
                                -->  
             |         |    0    |    1    |/
             |----x----|----x----|----x----|/
             |         |         |         |/
                  ^         ^    ^    ^    ^
 Up-upstream node_|         |    |    |    |_ Boundary
             Upstream node _|    |    |_ Downstream node 
                                |_ Current face

   A.3.) Neither up-upstream node nor down-downstream node :

                                N(0->1)
                                -->  
                     /|    0    |    1    |/
                     /|----x----|----x----|/
                     /|         |         |/        
                      ^    ^    ^    ^    ^
            Boundary _|    |    |    |    |_ Boundary
            Upstream node _|    |    |_ Downstream node 
                                |_ Current face

B.) from `PDE_LocalFEbound::' :

The notion of upstream direction is again a purely geometrical one, and 
must be understood as the opposite direction to the outward normal of 
the `PDE_LocalFEbound::' that is given by 
`PDE_LocalFEbound::outward_normal' and goes from the adjacent cell to
the bound. 


             |   UU    |    U    |/  N(outward normal)
             |----x----|----x----|/  -->
             |         |         |/
                  ^         ^    ^ 
 Up-upstream node_|         |    | 
             Upstream node _|    |___ Current boundary

Specific case:

   Neither up-upstream node nor down-downstream node :

                               
                     /|    U    |/ N(outward normal)
                     /|----x----|/ -->
                     /|         |/        
                           ^    ^
                           |    |_ Current boundary
                           |_ Upstream node 
  
  This first version of those servers needs a structured mesh.
  PUBLISHED  
*/

class MY_MUSCL_DataStructure : public PEL_Object
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      static MY_MUSCL_DataStructure const* object( 
                                            PDE_DomainAndFields const* dom ) ;

   //-- Side to field node connectivity

      // upstream global node from current `PDE_CursorFEside::'
      size_t upstream_node(
                                 PDE_CursorFEside const* fe,
                                 PDE_DiscreteField const* field ) const ;
      
      // distance from the face center to the upstream node
      double upstream_node_dist( PDE_CursorFEside const* fe,
			         PDE_DiscreteField const* field ) const ;

      // upstream global node from current `PDE_LocalFEbound::'
      size_t upstream_node(
                                 PDE_LocalFEbound const* fe,
                                 PDE_DiscreteField const* field ) const ;
      
      // distance from the face center to the upstream node
      double upstream_node_dist( PDE_LocalFEbound const* fe,
			         PDE_DiscreteField const* field ) const ;

      // downstream global node from current `PDE_CursorFEside::'
      size_t downstream_node(
                                 PDE_CursorFEside const* fe,
                                 PDE_DiscreteField const* field ) const ;
      
      // distance from the face center to the downstream node
      double downstream_node_dist(
                                 PDE_CursorFEside const* fe,
                                 PDE_DiscreteField const* field ) const ;

      // Check if the upstream node has an upstream node ?
      bool has_up_upstream_node( PDE_CursorFEside const* fe,
                                 PDE_DiscreteField const* field ) const ;

      // upstream global node of the upstream node 
      size_t up_upstream_node(
                                 PDE_CursorFEside const* fe,
                                 PDE_DiscreteField const* field ) const ;

      // distance from the upstream node to its upstream node
      double up_upstream_node_dist(
                                 PDE_CursorFEside const* fe,
                                 PDE_DiscreteField const* field ) const ;
      
      // Check if the upstream node has an upstream node ?
      bool has_up_upstream_node( PDE_LocalFEbound const* fe,
                                 PDE_DiscreteField const* field ) const ;

      // upstream global node of the upstream node 
      size_t up_upstream_node(
                                 PDE_LocalFEbound const* fe,
                                 PDE_DiscreteField const* field ) const ;

      // distance from the upstream node to its upstream node
      double up_upstream_node_dist(
                                 PDE_LocalFEbound const* fe,
                                 PDE_DiscreteField const* field ) const ;
      
      // Check if the downstream node has a downstream node ?
      bool has_down_downstream_node(
                                 PDE_CursorFEside const* fe,
                                 PDE_DiscreteField const* field ) const ;

      //downstream global node of the downstream node
      size_t down_downstream_node(
                                 PDE_CursorFEside const* fe,
                                 PDE_DiscreteField const* field ) const ;

      // distance from the downstream node to its downstream node
      double down_downstream_node_dist(
                                 PDE_CursorFEside const* fe,
                                 PDE_DiscreteField const* field ) const ;

      //-- Input - Output
      
      void display_connectivities(
                                 PDE_DiscreteField const* field,
                                 std::ostream& os,
                                 size_t indent_width ) const ;

      static size_t 
      side_normal_direction( PDE_CursorFEside const* fe, 
                                                      double eps=1.E-6  ) ;

      static size_t  
      bound_normal_direction( PDE_LocalFEbound const* fe, 
                                                      double eps=1.E-6  ) ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------
      
      MY_MUSCL_DataStructure( PDE_DomainAndFields const* dom ) ;

      MY_MUSCL_DataStructure( void ) ;
     ~MY_MUSCL_DataStructure( void ) ;
      MY_MUSCL_DataStructure( MY_MUSCL_DataStructure const& other ) ;
      MY_MUSCL_DataStructure& operator=(
                                MY_MUSCL_DataStructure const& other ) ;
      
      bool structured_mesh( PDE_DomainAndFields const* dom, double eps ) const ;

      MY_OneSideStreamCellNodes* connectivity(
         PDE_CursorFEside const* fe, PDE_DiscreteField const* field ) const ;
      
      MY_OneBoundStreamCellNodes* connectivity(
         PDE_LocalFEbound const* fe, PDE_DiscreteField const* field ) const ;
      
   //-- Class attributes

      static PEL_Vector* OBJECTS ;

   //-- Attributes

      PDE_DomainAndFields const* const DOM ;
      PEL_Vector* const CONNECTIVITY_SIDES ;
      PEL_Vector* const CONNECTIVITY_BOUNDS ;
} ;

#endif
