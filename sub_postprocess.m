function [displ, strain , stress, x_eval] = sub_postprocess ( Model , Mesh , d )
x_node = Mesh.x_node;
connectivity = Mesh.connectivity;
num_node = size(x_node,1);
num_element = size(connectivity,2);

%% Elasticity matrix (D matrix)
D  = sub_elasticity_matrix ( Model.material );

%% parametric coordinate to evaluate strain and stress
% You can plot at Gauss points. Then 
%GI = sub_Gauss_point_local;
%xi = GI.xi_2d;   % The code will plot at xi for all the elements.
%xi = [ -1/sqrt(3),-1 ; 1, -1/sqrt(3); 1/sqrt(3), 1; -1, 1/sqrt(3) ]; % 4 points in element

%xi = [-0.5,-1; 0,-1 ; 0.5,-1; 1, -1; 1,-0.5;1, 0; 1,0.5 ; 0, 0]; % 4 points in element
% xi = [-0.8,-1; 0,-1; 1,-1; 1,0; 1,0.8];

%xi = [-1,-1;0,0];
xi = [0,-1;-1,0];

%% evaluation
% local shape functions at local Gauss points
[N_local, dN_dxi_local, dN_deta_local] = sub_shape_local(xi);

% x_eval, displ   :  (num of evaluation points) - by - 2, Each row has the values for each evaluation point.
% strain, stress  :  (num of evaluation points) - by - 3, Each row has the values for each evaluation point.
displ  = []; % The size will be updated in the loop.
strain = [];
stress = [];
x_eval = [];

% loop over element
for idx_elem  =  1 : num_element
    
    global_node_index   =   connectivity ( : , idx_elem );
    x_node_local        =   x_node ( global_node_index , : );
    
    % associated indices of the global matrices
    index_K  =  zeros ( 2*4 , 1 );
    index_K(1:2:end) = 2*global_node_index-1;
    index_K(2:2:end) = 2*global_node_index;
%     index_K  =  sort ( [2*global_node_index-1 ; 2*global_node_index] );
   % loop over Gauss point
    for idx_eval  =  1 : size(xi,1)
        % shape functions
        N0      =   N_local ( idx_eval , : );
        dN_dxi  = [ dN_dxi_local( idx_eval , : ) ;
                    dN_deta_local( idx_eval , : )  ] ;
        % compute N and B matrices
        [ N , B ]  =  sub_get_N_and_B ( N0 , dN_dxi , x_node_local );
        
        % evaluate values
        u           =   N  *  d ( index_K ) ;   % 2-by-1
        epsilon     =   B  *  d ( index_K ) ;   % 3-by-1 in Voigt notation, 3rd component is gamma
        sigma       =   D  *  epsilon;          % 3-by-1 in Voigt notation, 3rd component is the shear stress
        
        % position
        x  =  N0 * x_node_local ;  % iso-parametric mapping, 1-by-2
        
        % store
        displ           =  [ displ ; u'];
        strain          =  [ strain ; epsilon' ];
        stress          =  [ stress ; sigma' ];
        x_eval          =  [ x_eval ; x ];
    end
    
end


%% Plot
% deformed mesh
% plot_mesh ( Mesh.x_node + [d(1:2:end), d(2:2:end)] , Mesh.connectivity ); title('Deformed Mesh, Approx')
% if Model.exact.use == 1
%     u_exact  =   Model.exact.displ ( Mesh.x_node(:,1), Mesh.x_node(:,2) );
%     plot_mesh ( Mesh.x_node + u_exact , Mesh.connectivity ); title('Deformed Mesh, Exact')
% end


% size(stress)
% x_eval

end