function [ K , f ]  =  sub_assembly_FIRI_HC( Model , Mesh , BC ,int_pts,kappa)

x_node = Mesh.x_node;
connectivity = Mesh.connectivity;
num_node    =   size(x_node,1);
num_element =   size(connectivity,2);
nen         =   size(connectivity,1);  % num of nodes per element
h = [1;-1;1;-1];
E       =  3e7;
%% Elasticity matrix (D matrix)
% D: 3-by-3
D  = sub_elasticity_matrix ( Model.material );

%% Gauss point
%int_pts = 1; %2 for full integration, 1 for reduced integration
GI = sub_Gauss_point_local_mod(int_pts);
GI1 = sub_Gauss_point_local_mod(1);
K  =  zeros ( 2*num_node , 2*num_node );
f  =  zeros ( 2*num_node , 1 );
Ks = zeros(8,8);
%% Domain integration
% As the same Gauss quadrature is used for all the elements, the FE shape
% functions in the parametric coordinate is the same. So, we compute them
% only once here and reuse for all the elements.
[ N_local , dN_dxi_local , dN_deta_local ]  =  sub_shape_local ( GI.xi_2d );
[ N_local1 , dN_dxi_local1 , dN_deta_local1 ]  =  sub_shape_local ( GI1.xi_2d );
% loop over elements
for idx_elem  =  1 : num_element
    
    global_node_index   =   connectivity ( : , idx_elem );    % the four nodes for (idx_elem)-th element
    x_node_element      =   x_node ( global_node_index , : );    % corresponding nodal positions
    % associated indices in the global matrices
    index_K             =   zeros ( 2 * nen  ,  1 );
    index_K(1:2:end)    =   2 * global_node_index  -  1  ;
    index_K(2:2:end)    =   2 * global_node_index        ;
    
    % loop over Gauss point
    for idx_Gauss  =  1 : 4
        % shape functions
        N0      =   N_local ( idx_Gauss , : );
        dN_dxi  = [ dN_dxi_local( idx_Gauss , : ) ;
                    dN_deta_local( idx_Gauss , : )  ] ;
                
        % compute N and B matrices
        [ N , B ]  =  sub_get_N_and_B ( N0 , dN_dxi , x_node_element );
        
        % body force
        x_gauss_point  =  N0 * x_node_element ;  % physical coordinate of the current gauss point, based on iso-parametric mapping
        b  =  Model.body_force ( x_gauss_point(1) , x_gauss_point(2) )';
        
        % get Jacobian
        F = dN_dxi * x_node_element ; % based on iso-parametric mapping
        J = det(F);

        % Stabilisation matrix 
        dN_dxi1  = [ dN_dxi_local1( idx_Gauss , : ) ;
                    dN_deta_local1( idx_Gauss , : )  ] ;
        J2 = det(dN_dxi1 * x_node_element);
        alphaxx = 1/8*kappa*J2*E*(sum(dN_dxi1(1,:).^2) + sum(dN_dxi1(2,:).^2));
        alphayy = 1/8*kappa*J2*E*(sum(dN_dxi1(1,:).^2) + sum(dN_dxi1(2,:).^2));
        gamma = (h - (h'*x_node_element(:,1))*dN_dxi1(1,:)' - (h'*x_node_element(:,2))*dN_dxi1(2,:)' );
        Ks(1:2:7,1:2:7) = alphaxx*(gamma*gamma') ;
        Ks(2:2:8,2:2:8) = alphayy*(gamma*gamma') ;

        % Assemble
        K ( index_K , index_K )  =  K ( index_K , index_K )  +  ...
                                    B' * D * B * J * GI.weight_2d(idx_Gauss) + Ks;
        f ( index_K )            =  f ( index_K )  +  ...
                                    N' * b * J * GI.weight_2d(idx_Gauss);
            
    end
    
    
end


%% Natural boundary
% local shape functions at local Gauss points
GI2 = sub_Gauss_point_local_mod(2);
N_local = sub_shape_1d_local ( GI2.xi_1d );
% loop over the NBC elements
num_nbc_segments = length ( BC.nbc_connectivity ) ;
for idx_nbc_seg  =  1 : num_nbc_segments
    num_nbc_element  =  size ( BC.nbc_connectivity{idx_nbc_seg} , 2 );
    for idx_nbc  =  1 : num_nbc_element

        global_node_index   =   BC.nbc_connectivity{idx_nbc_seg} ( : , idx_nbc );
        x_node_element        =   x_node ( global_node_index , : );

        % associated indices of the global matrices
        index_K  =  zeros ( 2*2 , 1 );
        index_K(1:2:end) = 2*global_node_index-1;
        index_K(2:2:end) = 2*global_node_index;

        for idx_Gauss  =  1 : 2
            % shape functions
            N0      =   N_local ( idx_Gauss , : );
            % traction
            x_gauss_point  =  N0 * x_node_element ;  % physical coordinate of the current gauss point, based on iso-parametric mapping
            h  =  Model.bc.nbc_value{idx_nbc_seg} ( x_gauss_point(1) , x_gauss_point(2) )';
            % Jacobian
            J  =  norm(x_node_element(2,:)-x_node_element(1,:)) / 2;
            % N matrix
            N = zeros ( 2 , 4 ) ;
            N ( 1 , 1:2:3 )  =  N0;
            N ( 2 , 2:2:4 )  =  N0;
            % Assemble
            f ( index_K ) =  f ( index_K )  +  ...
                                N' * h * J * GI2.weight_1d(idx_Gauss);

        end
    end
end


%% Make K sparse for fast inversion
K = sparse(K);

end