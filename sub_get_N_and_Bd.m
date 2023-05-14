function [ N , Bd ]  =  sub_get_N_and_Bd ( N0 , dN_dxi , x_node_local )


%% N matrix
N = zeros ( 2 , 8 ) ;
N ( 1 , 1:2:7 )  =  N0;
N ( 2 , 2:2:8 )  =  N0;


%% Bd matrix
% push forward
F = dN_dxi * x_node_local ; % based on iso-parametric mapping
dN_dx = F \ dN_dxi ;  % chain rule

% B matrix
Bd = zeros ( 4 , 8 ) ;
Bd ( 1 , 1:2:7 )  =  2/3*dN_dx ( 1 , : ) ;  % dN/dx
Bd ( 1 , 2:2:8 )  =  -1/3*dN_dx ( 2 , : ) ;  % dN/dy
Bd ( 2 , 1:2:7 )  =  -1/3*dN_dx ( 1 , : ) ;  % dN/dx
Bd ( 2 , 2:2:8 )  =  2/3*dN_dx ( 2 , : ) ;  % dN/dy
Bd ( 3 , 1:2:7 )  =  -1/3*dN_dx ( 1 , : ) ;  % dN/dx
Bd ( 3 , 2:2:8 )  =  -1/3*dN_dx ( 2 , : ) ;  % dN/dy
Bd ( 4 , 1:2:7 )  =  dN_dx ( 2 , : ) ;  % dN/dy
Bd ( 4 , 2:2:8 )  =  dN_dx ( 1 , : ) ;  % dN/dx




end