function [ Bv ]  =  sub_get_Bv ( N0 , dN_dxi , x_node_local )


%% N matrix
% N = zeros ( 2 , 8 ) ;
% N ( 1 , 1:2:7 )  =  N0;
% N ( 2 , 2:2:8 )  =  N0;


%% Bv matrix
% push forward
F = dN_dxi * x_node_local ; % based on iso-parametric mapping
dN_dx = F \ dN_dxi ;  % chain rule

% B matrix
Bv = zeros ( 4 , 8 ) ;
Bv ( 1 , 1:2:7 )  =  1/3*dN_dx ( 1 , : ) ;  % dN/dx
Bv ( 1 , 2:2:8 )  =  1/3*dN_dx ( 2 , : ) ;  % dN/dy
Bv ( 2 , 1:2:7 )  =  1/3*dN_dx ( 1 , : ) ;  % dN/dx
Bv ( 2 , 2:2:8 )  =  1/3*dN_dx ( 2 , : ) ;  % dN/dy
Bv ( 3 , 1:2:7 )  =  1/3*dN_dx ( 1 , : ) ;  % dN/dx
Bv ( 3 , 2:2:8 )  =  1/3*dN_dx ( 2 , : ) ;  % dN/dy
% Bd ( 4 , 1:2:7 )  =  dN_dx ( 2 , : ) ;  % dN/dy
% Bd ( 4 , 2:2:8 )  =  dN_dx ( 1 , : ) ;  % dN/dx




end