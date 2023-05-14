function [Model] = model_setup(nu,xint,yint)
% Set up the model.
Model = struct;
%% material properties
E       =  3e7;
lambda  =  E * nu / (1+nu) / (1-2*nu);
mu      =  E / 2 / (1+nu);
Model.material.Youngs_modulus         =  E;
Model.material.Poissons_ratio         =  nu;
Model.material.Lames_first_parameter  =  lambda;
Model.material.shear_modulus          =  mu;

%% geometry of domain
% (x_min, y_max) ----------------- (x_max, y_max)
%                |               |
%                |               |
%                |               |
%                |               |
%                |               |
% (x_min, y_min) ----------------- (x_max, y_min)
x_min = 0; x_max = 48;
y_min = -6; y_max = +6;
D = y_max - y_min;
L = x_max - x_min;
I = D^3/12;
E_b = E/(1-nu^2);
nu_b = nu/(1-nu);

Model.domain.minmax = [[x_min; x_max], [y_min; y_max]];

%% mesh info
interval_x = xint; interval_y = yint;
Model.mesh.interval = [interval_x, interval_y];

%% Define boundary conditions
% One boundary segment can have different type of BCs (essential or
% natural) in different directions (x or y) in this code.
% The essential BC (EBC) is seperately defined for each direction.
% The natural BC (NBC) should be specified for both directions, but put
% zero for the direction in which an EBC is specified.
P = 40*10^3;
%% Essential boundaries - x direction
num_ebcx =  1;  % num of essential boundaries.
Model.bc.ebcx_position    =  cell ( 1 , num_ebcx ) ;
Model.bc.ebcx_value       =  cell ( 1 , num_ebcx ) ;
idx_bc = 0;
% x essential (Dirichlet) boundary 1
% The actual (x,y) values are plugged in after the domain is meshed.
idx_bc = idx_bc + 1;
Model.bc.ebcx_position{idx_bc} =  @(x,y) find( abs(x-x_min) < 1e-8 );  % where the condition below will be applied.
Model.bc.ebcx_value{idx_bc}    =   @(x,y) [(P*y/(6*E_b*I).*((2+nu_b)*(y.^2 - D^2/4))) ].* ones(size(x));
% Model.bc.ebcx_value{idx_bc}    =  @(x,y) zeros(size(x));
% displacement in x direction
%% Essential boundaries - y direction
num_ebcy =  1;  % num of essential boundaries.
Model.bc.ebcy_position    =  cell ( 1 , num_ebcy ) ;
Model.bc.ebcy_value       =  cell ( 1 , num_ebcy ) ;

idx_bc = 0;
% y essential (Dirichlet) boundary 1
idx_bc = idx_bc + 1;
Model.bc.ebcy_position{idx_bc} =  @(x,y) find( abs(x-x_min) < 1e-8 );  % where the condition below will be applied.
Model.bc.ebcy_value{idx_bc}    =  @(x,y) [-(P/(6*E_b*I)*(3*nu_b*y.^2.*(L))) ].* ones(size(y));
% displacement in y direction

%% Natural boundaries - for both directions
num_nbc    =  1;  % num of natural boundaries.
Model.bc.nbc_position    =  cell ( 1 , num_nbc ) ;
Model.bc.nbc_value       =  cell ( 1 , num_nbc ) ;

idx_bc = 0;
% natural (Neumann) boundary 1
% The actual (x,y) values are plugged in after the domain is meshed.
idx_bc = idx_bc + 1;
% Model.bc.nbc_position{idx_bc} = @(x,y) find( abs(x-x_min) < 1e-8 );  % where the condition below will be applied.
Model.bc.nbc_position{idx_bc} = @(x,y) find( abs(x-x_max) < 1e-8 );  
Model.bc.nbc_value{idx_bc}    = @(x,y) [ 0  ,  -(P/2/I)*(D^2/4 - y^2) ] .* ones(size(x));  % traction in x and y directions. Put zero if an essential BC is applied to that direction.


%% Body force
Model.body_force   =  @(x,y) [ 0 , 0 ] .* ones(size(x));

%% Exact solution - if there exists
Model.exact.use         =   1;  % 0 - exact solution not used ,   1 - used
Model.exact.displ       =   @(x,y) [ (P*y/(6*E_b*I).*((6*L - 3*x).*x + (2+nu_b)*(y.^2 - D^2/4)))  , ...
                                      -(P/(6*E_b*I)*(3*nu_b*y.^2.*(L-x) +(4+5*nu_b)*D^2*x/4  +(3*L - x).*x.^2 )) ] .* ones(size(x));


end