clear all; close all

%% set up the model
nu = 0.4999; xint = 32; yint = 8;L = 48; H = 12;k = 0;
Model = model_setup(nu,xint,yint);  % Set up your model in "model_setup.m" to be read here.

%% discretize
Mesh = sub_discretization ( Model );

%% detect boundary nodes
BC = sub_get_boundary ( Model, Mesh );
C = {'r' 'b' '#EDB120' 'k'};
%% assembly loop for sigmaxx
fig = figure 
hold on
for i = 1:4

if i == 3
    [ K , f ]  =  sub_assembly_SRI ( Model , Mesh , BC);
    [ d ] = sub_solution ( K , f , BC );
elseif i ==4
    u_exact  =   Model.exact.displ ( Mesh.x_node(:,1), Mesh.x_node(:,2) );
    d= u_exact';
    d=d(:); 
else
int_pts = i; %% 1 for RI, 2 for FI
[ K , f ]  =  sub_assembly_FIRI_HC ( Model , Mesh , BC, int_pts,k);
[ d ] = sub_solution ( K , f , BC );
end
% postprocess
[displ, strain , stress, x_eval] = sub_postprocess ( Model , Mesh , d );

% sigmaxx
x0_ind = find(abs(x_eval(:,1) - L/2) < 1e-2);
plot(x_eval(x0_ind,2),stress(x0_ind,1),'color',C{i},"LineWidth",2)
xlabel("y-axis(inch)",'FontSize',14)
ylabel("\sigma_{xx} (Pa)",'FontSize',14)
title(sprintf('sigma_{xx} at x=L/2, %d * %d, nu = %0.4f ',xint,yint,nu),'FontSize',14)
% ylim([-0.3e8 0.3e8])

end
legend(["RI, ks = 0", "FI, ks= 0", "SRI",'Exact'],'FontSize',12)
% filename = sprintf('S11_%d_%d_nu_%0.4f.jpg',xint,yint,nu);
% saveas(fig,filename);

%% assembly loop for sigmaxy plot
fig= figure 
hold on
for i = 1:4

if i == 3
    [ K , f ]  =  sub_assembly_SRI ( Model , Mesh , BC);
    [ d ] = sub_solution ( K , f , BC );
elseif i ==4
    u_exact  =   Model.exact.displ ( Mesh.x_node(:,1), Mesh.x_node(:,2) );
    d= u_exact';
    d=d(:); 
else
int_pts = i; %% 1 for RI, 2 for FI
[ K , f ]  =  sub_assembly_FIRI_HC ( Model , Mesh , BC, int_pts,k);
[ d ] = sub_solution ( K , f , BC );
end
% postprocess
[displ, strain , stress, x_eval] = sub_postprocess ( Model , Mesh , d );

% sigmaxy
x0_ind = find(abs(x_eval(:,1) - L/2) < 1e-2);
plot(x_eval(x0_ind,2),stress(x0_ind,3),'color',C{i},"LineWidth",2)
end
xlabel("y-axis(inch)",'FontSize',14)
ylabel("\sigma_{xy} (Pa)",'FontSize',14)
title(sprintf('sigma_{xy} at x = L/2, %d * %d, nu = %0.4f ',xint,yint,nu),'FontSize',14)
legend(["RI, ks = 0", "FI, ks= 0", "SRI",'Exact'],'FontSize',12)
% filename = sprintf('S12_%d_%d_nu_%0.4f.jpg',xint,yint,nu);
% saveas(fig,filename);
% ylim([-0.5e4 inf])


%% assembly loop for uy plot
fig = figure 
hold on
for i = 1:4

if i == 3
    [ K , f ]  =  sub_assembly_SRI ( Model , Mesh , BC);
    [ d ] = sub_solution ( K , f , BC );
elseif i ==4
    u_exact  =   Model.exact.displ ( Mesh.x_node(:,1), Mesh.x_node(:,2) );
    d= u_exact';
    d=d(:); 
else
int_pts = i; %% 1 for RI, 2 for FI
[ K , f ]  =  sub_assembly_FIRI_HC ( Model , Mesh , BC, int_pts,k);
[ d ] = sub_solution ( K , f , BC );
end
% postprocess
[displ, strain , stress, x_eval] = sub_postprocess ( Model , Mesh , d );


% uy
y0_ind = find(abs(x_eval(:,2)) < 1e-2);
plot(x_eval(y0_ind,1),displ(y0_ind,2),'color',C{i},"LineWidth",2)
end
xlabel("x-axis(inch)",'FontSize',14)
ylabel("u_{y} (inch)",'FontSize',14)
title(sprintf('u_{y} at y=0, %d * %d, nu = %0.4f ',xint,yint,nu),'FontSize',14)
legend(["RI, ks = 0", "FI, ks= 0", "SRI",'Exact'],'FontSize',12)
% filename = sprintf('uy_%d_%d_nu_%0.4f.jpg',xint,yint,nu);
% saveas(fig,filename);

