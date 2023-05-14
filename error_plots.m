clear all; close all

%% set up the model
nu = 0.4999; xint = 16; yint = 4;L = 48; H = 12;

k = [0 0.004 0.01 0.02 0.05 0.1 0.2];
Model = model_setup(nu,xint,yint);  % Set up your model in "model_setup.m" to be read here.

%% discretize
Mesh = sub_discretization ( Model );

%% detect boundary nodes
BC = sub_get_boundary ( Model, Mesh );
% exact solution
u_exact  =   Model.exact.displ ( Mesh.x_node(:,1), Mesh.x_node(:,2) );
    d= u_exact';
    d=d(:); 
[displ, strain , stress, x_eval] = sub_postprocess ( Model , Mesh , d );
sigma_exact = stress;
u_exactsol = displ;


%% RI
for i = 1:7

int_pts = 1; 
[ K , f ]  =  sub_assembly_FIRI_HC ( Model , Mesh , BC, int_pts,k(i));
[ d ] = sub_solution ( K , f , BC );
[displ, strain , stress, x_eval] = sub_postprocess ( Model , Mesh , d );
sigmaRI{i} = stress;
uRI{i} = displ;
err_sigmaxx(i) = norm(stress(:,1) - sigma_exact(:,1));
err_sigmayy(i) = norm(stress(:,2) - sigma_exact(:,2));
err_sigmaxy(i) = norm(stress(:,3) - sigma_exact(:,3));
err_ux(i) = norm(displ(:,1) - u_exactsol(:,1));
err_uy(i) = norm(displ(:,2) - u_exactsol(:,2));
end

%% SRI
[ K , f ]  =  sub_assembly_SRI ( Model , Mesh , BC);
    [ d ] = sub_solution ( K , f , BC );

[displ, strain , stress, x_eval] = sub_postprocess ( Model , Mesh , d );
sigmaSRI = stress;
uSRI = displ;
errSRI_sigmaxx = norm(stress(:,1) - sigma_exact(:,1));
errSRI_sigmayy = norm(stress(:,2) - sigma_exact(:,2));
errSRI_sigmaxy = norm(stress(:,3) - sigma_exact(:,3));
errSRI_ux = norm(displ(:,1) - u_exactsol(:,1));
errSRI_uy = norm(displ(:,2) - u_exactsol(:,2));


figure 
hold on 
plot(k,err_ux,'r',"LineWidth",2);
plot(k,err_uy,'b',"LineWidth",2);
yline(errSRI_ux,"r--","Error norm of SRI for ux")
yline(errSRI_uy,"b--","Error norm of SRI for uy")
xlabel("stabilisation control parameter",'FontSize',14)
ylabel("error norm ",'FontSize',14)
title("displacement error norm at 16*4 mesh, nu = 0.4999",'FontSize',14)
legend(["for ux", "for uy"],'FontSize',12)


figure 
hold on 
plot(k,err_sigmaxx,'r',"LineWidth",2);
plot(k,err_sigmaxy,'b',"LineWidth",2);
yline(errSRI_sigmaxx,"r--","Error norm of SRI for sigmaxx")
yline(errSRI_sigmaxy,"b--","Error norm of SRI for sigmaxy")
xlabel("stabilisation control parameter",'FontSize',14)
ylabel("error norm ",'FontSize',14)
title("stress error norm at 16*4 mesh, nu = 0.4999",'FontSize',14)
legend(["for sigmaxx", "for sigmaxy"],'FontSize',12)


