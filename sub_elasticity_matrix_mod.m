function [ D ] = sub_elasticity_matrix_mod ( mat )

lambda  =   mat.Lames_first_parameter;
mu      =   mat.shear_modulus;
M       =   2*mu+lambda;

D       =   [      M , lambda , lambda,     0  ;
              lambda ,      M , lambda,     0  ;
              lambda , lambda ,      M,      0 ;
                   0 ,      0 ,      0,   mu  ] ;

end