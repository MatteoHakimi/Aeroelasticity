%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           EXERCISE 3.2 - ROOT LOCUS WITH UNSTEADY AERODYNAMICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The present code computes the root locus with unsteady aerodynamics. 
%  The system is written in first-order form using the approximated 
%  Theodorsen complex function (Padé) and introducing an aerodynamic state. 
%  The latter is associated to the proper part of the rational function. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace and close existing figures
clear all
close all
clc

%%%%%%%%%%%%%%%%%
% DATA
%%%%%%%%%%%%%%%%%

% % Mass ratio
% mass_ratio =5 ;
% % Frequency ratio
% freq_ratio =0.5 ;
% % Non-dimensional position of the elastic center
% xi_E =0.3 ;
% % Non-dimensional position of the mass center
% xi_G =0.45 ;
% % Non-dimensional radius of gyration
% r_alpha_2 =0.25 ;

% Mass ratio
mass_ratio = 3 ;
% Frequency ratio
freq_ratio = 0.8 ;
% Non-dimensional position of the elastic center
xi_E = 0.1 ;
% Non-dimensional position of the mass center
xi_G = 0.3 ;
% Non-dimensional radius of gyration
r_alpha_2 = 0.25 ;


% 3.2.1 SYSTEM WITH ADDED AERODYNAMIC STATES 

% Perfoming the polynomial division, Padé approximation 
% for Theodorsen function is rewritten as 
%                C(p) = C_0 + N_1(p)/D_2(P)
% where the rational part is PROPER. 
% In particular:
%                N_1(p) = coef_N_0 + p*coeff_N_1 and 
%                D_2(p) = coef_D_0 + p*coeff_D_1 + p^2*coef_D_2
% The coefficients must be inserted below. 

% Insert C_0
C_0 =0.5 ;
% Insert the coefficients of N_1(p)
coef_N_0 =0.022 ;
coef_N_1 =0.117 ;
% Insert the coefficients of D_2(p)
coef_D_0 =0.044 ;
coef_D_1 =0.552 ;
coef_D_2 =1 ;

% Minimum value of U ( non-dimensional )
U_min =0.1 ;
% Maximum value of U ( non-dimensional )
U_max = 10 ;
% Delta U ( non-dimensional )
dU =0.0001 ;
% Vector of values of U 
U = U_min:dU:U_max;
% Number of elements
nU = length(U);
% Matrix of roots ( 6 * nU )
roots = zeros(6,nU);
% Order of the differential system (written in second-order form)
n = 3;
% Identity matrix (n*n)
I = eye(n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF THE ROOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The new mass, damping and stiffness matrices are 3*3 

% The mass matrix is written as M = U^2*[m11 m12 0; m21 m22 0; m31 m32 m33]
% These elements do not depend on U
m11 = 1+1/mass_ratio;
m12 = xi_G+1/(2*mass_ratio);
m21 = xi_G+1/(2*mass_ratio);
m22 = r_alpha_2-xi_E^2+2*xi_E*xi_G+3/(8*mass_ratio);

% The damping matrix is written as D = U^2[d11 d12 0; 0 1 0; d31 d32 d33]/mass_ratio
% These elements do not depend on U
d11 = 2*C_0;
d12 = 1+2*C_0;

% The stiffness matrix is written as K = freq_ratio^2*[1 k12 k13;k21 k22 0;0 k32 k33]
% These elements do not depend on U
k21 = xi_E;
k22 = xi_E^2+r_alpha_2/freq_ratio^2;
k32 = -2*coef_N_0/(mass_ratio*freq_ratio^2); 
k33 = coef_D_0/freq_ratio^2; 
      
         for k = 1:nU

               % Non trivial elements of the new non-dimensional mass matrix (3*3) that depends on U
               m31 = -2*coef_N_1/(mass_ratio*U(k)^2);
               m32 = -2*coef_N_1/(mass_ratio*U(k)^2);
               m33 = coef_D_2/U(k)^2;
               % New non-dimensional mass matrix (3*3)
               M=U(k)^2*[m11 m12 0;m21 m22 0;m31 m32 m33];
               
               % Non trivial elements of the new non-dimensional damping matrix (3*3) that depend on U
               d31 = -2*coef_N_0/U(k)^2;
               d32 = -2*(coef_N_1+coef_N_0)/U(k)^2;
               d33 = coef_D_1*mass_ratio/U(k)^2;
               % New non-dimensional damping matrix (3*3)
               D = U(k)^2*[d11 d12 0;0 1 0;d31 d32 d33]/mass_ratio;
             
               % Non trivial elements of the new non-dimensional stiffness matrix (3*3) that depend on U
               k12 = xi_E+2*C_0*U(k)^2/(mass_ratio*freq_ratio^2);
               k13 = U(k)^2/freq_ratio^2;
               % New non-dimensional stiffness matrix (3*3)
               K=freq_ratio^2*[1 k12 k13;k21 k22 0;0 k32 k33];     
             
               % State matrix (6*6)
               A = [zeros(n,n) I(1:n,1:n);-M^(-1)*K -M^(-1)*D];
                
               % The current eigenvalues are stored in the k-th column
               roots(1:2*n,k) = eig(A);
               
         end
          
         % Root locus
         figure(1)
         plot(real(roots),imag(roots),'.b')
         grid on
         hold on
         xlabel('Real part')
         ylabel('Imaginary part')
         title('Root locus with unsteady aerodynamics')
%           axis ([-1.2 0.4 -2 2])  