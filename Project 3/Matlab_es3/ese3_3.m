%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       EXERCISE 3.3 - ROOT LOCUS WITH QUASI-STEADY AERODYNAMICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The present code computes the root locus with quasi aerodynamics. 
%  The computation is also done in the case of structural damping.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace and close existing figures
clear all
close all
clc

%%%%%%%%%%%%%%%%%
% DATA
%%%%%%%%%%%%%%%%%

% Mass ratio
mass_ratio = 3 ;
% Frequency ratio
freq_ratio = 0.8 ;
% Non-dimensional position of the elastic center
xi_E =0.1 ;
% Non-dimensional position of the mass center
xi_G =0.3 ;
% Non-dimensional radius of gyration
r_alpha_2 =0.25 ;

% 3.3.1 SYSTEM WITHOUTH AERODYNAMIC STATES FOR C(K) = 1 AND WITH NO STRUCTURAL DAMPING (ND)

% Minimum value of U ( non-dimensional )
U_min_ND =0.1 ;
% Maximum value of U ( non-dimensional )
U_max_ND =10 ;
% Delta U ( non-dimensional )
dU_ND =0.0001 ;
% Vector of values of U 
U_ND = U_min_ND:dU_ND:U_max_ND;
% Number of elements
nU_ND = length(U_ND);
% Matrix of roots ( 4 * nU )
roots_ND = zeros(4,nU_ND);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF THE ROOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         for i=1:nU_ND;
             
               % Non-dimensional mass matrix 
               M = U_ND(i)^2*[1 xi_G;xi_G r_alpha_2-xi_E^2+2*xi_E*xi_G];
               % Non-dimensional damping matrix
               D = U_ND(i)^2*[2/mass_ratio 3/mass_ratio;0 1/mass_ratio];
               % Non-dimensional stiffness matrix
               K = freq_ratio^2*[1 xi_E+2*U_ND(i)^2/(mass_ratio*freq_ratio^2);xi_E xi_E^2+r_alpha_2/freq_ratio^2];
               
               % Coefficients of the characteristic polynomial          
               a4 = M(1,1)*M(2,2)-M(1,2)*M(2,1);
               a3 = D(2,2)*M(1,1)+D(1,1)*M(2,2)-D(1,2)*M(2,1)-D(2,1)*M(1,2);
               a2 = M(1,1)*K(2,2)+M(2,2)*K(1,1)-M(1,2)*K(2,1)-M(2,1)*K(1,2)+D(1,1)*D(2,2)-D(2,1)*D(1,2);
               a1 = D(2,2)*K(1,1)+D(1,1)*K(2,2)-D(2,1)*K(1,2)-D(1,2)*K(2,1);
               a0 = K(1,1)*K(2,2)-K(1,2)*K(2,1);
               polynomial = [a4 a3 a2 a1 a0]; 
               
               roots_ND(1:4,i) = roots(polynomial);
         
         end 

         % Root locus
         figure(1)
         plot(real(roots_ND),imag(roots_ND),'.b')
         grid on
         xlabel('Real part')
         ylabel('Imaginary part')
         title('Root locus with quasi-steady aerodynamics')
 
% 3.3.2 SYSTEM WITHOUTH AERODYNAMIC STATES WITH C(K) = 1 BUT WITH STRUCTURAL DAMPING (WD) 
 
% Parameter for structural damping
damp =0.005 ;

% Minimum value of U ( non-dimensional )
U_min_WD =0.6087 ;
% Maximum value of U ( non-dimensional )
U_max_WD =1.5 ;
% Delta U ( non-dimensional )
dU_WD =0.0001 ;
% Vector of values of U 
U_WD = U_min_WD:dU_WD:U_max_WD;
% Number of elements
nU_WD = length(U_WD);
% Matrix of roots ( 4 * nU )
roots_WD = zeros(4,nU_WD);

         % Non-dimensional structural damping matrix 
         SD = freq_ratio*2*damp*[1 0;0 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF THE ROOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         for i = 1:nU_WD
             
               % Non-dimensional mass matrix
               M = U_WD(i)^2*[1 xi_G;xi_G r_alpha_2-xi_E^2+2*xi_E*xi_G];
               % Non-dimensional aerodynamic damping matrix
               D = U_WD(i)^2*[2/mass_ratio 3/mass_ratio;0 1/mass_ratio];
               % Non-dimensional stiffness matrix
               K = freq_ratio^2*[1 xi_E+2*U_WD(i)^2/(mass_ratio*freq_ratio^2);xi_E xi_E^2+r_alpha_2/freq_ratio^2];
               % Overall non-dimensional damping matrix
               D_tot = D + SD; 
              
               % Coefficients of the characteristic polynomial
               a4 = M(1,1)*M(2,2)-M(1,2)*M(2,1);
               a3 = D_tot(2,2)*M(1,1)+D_tot(1,1)*M(2,2)-D_tot(1,2)*M(2,1)-D_tot(2,1)*M(1,2);
               a2 = M(1,1)*K(2,2)+M(2,2)*K(1,1)-M(1,2)*K(2,1)-M(2,1)*K(1,2)+D_tot(1,1)*D_tot(2,2)-D_tot(2,1)*D_tot(1,2);
               a1 = D_tot(2,2)*K(1,1)+D_tot(1,1)*K(2,2)-D_tot(2,1)*K(1,2)-D_tot(1,2)*K(2,1);
               a0 = K(1,1)*K(2,2)-K(1,2)*K(2,1);
               polynomial_WD = [a4 a3 a2 a1 a0]; 
               
               % The current roots are stored in the i-th row
               roots_WD(1:4,i) = roots(polynomial_WD);
               
         end

         % Root locus
         figure(2)
         plot(real(roots_WD),imag(roots_WD),'.b')
         grid on
         xlabel('Real part')
         ylabel('Imaginary part')
         title('Root locus with quasi-steady aerodynamics and structural damping')