%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXERCISE 2.1 - LINEAR STABILITY WITH STEADY AND QUASI-STEADY AERODYNAMICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The present code computes the non-dimensional divergence velocity of the 
% typical section and the variation of the roots of the characteristic 
% polynomial with respect to U in the cases of steady and quasi-steady 
% aerodynamics. 
% The results are presented in the complex plane (root locus). 
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
mass_ratio =5 ;                        %velivolo aviazione generale 8.084
% Frequency ratio
freq_ratio = 0.5;                      %velivolo aviazione generale 0.5
% Non-dimensional position of the elastic center
xi_E =0.3 ;                            %velivolo aviazione generale 0.5
% Non-dimensional position of the mass center
xi_G =0.45 ;                           %velivolo aviazione generale 0.167
% Non-dimensional radius of gyration
r_alpha_2 =0.25 ;                      %velivolo aviazione generale 0.1365

% 2.1.1 SYSTEM MATRICES WITH STEADY AERODYNAMICS

       % Non-dimensional mass matrix
         M = [1 xi_G;xi_G r_alpha_2-xi_E^2+2*xi_E*xi_G]
       % Non-dimensional stiffness matrix
         K = freq_ratio^2*[1 xi_E;xi_E xi_E^2+r_alpha_2/(freq_ratio^2)]

       % Non-dimensional matrix due to the steady aerodynamic loads
         SA = [0 2/mass_ratio;0 0]
         
% 2.1.2 NON-DIMENSIONAL DIVERGENCE VELOCITY         

         % The divergence condition is in the form a*U_D^4 + b*U_D^2 + c = 0 
         a =  SA(1,1)*SA(2,2)-SA(1,2)*SA(2,1);                              
         b = (SA(1,1)*K(2,2)+SA(2,2)*K(1,1)-SA(1,2)*K(2,1)-SA(2,1)*K(1,2)); 
         c = (K(1,1)*K(2,2)-K(1,2)*K(2,1));                                  
         
         polynomial_D = [a b c];
         
         % Divergence velocity
         U_D = (roots(polynomial_D))^0.5

% 2.1.3 ROOT LOCUS WITH STEADY AERODYNAMICS

% Maximum value of the parameter U_SA ( SA = Steady Aerodynamics )
U_SA_max =2 ;
% Delta U_SA
d_U_SA =0.0001 ;
% Vector of values of the parameter
U_SA = 0:d_U_SA:U_SA_max;
% Number of elements
n_U_SA = length(U_SA);
% Matrix of roots ( 4 roots * n U_SA)
s_SA = zeros(4,n_U_SA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF THE ROOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These coefficients of the characteristic polynomial do not depend on U_SA
a4 = M(1,1)*M(2,2)-M(1,2)*M(2,1);
a3 = 0;
a1 = 0;

         for i = 1:n_U_SA
             
             % These coefficients of the characteristic polynomial depend on U_SA
             a2 = M(1,1)*(K(2,2)+U_SA(i)^2*SA(2,2))+M(2,2)*(K(1,1)+U_SA(i)^2*SA(1,1))-M(1,2)*(K(2,1)+U_SA(i)^2*SA(2,1))-M(2,1)*(K(1,2)+U_SA(i)^2*SA(1,2));
             a0 = (K(1,1)+U_SA(i)^2*SA(1,1))*(K(2,2)+U_SA(i)^2*SA(2,2))-(K(1,2)+U_SA(i)^2*SA(1,2))*(K(2,1)+U_SA(i)^2*SA(2,1));
             
             polynomial_AS = [a4 a3 a2 a1 a0];
             
             % The current roots are stored in the i-th column of the matrix
             s_SA(1:4,i) = roots(polynomial_AS);

         end                                                                       
         
         % Root locus with steady aerodynamics
         figure(1)
         plot(real(s_SA),imag(s_SA),'.b')
         grid on 
         hold on
         % Roots with no air ( U_SA = 0 )
         plot(real(s_SA(:,1)),imag(s_SA(:,1)),'ko','MarkerSize',8,'LineWidth',2)
         xlabel('Real part')
         ylabel('Imaginary part')
         title('Root locus with steady aerodynamics and no structural damping')
     
          
%  2.1.4 ROOT LOCUS WITH QUASI-STEADY AERODYNAMICS AND STRUCTURAL DAMPING 

% Parameter for structural damping
damp =0.005 ;

         % Non-dimensional structural damping matrix
         D = freq_ratio*2*damp*[1 0;0 1]

         % Insert the elements of the non-dimensional matrix due to the quasi-steady aerodynamic loads
         QSA = [2/mass_ratio,0;0 0]

% Maximum value of the parameter U_QSA ( QSA = Quasi-Steady Aerodynamics )
U_QSA_max =2 ;
% Delta U_QSA
d_U_QSA =0.0001 ;
% Vector of values of the parameter
U_QSA = 0:d_U_QSA:U_QSA_max;
% Number of elements
n_U_QSA = length(U_QSA);
% Matrix of roots ( 4 roots * n U_QSA)
s_QSA = zeros(4,n_U_QSA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF THE ROOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This coefficient of the characteristic polynomial does not depend on U_QSA
a4 = M(1,1)*M(2,2)-M(1,2)*M(2,1);
  
           for i=1:n_U_QSA;
             
             % These coefficients of the characteristic polynomial depend on U_QSA
             a3 = (U_QSA(i)*QSA(2,2)+D(2,2))*M(1,1)+(U_QSA(i)*QSA(1,1)+D(1,1))*M(2,2)-(U_QSA(i)*QSA(1,2)+D(1,2))*M(2,1)-(U_QSA(i)*QSA(2,1)+D(2,1))*M(1,2);
             a2 = M(1,1)*(K(2,2)+U_QSA(i)^2*SA(2,2))+M(2,2)*(K(1,1)+U_QSA(i)^2*SA(1,1))-M(1,2)*(K(2,1)+U_QSA(i)^2*SA(2,1))-M(2,1)*(K(1,2)+U_QSA(i)^2*SA(1,2))+(U_QSA(i)*QSA(1,1)+D(1,1))*(U_QSA(i)*QSA(2,2)+D(2,2))-(U_QSA(i)*QSA(2,1)+D(2,1))*(U_QSA(i)*QSA(1,2)+D(1,2));
             a1 = (U_QSA(i)*QSA(2,2)+D(2,2))*(K(1,1)+U_QSA(i)^2*SA(1,1))+(U_QSA(i)*QSA(1,1)+D(1,1))*(K(2,2)+U_QSA(i)^2*SA(2,2))-(U_QSA(i)*QSA(2,1)+D(2,1))*(K(1,2)+U_QSA(i)^2*SA(1,2))-(U_QSA(i)*QSA(1,2)+D(1,2))*(K(2,1)+U_QSA(i)^2*SA(2,1));
             a0 = (K(1,1)+U_QSA(i)^2*SA(1,1))*(K(2,2)+U_QSA(i)^2*SA(2,2))-(K(1,2)+U_QSA(i)^2*SA(1,2))*(K(2,1)+U_QSA(i)^2*SA(2,1));
             
             polynomial_QSA = [a4 a3 a2 a1 a0]; 
             
             % The current roots are stored in the i-th column
             s_QSA(1:4,i) = roots(polynomial_QSA);
          
           end 

         % Root locus with quasi-steady aerodynamics and structural damping
         figure(2)
         plot(real(s_QSA),imag(s_QSA),'.b')
         grid on 
         hold on
         % Roots with no air ( U_QSA = 0 )
         plot(real(s_QSA(:,1)),imag(s_QSA(:,1)),'ko','MarkerSize',8,'LineWidth',2)
         xlabel('Real part')
         ylabel('Imaginary part')
         title('Root locus with quasi-steady aerodynamics and structural damping')
         