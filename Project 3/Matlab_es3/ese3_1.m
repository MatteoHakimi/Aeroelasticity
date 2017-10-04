%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             EXERCISE 3.1 - STABILITY MARGIN USING U^2=f(k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The present code compares Theodorsen complex function (Bessel functions)
%  with the rational polynomial approximation (Padé). 
%  Both forms of Theodorsen function are used to compute the stability
%  margin using the characteristic equation of the system in the form
%  U^2=f(k). 
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
mass_ratio = 5 ;
% Frequency ratio
freq_ratio =0.5 ;
% Non-dimensional position of the elastic center
xi_E = 0.3;
% Non-dimensional position of the mass center
xi_G = 0.45;
% Non-dimensional radius of gyration
r_alpha_2 =0.25 ;

% 3.1.1 THEODORSEN FUNCTION

% Minimum value of the reduced frequency
k_min = 0;
% Maximum value of the reduced frequency
k_max =2 ;
% Delta k
dk =0.01 ;
% Vector of values of k
k = k_min:dk:k_max;
% Number of elements
nk = length(k);
% Vectors to store Theodorsen complex function
C_approx = zeros(nk,1);
C_exact = zeros(nk,1);
% Vector for Hankel functions
H = zeros(2,nk);

         % Computation of Theodorsen function
         for m = 1:nk
             % Padé approximation
             C_approx(m) = 0.5*(1i*k(m)+0.135)*(1i*k(m)+0.651)/((1i*k(m)+0.0965)*(1i*k(m)+0.4555));
             % Hankel functions
             for n=1:2
                 H(n,m) = besselh(n-1,2,k(m));
             end;
             % Exact Theodorsen function
             C_exact(m) = H(2,m)/(H(2,m)+H(1,m)*1i);
         end;
         C_exact(1) = 1;

         % Comparison of the real parts
         figure(1)
         plot(k,real(C_exact),'ro',k,real(C_approx),'bx')
         grid on 
         legend('Bessel functions','Padé approximation')      
         xlabel('k [-]')
         ylabel('Re (C) [-]')
         title('Real part of Theodorsen complex function vs reduced frequency')
         
         % Comparison of the imaginary parts
         figure(2)
         plot(k,imag(C_exact),'ro',k,imag(C_approx),'bx')
         grid on
         legend('Bessel functions','Padé approximation')      
         xlabel('k [-]')
         ylabel('Im (C) [-]')
         title('Imaginary part of Theodorsen complex function vs reduced frequency')

% 3.1.2 STABILITY MARGINS

       % Non-dimensional mass matrix 
       M = [1+1/mass_ratio xi_G+1/(2*mass_ratio);xi_G+1/(2*mass_ratio) r_alpha_2-xi_E^2+2*xi_E*xi_G+3/(8*mass_ratio)];
       % Non-dimensional stiffness matrix
       K = freq_ratio^2*[1 xi_E;xi_E xi_E^2+r_alpha_2/(freq_ratio^2)];
       % Non-dimensional damping matrix 
       D = [0 1/mass_ratio;0 1/mass_ratio];

% The computation of U^2 = f(k) is performed for the values of k defined above (apart from k = 0)
% Two solutions f_1(k) and f_2(k) are obtained. 

% Vectors to store the exact values of f(k) ( using C_exact )
f_1_exact = zeros(nk-1,1);
f_2_exact = zeros(nk-1,1);
% Vectors to store the approximated values of f(k) ( using C_approx )
f_1_approx = zeros(nk-1,1);
f_2_approx = zeros(nk-1,1);

         for m = 1:nk
             
             % Restriction to the imaginary axis
             p = k(m)*1i;
             
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % Case 1 - exact Theodorsen function
             C = C_exact(m);
             % Non-dimensional matrix due to the unsteady aerodynamic loads
             FUA = [C*p*2/mass_ratio (1+p)*C*2/mass_ratio;0 0];
             
             % Computation of the solution f_1 
             c4 = (M(1,1)*p^2+D(1,1)*p+FUA(1,1))*(M(2,2)*p^2+D(2,2)*p+FUA(2,2))-(M(1,2)*p^2+D(1,2)*p+FUA(1,2))*(M(2,1)*p^2+D(2,1)*p+FUA(2,1));
             c2 = (M(1,1)*p^2+D(1,1)*p+FUA(1,1))*K(2,2)+(M(2,2)*p^2+D(2,2)*p+FUA(2,2))*K(1,1)-(M(2,1)*p^2+D(2,1)*p+FUA(2,1))*K(1,2)-(M(1,2)*p^2+D(1,2)*p+FUA(1,2))*K(2,1);
             c0 = K(1,1)*K(2,2)-K(2,1)*K(1,2);
             f_1_exact(m) = (-c2+(c2^2-4*c4*c0)^0.5)/(2*c4);
             
             % Computation of the solution f_2
             c4 = (M(1,1)*p^2+D(1,1)*p+FUA(1,1))*(M(2,2)*p^2+D(2,2)*p+FUA(2,2))-(M(1,2)*p^2+D(1,2)*p+FUA(1,2))*(M(2,1)*p^2+D(2,1)*p+FUA(2,1));
             c2 = (M(1,1)*p^2+D(1,1)*p+FUA(1,1))*K(2,2)+(M(2,2)*p^2+D(2,2)*p+FUA(2,2))*K(1,1)-(M(2,1)*p^2+D(2,1)*p+FUA(2,1))*K(1,2)-(M(1,2)*p^2+D(1,2)*p+FUA(1,2))*K(2,1);
             c0 = K(1,1)*K(2,2)-K(2,1)*K(1,2);
             f_2_exact(m) = (-c2-(c2^2-4*c4*c0)^0.5)/(2*c4);
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % Case 2 - Padé approximation
             C = C_approx(m);
             % Non-dimensional matrix due to the unsteady aerodynamic loads
             FUA = [C*p*2/mass_ratio (1+p)*C*2/mass_ratio;0 0];
             
             % Computation of the solution f_1 
             c4 = (M(1,1)*p^2+D(1,1)*p+FUA(1,1))*(M(2,2)*p^2+D(2,2)*p+FUA(2,2))-(M(1,2)*p^2+D(1,2)*p+FUA(1,2))*(M(2,1)*p^2+D(2,1)*p+FUA(2,1));
             c2 = (M(1,1)*p^2+D(1,1)*p+FUA(1,1))*K(2,2)+(M(2,2)*p^2+D(2,2)*p+FUA(2,2))*K(1,1)-(M(2,1)*p^2+D(2,1)*p+FUA(2,1))*K(1,2)-(M(1,2)*p^2+D(1,2)*p+FUA(1,2))*K(2,1);
             c0 = K(1,1)*K(2,2)-K(2,1)*K(1,2);
             f_1_approx(m) = (-c2+(c2^2-4*c4*c0)^0.5)/(2*c4);
             
             % Computation of the solution f_2
             c4 = (M(1,1)*p^2+D(1,1)*p+FUA(1,1))*(M(2,2)*p^2+D(2,2)*p+FUA(2,2))-(M(1,2)*p^2+D(1,2)*p+FUA(1,2))*(M(2,1)*p^2+D(2,1)*p+FUA(2,1));
             c2 = (M(1,1)*p^2+D(1,1)*p+FUA(1,1))*K(2,2)+(M(2,2)*p^2+D(2,2)*p+FUA(2,2))*K(1,1)-(M(2,1)*p^2+D(2,1)*p+FUA(2,1))*K(1,2)-(M(1,2)*p^2+D(1,2)*p+FUA(1,2))*K(2,1);
             c0 = K(1,1)*K(2,2)-K(2,1)*K(1,2);
             f_2_approx(m) = (-c2-(c2^2-4*c4*c0)^0.5)/(2*c4);
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
         end
         
%%%%%%%%%%%%%%%%%
% f(k) GRAPHS 
%%%%%%%%%%%%%%%%%

         % Real part of f_1 
         figure(3)
         plot(k(1:end),real(f_1_exact),'ro',k(1:end),real(f_1_approx),'bx')
         grid on 
         axis([0 k_max 0 5])
         legend('Bessel functions','Padé approximation')      
         xlabel('k [-]')
         ylabel('Re f_1(k) [-]')
         title('Real part of f_1(k) vs reduced frequency')
         
         % Real part of f_2 
         figure(4)
         plot(k(1:end),real(f_2_exact),'ro',k(1:end),real(f_2_approx),'bx')
         grid on
         legend('Bessel functions','Padé approximation')      
         xlabel('k [-]')
         ylabel('Re f_2(k) [-]')
         title('Real part of f_2(k) vs reduced frequency')

         % Imaginary part of f_1 
         figure(5)
         plot(k(1:end),imag(f_1_exact),'ro',k(1:end),imag(f_1_approx),'bx')
         grid on
         axis([0 k_max -1 0.5])
         legend('Bessel functions','Padé approximation')      
         xlabel('k [-]')
         ylabel('Im f_1(k) [-]')
         title('Imaginary part of f_1(k) vs reduced frequency')
         
         % Imaginary part of f_2 
         figure(6)
         plot(k(1:end),imag(f_2_exact),'ro',k(1:end),imag(f_2_approx),'bx')
         grid on 
         legend('Bessel functions','Padé approximation')      
         xlabel('k [-]')
         ylabel('Im f_2(k) [-]')
         title('Imaginary part of f_2(k) vs reduced frequency')
