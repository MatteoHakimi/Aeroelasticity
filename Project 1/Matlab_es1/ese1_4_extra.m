%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       EXERCISE 1.4.2 - DAMPED FREE RESPONSE IN CRITICAL CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The present code computes the components of the damped free response
%  for Lambda = Lambda_cr. For the procedure, see EXERCISE 1.2 and 1.4.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace and close figures
clear all
close all
clc

%%%%%%%%%%%%%%%%%
% DATA
%%%%%%%%%%%%%%%%%

% Natural frequency mode 1 ( Hz )
f_1 =15 ;
% Natural frequency mode 2 ( Hz )
f_2 =30 ;
% Damping mode 1
zeta_1 =0.03 ;
% Damping mode 2
zeta_2 =0.01 ;

% Natural frequencies ( rad/s )
omega_1 = f_1*2*pi;
omega_2 = f_2*2*pi;

%%%%%%%%%%%%%%%%%%%%%%
% INITIAL CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%

% Initial condition w1(0)
IC_1 =1 ;
% Initial condition w2(0)
IC_2 =0 ;
% Initial condition w1'(0)
IC_3 =0 ;
% Initial condition w2'(0)
IC_4 =0 ;

% 1.5.1 UNDAMPED FREE RESPONSE

       % Critical value of the parameter
       lambda_WD_cr = abs(((zeta_1*zeta_2*(omega_2^2-omega_1^2)^2)/((zeta_1+zeta_2)^2)+zeta_1*zeta_2*(zeta_1*omega_2^2+zeta_2*omega_1^2)/(zeta_1+zeta_2))^(1/2))

       % Coefficients of the characteristic polynomial
       a4 = 1;
       a3 = zeta_1+zeta_2;
       a2 = omega_1^2+omega_2^2+zeta_1*zeta_2;
       a1 = zeta_1*(omega_2^2)+zeta_2*(omega_1^2);
       a0 = omega_1^2*omega_2^2+(lambda_WD_cr)^2;
       
       polinomio_WD = [a4 a3 a2 a1 a0];
       
       % Roots of the characteristic polynomial
       s_WD = roots(polinomio_WD)

       % Matrix of eigenvectors ( 2 components * 4 roots )
       eigenvectors_WD = zeros(2,4);
       
       % Computation of the eigenvectors
       for i=1:4
           eigenvectors_WD(1:2,i) = [-lambda_WD_cr/(s_WD(i)^2+zeta_1*s_WD(i)+omega_1^2);1];
       end;
       eigenvectors_WD

       % Computation of the constans depending on the initial conditions
       matrix_WD(1:2,1:4) = eigenvectors_WD(1:2,1:4);
       matrix_WD(3:4,1:4) = [eigenvectors_WD(1,1)*s_WD(1,1) eigenvectors_WD(1,2)*s_WD(2,1) eigenvectors_WD(1,3)*s_WD(3,1) eigenvectors_WD(1,4)*s_WD(4,1);s_WD(1,1) s_WD(2,1) s_WD(3,1) s_WD(4,1)];
       const_WD = matrix_WD^(-1)*[IC_1;IC_2;IC_3;IC_4]

% Minimum value of time ( seconds )
t_min =0 ;
% Maximum value of time ( seconds ) 
t_max =10 ;
% Delta t ( seconds )
dt =0.01 ;
% Time discretization
t = t_min:dt:t_max;
% Number of roots 
n_roots = 4; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF THE COMPONENTS OF THE FREE RESPONSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       % Components of the free response
       w_1 = 0;
       w_2 = 0;
       for i = 1:n_roots
         w_1 = w_1 + const_WD(i)*eigenvectors_WD(1,i)*exp(s_WD(i)*t);
         w_2 = w_2 + const_WD(i)*exp(s_WD(i)*t);
       end 

       % First component of the free response 
       figure(1)
       plot(t,w_1,'-bx')
       grid on
       xlabel('t [s]')
       ylabel('w_1 (t)')
       title('Component of the damped free response w_1(t) ( \Lambda = \Lambda_{cr} )')
       
       % Second component of the free response
       figure(2)
       plot(t,w_2,'-rx')
       grid on 
       xlabel('t [s]');
       ylabel('w_2 (t)');
       title('Component of the damped free response w_2(t) ( \Lambda = \Lambda_{cr} )')