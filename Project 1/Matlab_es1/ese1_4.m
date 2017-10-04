%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       EXERCISE 1.4 - UNDAMPED FREE RESPONSE IN CRITICAL CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The present code computes the components of the undamped free response
%  for Lambda --> Lambda_cr. For the procedure, see EXERCISE 1.2.
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

% No damping

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

% 1.4.1 UNDAMPED FREE RESPONSE 

       % Critical value of the parameter
       lambda_ND_cr = abs((omega_1^2-omega_2^2)/2);
       
% The components of the free response are computed for Lambda = perc*Lambda_cr ( choose a value close to 1 ) 
perc =0.999 ;
% Current value of Lambda 
lambda_ND = perc*lambda_ND_cr

       % Coefficients of the characteristic polynomial (no damping)
       a4 = 1;
       a3 = 0;
       a2 = omega_1^2+omega_2^2;
       a1 = 0;
       a0 = omega_1^2*omega_2^2+lambda_ND^2;
         
       polinomio_ND = [a4 a3 a2 a1 a0];
       
       % Roots of the characteristic polynomial
       s_ND = roots(polinomio_ND)

       % Matrix of eigenvectors ( 2 components * 4 roots )
       eigenvectors_ND = zeros(2,4);
       
       % Computation of the eigenvectors
       for i = 1:4
           eigenvectors_ND(1:2,i) = [-lambda_ND/(s_ND(i)^2+omega_1^2);1];
       end;
       eigenvectors_ND 
       
       % Computation of the constants depending on the initial conditions 
       matrix(1:2,1:4) = eigenvectors_ND(1:2,1:4);
       matrix(3:4,1:4) = [eigenvectors_ND(1,1)*s_ND(1,1) eigenvectors_ND(1,2)*s_ND(2,1) eigenvectors_ND(1,3)*s_ND(3,1) eigenvectors_ND(1,4)*s_ND(4,1);s_ND(1,1) s_ND(2,1) s_ND(3,1) s_ND(4,1)];
       const_ND = matrix^(-1)*[IC_1;IC_2;IC_3;IC_4]

% Minimum value of time ( seconds )
t_min =0 ;
% Maximum value of time ( seconds ) 
t_max =0.2 ;
% Delta t ( seconds )
dt =0.0002 ;
% Time discretization
t = t_min:dt:t_max;
% Number of roots 
n_roots = 4; 
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF THE COMPONENTS OF THE FREE RESPONSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
       w_1 = 0;
       w_2 = 0;
       for i = 1:n_roots
         % Components of the free response  
         w_1 = w_1 + const_ND(i)*eigenvectors_ND(1,i)*exp(s_ND(i)*t);
         w_2 = w_2 + const_ND(i)*exp(s_ND(i)*t);
       end
       
       % First component of the free response 
       figure(1)
       plot(t,w_1,'b')
       grid on
       xlabel('t [s]');
       ylabel('w_{1} (t)');
       title('Component of the undamped free response w_{1}(t) ( \Lambda --> \Lambda_{cr} )');
       
       % Second component of the free response
       figure(2)
       plot(t,44.06*t,'g')
       hold on
       plot(t,w_2,'r')
       grid on 
       xlabel('t [s]');
       ylabel('w_{2} (t)');
       title('Component of the undamped free response w_{2}(t) ( \Lambda --> \Lambda_{cr} )');
       
       figure(3)
       plot(t,w_1,'b')
       hold on
       plot(t,w_2,'r')
       grid on
       xlabel('t [s]');
       ylabel('w_ (t)');
       legend('w_{1}(t)','w_{2}(t)')
       title('Component of the undamped free response w_{1}(t) and w_{2}(t) ( \Lambda --> \Lambda_{cr} )');