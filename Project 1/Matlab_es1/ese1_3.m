%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         EXERCISE 1.3 - DRIVEN RESPONSE BELOW THE FLUTTER SPEED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The present code computes the driven response ( no I.C. ) to an 
%  impulsive gust load for Lambda < Lambda_cr.
%  The roots of the characteristic polynomial and the constant to write
%  the driven response are first evaluated. Then, the driven response is 
%  computed and presented in 3D and 2D plots. 
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

% 1.3.1 QUANTITIES TO COMPUTE THE DRIVEN RESPONSE (DR)

       if (zeta_1>0 && zeta_2>0)           
           % Stability margin with damping    
           lambda_cr = abs(((zeta_1*zeta_2*(omega_2^2-omega_1^2)^2)/((zeta_1+zeta_2)^2)+zeta_1*zeta_2*(zeta_1*omega_2^2+zeta_2*omega_1^2)/(zeta_1+zeta_2))^(1/2))      
       % To prevent divide-by-zero errors
       else           
           % Stability margin with no damping
           lambda_cr = abs((omega_1^2-omega_2^2)/2)      
       end
       
% The driven response is computed for Lambda = perc*Lambda_cr ( 0 < perc < 1)
perc = 0; 
% Current value of Lambda 
lambda_DR = perc*lambda_cr 

       % Coefficients of the characteristic polynomial
       a4 = 1;
       a3 = zeta_1+zeta_2;
       a2 = omega_1^2+omega_2^2+zeta_1*zeta_2;
       a1 = zeta_1*(omega_2^2)+zeta_2*(omega_1^2);
       a0 = omega_1^2*omega_2^2+(lambda_DR)^2;
         
       polynomial_LWD_DR = [a4 a3 a2 a1 a0];
         
       % Roots of the characteristic polynomial
       s_DR = roots(polynomial_LWD_DR)

       % Constants to compute the first component of the driven response
       A(1,1) = (s_DR(1,1)^2+zeta_2*s_DR(1,1)+omega_2^2)/((s_DR(1,1)-s_DR(2,1))*(s_DR(1,1)-s_DR(3,1))*(s_DR(1,1)-s_DR(4,1)));
       A(2,1) = (s_DR(2,1)^2+zeta_2*s_DR(2,1)+omega_2^2)/((s_DR(2,1)-s_DR(1,1))*(s_DR(2,1)-s_DR(3,1))*(s_DR(2,1)-s_DR(4,1)));
       A(3,1) = (s_DR(3,1)^2+zeta_2*s_DR(3,1)+omega_2^2)/((s_DR(3,1)-s_DR(2,1))*(s_DR(3,1)-s_DR(1,1))*(s_DR(3,1)-s_DR(4,1)));
       A(4,1) = (s_DR(4,1)^2+zeta_2*s_DR(4,1)+omega_2^2)/((s_DR(4,1)-s_DR(2,1))*(s_DR(4,1)-s_DR(3,1))*(s_DR(4,1)-s_DR(1,1)));

       % Constants to compute the second component of the driven response
       B(1,1) = 1/((s_DR(1,1)-s_DR(2,1))*(s_DR(1,1)-s_DR(3,1))*(s_DR(1,1)-s_DR(4,1)));
       B(2,1) = 1/((s_DR(2,1)-s_DR(1,1))*(s_DR(2,1)-s_DR(3,1))*(s_DR(2,1)-s_DR(4,1)));
       B(3,1) = 1/((s_DR(3,1)-s_DR(2,1))*(s_DR(3,1)-s_DR(1,1))*(s_DR(3,1)-s_DR(4,1)));
       B(4,1) = 1/((s_DR(4,1)-s_DR(2,1))*(s_DR(4,1)-s_DR(3,1))*(s_DR(4,1)-s_DR(1,1)));

% 1.3.2 3D DRIVEN RESPONSE

% Minimum value of time ( seconds )
t_min =0 ;
% Maximum value of time ( seconds )
t_max =300 ;
% Delta t ( seconds )
dt =0.4 ;
% Delta x
dx =0.001 ;
% Space and time discretization for 3D plot
[x,t] = meshgrid(0:dx:1,t_min:dt:t_max);
% Eigenfunctions for 3D plot
eigenfunction_1_3D = sin(pi*x);
eigenfunction_2_3D = sin(2*pi*x);
% Number of roots 
n_roots = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF THE DRIVEN RESPONSE FOR 3D PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       w_1 = 0;
       w_2 = 0;
       for i = 1:n_roots
         % Components of the driven response
         w_1 = w_1 + (2/pi)*A(i)*exp(s_DR(i)*t);
         w_2 = w_2 + (2*lambda_DR/pi)*B(i)*exp(s_DR(i)*t);
       end
       d_r_3D = w_1.*eigenfunction_1_3D + w_2.*eigenfunction_2_3D;

       % First component of the driven response
       figure(1)
       plot(t,w_1,'b')
       grid on
       xlabel('t [s]')
       ylabel('w_{1} (t)')
       title('Component of the driven response w_1(t) ( \Lambda < \Lambda_{cr} )')
         
       % Second component of the driven response
       figure(2)
       plot(t,w_2,'r')
       grid on 
       xlabel('t [s]')
       ylabel('w_{2}(t)')
       title('Component of the driven response w_2(t) ( \Lambda < \Lambda_{cr} )')
         
       % 3D plot
       figure(3)
       mesh(x,t,real(d_r_3D))
       grid on
       xlabel('x/L ')
       ylabel('t [s]')
       zlabel('w (x,t)')
       title('Driven response in space and time ( \Lambda < \Lambda_{cr} )')

% 1.3.3 2D DRIVEN RESPONSE

% Initial time ( seconds )
t_min =0 ; 
% Delta t ( seconds )
dt =30 ;
% Number of times 
n_t = 6;
% Maximum time ( seconds )
t_max = t_min + (n_t-1)*dt;
% Time discretization for 2D plot
t = t_min:dt:t_max;
% Space discretization for 2D plot
x_2D = 0:dx:1;
n_x = length(x_2D);
% Eigenfunctions for 2D plot
eigenfunction_1_2D = sin(pi*x_2D);
eigenfunction_2_2D = sin(2*pi*x_2D);
% Matrix to save the driven response
d_r = zeros(n_t,n_x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF THE DRIVEN RESPONSE FOR 2D PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         for i = 1:n_t
             w_1 = 0;
             w_2 = 0;
             for j = 1:n_roots
                % Components of the driven response at the current time
                w_1 = w_1 + (2/pi)*A(j)*exp(s_DR(j)*t(i));
                w_2 = w_2 + (2*lambda_DR/pi)*B(j)*exp(s_DR(j)*t(i));
             end
             % The results at t = t_min + (i-1)*dt are save in the i-th row
             d_r(i,:) = w_1*eigenfunction_1_2D + w_2*eigenfunction_2_2D;
         end;
       
       % 2D plot
       figure(4)
       plot(x_2D,d_r(1,:),'b',x_2D,d_r(2,:),'r',x_2D,d_r(3,:),'g',x_2D,d_r(4,:),'m',x_2D,d_r(5,:),'y',x_2D,d_r(6,:),'c')
       grid on
       legend('t_{0}','t_{1}','t_{2}','t_{3}','t_{4}','t_{5}')
       xlabel('x/L ')
       ylabel('w (x)')
       title('Driven response in space at different times ( \Lambda < \Lambda_{cr} )')
