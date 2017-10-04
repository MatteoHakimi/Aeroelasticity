%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         EXERCISE 1.2 - FREE RESPONSE BELOW THE FLUTTER SPEED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The present code computes the free response for Lambda < Lambda_cr. 
%  The roots of the characteristic polynomial, the eigenvectors, and the 
%  constants depending on the initial conditions are first computed. 
%  Then, the free response is obtained and presented in 3D and 2D plots. 
%  The code also computes the free response to the same initial conditions
%  in the case of no damping and no flow ( Lambda = 0 ). 
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
zeta_1 =0.01 ;
% Damping mode 2
zeta_2 =0.03 ;

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

% 1.2.1 STABILITY MARGIN AND CURRENT VALUE OF LAMBDA 

       if (zeta_1>0 && zeta_2>0)           
           % Stability margin with damping    
           lambda_cr = abs(((zeta_1*zeta_2*(omega_2^2-omega_1^2)^2)/((zeta_1+zeta_2)^2)+zeta_1*zeta_2*(zeta_1*omega_2^2+zeta_2*omega_1^2)/(zeta_1+zeta_2))^(1/2))      
       % To prevent divide-by-zero errors
       else           
           % Stability margin with no damping
           lambda_cr = abs((omega_1^2-omega_2^2)/2)      
       end

% The free response will be computed for Lambda = perc*Lambda_cr ( 0 < perc < 1)
perc =0.9 ; 
% Current value of Lambda 
lambda_FR = perc*lambda_cr 

% 1.2.2 DAMPED FREE RESPONSE 

       % Coefficients of the characteristic polynomial
       a4 = 1;
       a3 = zeta_1+zeta_2;
       a2 = omega_1^2+omega_2^2+zeta_1*zeta_2;
       a1 = zeta_1*(omega_2^2)+zeta_2*(omega_1^2);
       a0 = omega_1^2*omega_2^2+(lambda_FR)^2;
         
       polynomial_WD_FR = [a4 a3 a2 a1 a0];
         
       % Roots of the characteristic polynomial
       s_FR = roots(polynomial_WD_FR)

       % Matrix of eigenvectors ( 2 components * 4 roots )
       eigenvectors = zeros(2,4);
       
       % Computation of the eigenvectors
       for i=1:2
             eigenvectors(1:2,i) = [-lambda_FR/(s_FR(i,1)^2+zeta_1*s_FR(i,1)+omega_1^2);1];
             eigenvectors(1:2,i+2) = [1;-lambda_FR/(s_FR(i+2,1)^2+zeta_2*s_FR(i+2,1)+omega_2^2)];
       end;
         
       eigenvectors
          
       % Computation of the constants depending on the initial conditions 
       matrix(1:2,1:4) = eigenvectors(1:2,1:4);
       matrix(3:4,1:4) = [eigenvectors(1,1)*s_FR(1,1) eigenvectors(1,2)*s_FR(2,1) eigenvectors(1,3)*s_FR(3,1) eigenvectors(1,4)*s_FR(4,1);s_FR(1,1) s_FR(2,1) s_FR(3,1) s_FR(4,1)];
       const = matrix^(-1)*[IC_1;IC_2;IC_3;IC_4]
       
       % 1.2.2.1 3D FREE RESPONSE

% Minimum value of time ( seconds )
t_min =0 ;
% Maximum value of time ( seconds )
t_max =1 ;
% Delta t ( seconds )
dt =1 ;
% Delta x ( seconds )
dx =0.001 ;
% Space and time discretization for 3D plot
[x,t] = meshgrid(0:dx:1,t_min:dt:t_max);
% Eigenfunctions for 3D plot
eigenfunction_1_3D = sin(pi*x);
eigenfunction_2_3D = sin(2*pi*x);
% Number of roots
n_roots = 4;
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF THE FREE RESPONSE FOR 3D PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       w_1 = 0;
       w_2 = 0;     
       for i = 1:n_roots
           % Components of the free response
           w_1 = w_1+const(i)*eigenvectors(1,i)*exp(s_FR(i)*t);
           w_2 = w_2+const(i)*eigenvectors(2,i)*exp(s_FR(i)*t);
       end
       % 3D Free response 
       f_r_3D = w_1.*eigenfunction_1_3D + w_2.*eigenfunction_2_3D;
         
       % First component of the free response 
       figure(1)
       plot(t,w_1,'-b')
       grid on
       xlabel('t [s]')
       ylabel('w_1 (t)')
       title('Component of the free response w_1(t) ( \Lambda < \Lambda_{cr} )')
       
       % Second component of the free response
       figure(2)
       plot(t,w_2,'-r')
       grid on 
       xlabel('t [s]')
       ylabel('w_2 (t)')
       title('Component of the free response w_2(t) ( \Lambda < \Lambda_{cr} )')

       % 3D plot
       figure(3)
       mesh(x,t,real(f_r_3D))
       xlabel('x/L ')
       ylabel('t [s]')
       zlabel('w (x,t)')
       title('Free response in space and time ( \Lambda < \Lambda_{cr} )')

       % 1.2.2.2 2D FREE RESPONSE

% Initial time ( seconds )
t_min = 0 ; 
% Delta t ( seconds )
dt = 0.02 ;
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
% Matrix to save the free response at different times
f_r = zeros(n_t,n_x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF THE FREE RESPONSE FOR 2D PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         for i = 1:n_t
             % Components of the free response at the current time
             w_1 = 0;
             w_2 = 0;
             for j = 1:n_roots  
               w_1 = w_1 + const(j)*eigenvectors(1,j)*exp(s_FR(j)*t(i));
               w_2 = w_2 + const(j)*eigenvectors(2,j)*exp(s_FR(j)*t(i));
             end 
             % The results for t = (i-1)*dt are stored in the i-th row 
             f_r(i,:) = w_1*eigenfunction_1_2D + w_2*eigenfunction_2_2D;
         end;         
                                                         
         % 2D plot
         figure(4)
         plot(x_2D,f_r(1,:),'-m',x_2D,f_r(2,:),'-r',x_2D,f_r(3,:),'-g',x_2D,f_r(4,:),'-c',x_2D,f_r(5,:),'-b',x_2D,f_r(6,:),'-y')
         grid on
         legend('t_{0}','t_{1}','t_{2}','t_{3}','t_{4}','t_{5}')
         xlabel('x/L')
         ylabel('w (x)')
         title('Free response in space at different times ( \Lambda < \Lambda_{cr} )')

% 1.2.3 FREE RESPONSE WITH NO DAMPING AND NO FLOW (ND_NF)

       % 1.2.3.1 3D FREE RESPONSE

% Minimum value of time ( seconds )
t_min =0 ;
% Maximum value of time ( seconds )
t_max =0.2 ;
% Delta t ( seconds )
dt = 0.0001; 
% Space and time discretization for 3D plot
[x,t] = meshgrid(0:dx:1,t_min:dt:t_max);
% Eigenfunctions for 3D plot
eigenfunction_1_3D = sin(pi*x);
eigenfunction_2_3D = sin(2*pi*x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF THE FREE RESPONSE FOR 3D PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       % Components of the free response
       w_1_ND_NF = cos(omega_1*t);
       w_2_ND_NF = cos(omega_2*t); 
       % 3D Free response   
       f_r_3D_ND_NF = IC_1*w_1_ND_NF.*eigenfunction_1_3D + IC_2*w_2_ND_NF.*eigenfunction_2_3D;
                
       % 3D plot
       figure(5)
       mesh(x,t,real(f_r_3D_ND_NF))
       grid on
       xlabel('x/L ')
       ylabel('t [s]')
       zlabel('w (x,t)')
       title('Undamped free response in space and time with no flow ( \Lambda = 0 )')
       
       % 1.2.3.2 2D FREE RESPONSE

% Initial time ( seconds )
t_min =0 ; 
% Delta t ( seconds )
dt =0.02 ;
% Number of times 
n_t = 6;
% Maximum time ( seconds )
t_max = t_min + (n_t-1)*dt;
% Time discretization for 2D plot
t = t_min:dt:t_max;
% Space discretization is the same than before
% Matrix to save the free response at different times
f_r_ND_NF = zeros(n_t,n_x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF THE FREE RESPONSE FOR 2D PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         for i = 1:n_t
             % Components of the free response at the current time
             w_1_ND_NF = cos(omega_1*t(i));
             w_2_ND_NF = cos(omega_2*t(i));
             % The results for t = (i-1)*dt are stored in the i-th row 
             f_r_ND_NF(i,:) = IC_1*w_1_ND_NF*eigenfunction_1_2D + IC_2*w_2_ND_NF*eigenfunction_2_2D;
         end;

         % 2D plot
         figure(6)
         plot(x_2D,f_r_ND_NF(1,:),'-m',x_2D,f_r_ND_NF(2,:),'-r',x_2D,f_r_ND_NF(3,:),'-g',x_2D,f_r_ND_NF(4,:),'-c',x_2D,f_r_ND_NF(5,:),'-b',x_2D,f_r_ND_NF(6,:),'-y')
         legend('t_{0}','t_{1}','t_{2}','t_{3}','t_{4}','t_{5}')
         grid on
         xlabel('x/L')
         ylabel('w (x)')
         title('Undamped free response in space at different times with no flow ( \Lambda = 0 )')
