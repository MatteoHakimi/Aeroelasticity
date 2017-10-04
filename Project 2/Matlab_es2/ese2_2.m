%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXERCISE 2.2 - AEROELASTIC FREE RESPONSE WITH QUASI-STEADY AERODYNAMICS
%                USING THE EIGENVALUES AND EIGENVECTORS OF THE STATE 
%                MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The present code computes the free response of the typical section with
%  quasi-steady aerodynamics using the eigenvalues and eigenvectors of the 
%  state matrix of the system.
%  The eigenvalues, eigenvectors, and the constants depending on the 
%  applied initial conditions are first computed for a given value of the 
%  parameter U (below or above the flutter speed). 
%  The solution is evaluated and presented in 2D and 3D plots. 
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
mass_ratio =5 ;
% Frequency ratio
freq_ratio =0.5 ;
% Non-dimensional position of the elastic center
xi_E =0.3 ;
% Non-dimensional position of the mass center
xi_G =0.45 ;
% Non-dimensional radius of gyration
r_alpha_2 =0.25 ;
% Parameter for structural damping
damp =0.005 ;

%%%%%%%%%%%%%%%%%%%%%%
% INITIAL CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%

% Initial condition x1(0)
IC_1 =0.01 ;
% Initial condition x2(0)
IC_2 =0 ;
% Initial condition x3(0)
IC_3 =0 ;
% Initial condition x4(0)
IC_4 =0 ;

% 2.2.1 SYSTEM MATRICES WITH QUASI-STEADY AERODYNAMICS

       % Non-dimensional mass matrix
       M = [1 xi_G;xi_G r_alpha_2-xi_E^2+2*xi_E*xi_G];
       % Non-dimensional stiffness matrix
       K = freq_ratio^2*[1 xi_E;xi_E xi_E^2+r_alpha_2/(freq_ratio^2)];
       % Non-dimensional structural damping matrix
       D = freq_ratio*2*damp*[1 0;0 1];

       % Non-dimensional matrix due to the steady aerodynamic loads
       SA = [0 2/mass_ratio;0 0];

       % Insert the elements of the non-dimensional matrix due to the quasi-steady aerodynamic loads
       QSA = [2/mass_ratio 0; 0 0];

% 2.2.2 STATE MATRIX A 

% Insert the flutter speed estimated from the root locus (quasi-steady aerodynamics)
U_QSA_flutter =0.6618 ;
% The free response will be computed for U = U_QSA_flutter*perc
perc =0.8 ;
% Current value of U 
U = perc*U_QSA_flutter;
% Order of the differential system (written in second-order form)
n =2 ;
% Identity matrix ( n*n ) 
I = eye(n);

       % State matrix
       A = [zeros(n,n) I;-M^(-1)*(K+SA*U^2) -M^(-1)*(D+QSA*U)]

% 2.2.2 QUANTITIES TO COMPUTE THE AEROELASTIC RESPONSE

       % Eigenvalues of the state matrix
       eigenvalues = eig(A)
     
       % Eigenvectors of the state 
       [u,diag_A] = eig(A);
       u

       % Constants depending on the initial conditions
       const = u^(-1)*[IC_1;IC_2;IC_3;IC_4]

% 2.2.3 AEROELASTIC FREE RESPONSE

% Minimum value of time ( non-dimensional )
t_min =0 ;
% Maximum value of time ( non-dimensional )
t_max =100 ;
% Delta t ( non-dimensional ) 
dt =0.01 ;
% Vector of times
t = t_min:dt:t_max;
% Number of times
nt = length(t);
% Matrix to store the solution
free_response = zeros(2*n,nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME-HISTORIES OF THE STATE VARIABLES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for i=1:2*n
              % Computation of the time-history of the i-th state variable
              % 1 = vertical traslation of the quarter-chord point
              % 2 = angle of attack
              % 3 = vertical velocity of the quarter-chord point
              % 4 = pitch rate 
              sum = 0;
              for j=1:2*n
                  sum = sum + const(j)*u(i,j)*exp(eigenvalues(j)*t);
              end
              % The time-history for the i-th state variable is stored in the i-th row of the matrix 
              free_response(i,:) = sum;
        end
               
        % Non-dimensional vertical translation of the quarter-chord point and angle of attack
        figure(1)
        plot(t,free_response(1,:),'b-',t,free_response(2,:),'r-')
        grid on 
        legend('Vertical translation','Angle of attack')
        xlabel('t [-]')
        ylabel('Vertical translation and angle of attack [-]')
        title('Time-histories of the vertical translation and angle of attack')
        
        % Non-dimensional vertical velocity of the quarter-chord point and pitch rate
        figure(2)
        plot(t,free_response(3,:),'b-',t,free_response(4,:),'r-')
        grid on 
        legend('Vertical velocity','Pitch rate')
        xlabel('t [-]')
        ylabel('Vertical velocity and pitch rate [-]')
        title('Time-histories of the vertical velocity and pitch rate')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERTICAL DISPLACEMENT FIELD FOR 3D PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Minimum value of x (leading edge) ( non-dimensional )
x_min = -1/4;
% Maximum value of x (trailing edge) ( non-dimensional )
x_max = 3/4;
% Space discretization ( non-dimensional ) 
dx = 0.01;
% Minimum value of time ( non-dimensional )
t_min_3D =0 ;
% Maximum value of time ( non-dimensional ) 
t_max_3D =100 ;
% Delta t ( non-dimensional ) 
dt_3D =0.01 ;
% Discretization of space and time for 3D plot
[x,t] = meshgrid(x_min:dx:x_max,t_min_3D:dt_3D:t_max_3D);

          h = 0;
          alpha = 0;      
          % Computation of h and alpha for 3D plot
          for i=1:2*n
             h = h + const(i)*u(1,i)*exp(eigenvalues(i)*t);
             alpha = alpha + const(i)*u(2,i)*exp(eigenvalues(i)*t);
          end
          % Non-dimensional vertical displacement field for 3D plot
          w = h+alpha.*x;

          % 3D plot 
          figure(3)
          mesh(x,t,real(w))
          xlabel('x [-]')
          ylabel('t [-]')
          zlabel('w (x,t) [-]')
          grid on
          title('Vertical displacement w(x,t) in space and time')
       