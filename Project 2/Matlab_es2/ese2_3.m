%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXERCISE 2.3 - AEROELASTIC FREE RESPONSE WITH QUASI-STEADY AERODYNAMICS
%                USING NUMERICAL METHODS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The present code computes the free response of the typical section with
%  quasi-steady aerodynamics using numerical methods. The results are 
%  presented in terms of the time-histories of the state variables. 
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
mass_ratio = ;
% Frequency ratio
freq_ratio = ;
% Non-dimensional position of the elastic center
xi_E = ;
% Non-dimensional position of the mass center
xi_G = ;
% Non-dimensional radius of gyration
r_alpha_2 = ;
% Parameter for structural damping
damp = ;

%%%%%%%%%%%%%%%%%%%%%%
% INITIAL CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%

% Initial condition x1(0)
IC_1 = ;
% Initial condition x2(0)
IC_2 = ;
% Initial condition x3(0)
IC_3 = ;
% Initial condition x4(0)
IC_4 = ;

% 2.3.1 STATE MATRIX A 

       % Non-dimensional mass matrix
       M = [1 xi_G;xi_G r_alpha_2-xi_E^2+2*xi_E*xi_G];
       % Non-dimensional stiffness matrix
       K = freq_ratio^2*[1 xi_E;xi_E xi_E^2+r_alpha_2/(freq_ratio^2)];
       % Non-dimensional damping matrix
       D = freq_ratio*2*damp*[1 0;0 1];

       % Non-dimensional matrix due to the steady aerodynamic loads
       SA = [0 2/mass_ratio;0 0];

       % Insert the elements of the non-dimensional matrix due to the quasi-steady aerodynamic loads
       QSA = [];
       
% Insert the flutter speed estimated from the root locus (quasi-steady aerodynamics)
U_QSA_flutter = ;
% The free response will be computed for U = U_QSA_flutter*perc 
perc = ;
% Current value of U 
U = perc*U_QSA_flutter;
% Order of the differential system (written in second-order form)
n = ;
% Identity matrix ( n*n ) 
I = eye(2*n);
       
       % State matrix
       A = [zeros(n,n) eye(n);-M^(-1)*(K+SA*U^2) -M^(-1)*(D+QSA*U)];
       
     
% 2.3.2 NUMERICAL METHODS
  
% Choose a numerical method
% alpha = 1 for an esplicit method
% alpha = 0.5 for Crank-Nicholson method 
% alpha = 0 for an implicit method
alpha = ;
% Minimum value of time ( non-dimensional )
t_min = ;
% Maximum value of time ( non-dimensional )
t_max = ;
% Delta t ( non-dimensional ) 
dt = ;
% Time discretization
t = t_min:dt:t_max;
% Number of times
nt = length(t);
% Matrix to store the solution
free_response_NM = zeros(2*n,nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF THE FREE RESPONSE FOR 2D PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         % Initial time
         free_response_NM(:,1)=[IC_1;IC_2;IC_3;IC_4];
         % Integration in time
         for j=1:nt-1
            % The i-th state variable is stored in the i-th row
            % 1 = vertical translation of the quarter-chord point 
            % 2 = angle of attack
            % 3 = vertical velocity of the quarter-chord point
            % 4 = pitch rate 
            free_response_NM(:,j+1)=(I-(1-alpha)*dt*A)^(-1)*(I+alpha*dt*A)*free_response_NM(:,j);
         end


        % Non-dimensional vertical translation of the quarter-chord point and angle of attack
        figure(1)
        plot(t,free_response_NM(1,:),'b-',t,free_response_NM(2,:),'r-')
        grid on 
        legend('Vertical translation','Angle of attack')
        xlabel('t [-]')
        ylabel('Vertical translation and angle of attack [-]')
        title('Time-histories of the vertical translation and angle of attack')
        
        % Non-dimensional vertical velocity of the quarter-chord point and pitch rate
        figure(2)
        plot(t,free_response_NM(3,:),'b-',t,free_response_NM(4,:),'r-')
        grid on 
        legend('Vertical velocity','Pitch rate')
        xlabel('t [-]')
        ylabel('Vertical velocity and pitch rate [-]')
        title('Time-histories of the vertical velocity and pitch rate')