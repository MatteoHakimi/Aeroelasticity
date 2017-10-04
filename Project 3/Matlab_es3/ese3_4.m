%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    EXERICSE 3.4 - AEROELASTIC FREE RESPONSE WITH UNSTEADY AERODYNAMICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The present code computes the aeroelastic free response of the typical
%  section in the case of unsteady aerodynamics using the eigenvalues
%  and eigenvectors of the state matrix ( finite-state aerodynamics ).
%  The system of ODEs is obtained using Padé approximation to C(k).
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
mass_ratio = 5;
% Frequency ratio
freq_ratio =0.5;
% Non-dimensional position of the elastic center
xi_E =0.3 ;
% Non-dimensional position of the mass center
xi_G =0.45 ;
% Non-dimensional radius of gyration
r_alpha_2 =0.25 ;

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
% Initial condition x5(0)
IC_5 =0 ;
% Initial condition x6(0)
IC_6 =0 ;


% 3.4.1 STATE MATRIX WITH ADDED AERODYNAMIC STATES  

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
coef_D_0 = 0.044;
coef_D_1 =0.552 ;
coef_D_2 =1 ;

% Insert the non-dimensional flutter velocity
U_F =1.17 ;
% The free response is evaluated for U = perc*U_F 
perc =1.2 ;
% Current non-dimensional velocity
U = perc*U_F;

         % New non-dimensional mass matrix (3*3)
         M = U^2*[1+1/mass_ratio xi_G+1/(2*mass_ratio) 0;xi_G+1/(2*mass_ratio) r_alpha_2-xi_E^2+2*xi_E*xi_G+3/(8*mass_ratio) 0;-2*coef_N_1/(mass_ratio*U^2) -2*coef_N_1/(mass_ratio*U^2) 1/U^2];
         % New non-dimensional damping matrix (3*3)
         D = U^2*[2*C_0 1+2*C_0 0;0 1 0;-2*coef_N_0/U^2 -2*(coef_N_1+coef_N_0)/U^2 coef_D_1*mass_ratio/U^2]/mass_ratio;
         % New non-dimensional stiffness matrix (3*3)
         K = freq_ratio^2*[1 xi_E+2*C_0*U^2/(mass_ratio*freq_ratio^2) U^2/freq_ratio^2;xi_E xi_E^2+r_alpha_2/freq_ratio^2 0;0 -2*coef_N_0/(mass_ratio*freq_ratio^2) coef_D_0/freq_ratio^2];           
 
% Order of the differential system (written in second-order form)
n = 3;
% Identity matrix (n*n)
I = eye(n);

         % State matrix (6*6)
         A = [zeros(n,n) I;-M^(-1)*K -M^(-1)*D]
       
         
% 3.4.2 QUANTITIES TO OBTAIN THE AEROELASTIC FREE RESPONSE

         % Eigenvalues of the state matrix
         eigenvalues(1:2*n)=eig(A)
         
         % Eigenvectors of the state matrix
         [u,diag_A] = eig(A);
         u

         % Constants depending on the initial conditions
         const = u^(-1)*[IC_1; IC_2;IC_3;IC_4;IC_5;IC_6]

% 3.4.3 AEROELASTIC FREE RESPONSE 

% Minimum value of time ( non-dimensional )
t_min =0 ;
% Maximum value of time ( non-dimensional )
t_max =100 ;
% Delta t ( non-dimensional ) 
dt =0.001 ;
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
              % 3 = added aerodynamic state
              % 4 = vertical velocity of the quarter-chord point
              % 5 = pitch rate
              % 6 = time derivative of the added aerodynamic state
              sum = 0;
              for j=1:2*n
                  sum = sum + const(j)*u(i,j)*exp(eigenvalues(j)*t);
              end;
              
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
        plot(t,free_response(4,:),'b-',t,free_response(5,:),'r-')
        grid on 
        legend('Vertical velocity','Pitch rate')
        xlabel('t [-]')
        ylabel('Vertical velocity and pitch rate [-]')
        title('Time-histories of the vertical velocity and pitch rate')       
    
        % Added aerodynamic state and its time-derivative
        figure(3)
        plot(t,free_response(3,:),'b-',t,free_response(6,:),'r-')
        grid on 
        legend('Aerodynamic state','Time-derivative of the aerodynamic state')
        xlabel('t [-]')
        ylabel('Aerodynamic state and its time-derivative [-]')
        title('Time-histories of the aerodynamic state and its time-derivative')         