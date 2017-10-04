%% Luogo delle radici della risposta libera.
clear 
clc


f_1=15; % frequenza modo 1 [Hz]
f_2=30; % frequenza modo 2 [Hz]

zita_1=0.03; % smorzamento zita 1
zita_2=0.03; % smorzamento zita 2

omega_1=2*pi*f_1; % pulsazione modo 1
omega_2=2*pi*f_2; % pulsazione modo 2

a4=1;                                 % coefficiente a4 del polinomio caratteristico
a3=(zita_1+zita_2);                   % coefficiente a3 del polinomio caratteristico
a2=omega_1^2+omega_2^2+zita_1*zita_2; % coefficiente a2 del polinomio caratteristico
a1=zita_2*omega_1^2+zita_1*omega_2^2; % coefficiente a1 del polinomio caratteristico
a0=(omega_1)^2*(omega_2)^2;           % coefficiente a0 del polinomio caratteristico

num=1;              % numeratore funzione di trasferimento
den=[a4,a3,a2,a1,a0]; % denominatore funzione di trasferimento

Hs=tf(num,den); %  funzione di trasferimento H(s)
rlocus(Hs); % luogo delle radici
axis([-0.2 0.2 200 -200])
