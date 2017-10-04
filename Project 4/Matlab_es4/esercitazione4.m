clear all
clc
close all

P1=textread('polo1.txt');
P2=textread('polo2.txt');
P3=textread('polo3.txt');

real_1=P1(:,6);
imag_1=P1(:,7);
real_2=P2(:,6);
imag_2=P2(:,7);
real_3=P3(:,6);
imag_3=P3(:,7);

vel_1=P1(:,3);
damp_1=P1(:,4);
vel_2=P2(:,3);
damp_2=P2(:,4);
vel_3=P3(:,3);
damp_3=P3(:,4);




figure(1)
hold on
plot(real_1,imag_1,'b.')
plot(real_2,imag_2,'g.')
plot(real_3,imag_3,'r.')
plot(real_1(1,1),imag_1(1,1),'kx')
plot(real_2(1,1),imag_2(1,1),'kx')
plot(real_3(1,1),imag_3(1,1),'kx')
xlabel('Re')
ylabel('Im')
grid on
title('Root Locus')

figure(2)
plot(vel_1,damp_1,'b')
hold on
plot(vel_2,damp_2,'g.')
plot(vel_3,damp_3,'r')
xlabel('Velocity [m/s]')
ylabel('Damping [-]')
grid on
title('Damping vs velocity')



U_D1=P1(find(P1(:,6)>0)-1,3);
U_D=min(U_D1)

U_F1=P3(find(P3(:,6)>0),3);
U_F=min(U_F1)
k_F1=P3(find(P3(:,6)>0),1);
k_F=k_F1(1,1)