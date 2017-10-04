clear all
clc
close all

%M=0.5566 rho=0.9046 t=1.02
P1=textread('f_1.txt');
P2=textread('f_2.txt');
P3=textread('f_3.txt');
P4=textread('f_4.txt');
P5=textread('f_5.txt');
P6=textread('f_6.txt');

%M=0.0001 rho=0.9046 t=1.02
P1M=textread('f_1M.txt');
P2M=textread('f_2M.txt');
P3M=textread('f_3M.txt');
P4M=textread('f_4M.txt');
P5M=textread('f_5M.txt');
P6M=textread('f_6M.txt');

%M=0.5566 rho=0.19475 t=1.02
P1R=textread('f_1R.txt');
P2R=textread('f_2R.txt');
P3R=textread('f_3R.txt');
P4R=textread('f_4R.txt');
P5R=textread('f_5R.txt');
P6R=textread('f_6R.txt');

%M=0.5566 rho=0.9046 t=1.52
P1t=textread('f_1t.txt');
P2t=textread('f_2t.txt');
P3t=textread('f_3t.txt');
P4t=textread('f_4t.txt');
P5t=textread('f_5t.txt');
P6t=textread('f_6t.txt');



real_1=P1(:,6);
imag_1=P1(:,7);
real_2=P2(:,6);
imag_2=P2(:,7);
real_3=P3(:,6);
imag_3=P3(:,7);
real_4=P4(:,6);
imag_4=P4(:,7);
real_5=P5(:,6);
imag_5=P5(:,7);
real_6=P6(:,6);
imag_6=P6(:,7);

real_1M=P1M(:,6);
imag_1M=P1M(:,7);
real_2M=P2M(:,6);
imag_2M=P2M(:,7);
real_3M=P3M(:,6);
imag_3M=P3M(:,7);
real_4M=P4M(:,6);
imag_4M=P4M(:,7);
real_5M=P5M(:,6);
imag_5M=P5M(:,7);
real_6M=P6M(:,6);
imag_6M=P6M(:,7);

real_1R=P1R(:,6);
imag_1R=P1R(:,7);
real_2R=P2R(:,6);
imag_2R=P2R(:,7);
real_3R=P3R(:,6);
imag_3R=P3R(:,7);
real_4R=P4R(:,6);
imag_4R=P4R(:,7);
real_5R=P5R(:,6);
imag_5R=P5R(:,7);
real_6R=P6R(:,6);
imag_6R=P6R(:,7);

real_1t=P1t(:,6);
imag_1t=P1t(:,7);
real_2t=P2t(:,6);
imag_2t=P2t(:,7);
real_3t=P3t(:,6);
imag_3t=P3t(:,7);
real_4t=P4t(:,6);
imag_4t=P4t(:,7);
real_5t=P5t(:,6);
imag_5t=P5t(:,7);
real_6t=P6t(:,6);
imag_6t=P6t(:,7);



%plot M=0.5566 rho=0.9046 t=1.02
figure(1)
hold on
plot(real_1,imag_1,'b.')
plot(real_2,imag_2,'g.')
plot(real_3,imag_3,'r.')
plot(real_4,imag_4,'c.')
plot(real_5,imag_5,'m.')
plot(real_6,imag_6,'y.')
plot(real_1(1,1),imag_1(1,1),'kx')
plot(real_2(1,1),imag_2(1,1),'kx')
plot(real_3(1,1),imag_3(1,1),'kx')
plot(real_4(1,1),imag_4(1,1),'kx')
plot(real_5(1,1),imag_5(1,1),'kx')
plot(real_6(1,1),imag_6(1,1),'kx')
xlabel('Re')
ylabel('Im')
grid on
title('Root Locus')

%plot M=0.0001 rho=0.9046 t=1.02
figure(2)
hold on
plot(real_1M,imag_1M,'b.')
plot(real_2M,imag_2M,'g.')
plot(real_3M,imag_3M,'r.')
plot(real_4M,imag_4M,'c.')
plot(real_5M,imag_5M,'m.')
plot(real_6M,imag_6M,'y.')
plot(real_1M(1,1),imag_1M(1,1),'kx')
plot(real_2M(1,1),imag_2M(1,1),'kx')
plot(real_3M(1,1),imag_3M(1,1),'kx')
plot(real_4M(1,1),imag_4M(1,1),'kx')
plot(real_5M(1,1),imag_5M(1,1),'kx')
plot(real_6M(1,1),imag_6M(1,1),'kx')
xlabel('Re')
ylabel('Im')
grid on
title('Root Locus')

%plot M=0.5566 rho=0.19475 t=1.02
figure(3)
hold on
plot(real_1R,imag_1R,'b.')
plot(real_2R,imag_2R,'g.')
plot(real_3R,imag_3R,'r.')
plot(real_4R,imag_4R,'c.')
plot(real_5R,imag_5R,'m.')
plot(real_6R,imag_6R,'y.')
plot(real_1R(1,1),imag_1R(1,1),'kx')
plot(real_2R(1,1),imag_2R(1,1),'kx')
plot(real_3R(1,1),imag_3R(1,1),'kx')
plot(real_4R(1,1),imag_4R(1,1),'kx')
plot(real_5R(1,1),imag_5R(1,1),'kx')
plot(real_6R(1,1),imag_6R(1,1),'kx')
xlabel('Re')
ylabel('Im')
grid on
title('Root Locus')

%plot M=0.5566 rho=0.9046 vs rho=0.19475 t=1.02
figure(4)
hold on
plot(real_1,imag_1,'b.')
plot(real_1R,imag_1R,'r.')
plot(real_2,imag_2,'b.')
plot(real_2R,imag_2R,'r.')
plot(real_3,imag_3,'b.')
plot(real_3R,imag_3R,'r.')
plot(real_4,imag_4,'b.')
plot(real_4R,imag_4R,'r.')
plot(real_5,imag_5,'b.')
plot(real_5R,imag_5R,'r.')
plot(real_6,imag_6,'b.')
plot(real_6R,imag_6R,'r.')
plot(real_1R(1,1),imag_1R(1,1),'kx')
plot(real_2R(1,1),imag_2R(1,1),'kx')
plot(real_3R(1,1),imag_3R(1,1),'kx')
plot(real_4R(1,1),imag_4R(1,1),'kx')
plot(real_5R(1,1),imag_5R(1,1),'kx')
plot(real_6R(1,1),imag_6R(1,1),'kx')
plot(real_1(1,1),imag_1(1,1),'kx')
plot(real_2(1,1),imag_2(1,1),'kx')
plot(real_3(1,1),imag_3(1,1),'kx')
plot(real_4(1,1),imag_4(1,1),'kx')
plot(real_5(1,1),imag_5(1,1),'kx')
plot(real_6(1,1),imag_6(1,1),'kx')
legend('\rho_{\infty}=0.9046 kg/m^3','\rho_{\infty}=0.19475 kg/m^3')
xlabel('Re')
ylabel('Im')
grid on
title('Root Locus')

%plot M=0.5566 rho=0.9046 t=1.52
figure(5)
hold on
plot(real_1t,imag_1t,'b.')
plot(real_2t,imag_2t,'g.')
plot(real_3t,imag_3t,'r.')
plot(real_4t,imag_4t,'c.')
plot(real_5t,imag_5t,'m.')
plot(real_6t,imag_6t,'y.')
plot(real_1t(1,1),imag_1t(1,1),'kx')
plot(real_2t(1,1),imag_2t(1,1),'kx')
plot(real_3t(1,1),imag_3t(1,1),'kx')
plot(real_4t(1,1),imag_4t(1,1),'kx')
plot(real_5t(1,1),imag_5t(1,1),'kx')
plot(real_6t(1,1),imag_6t(1,1),'kx')
xlabel('Re')
ylabel('Im')
grid on
title('Root Locus')

%plot M=0.5566 rho=0.9046 t=1.02 vs t=1.52
figure(6)
hold on
plot(real_1,imag_1,'b.')
plot(real_1t,imag_1t,'r.')
plot(real_2,imag_2,'b.')
plot(real_2t,imag_2t,'r.')
plot(real_3,imag_3,'b.')
plot(real_3t,imag_3t,'r.')
plot(real_4,imag_4,'b.')
plot(real_4t,imag_4t,'r.')
plot(real_5,imag_5,'b.')
plot(real_5t,imag_5t,'r.')
plot(real_6,imag_6,'b.')
plot(real_6t,imag_6t,'r.')
plot(real_1t(1,1),imag_1t(1,1),'kx')
plot(real_2t(1,1),imag_2t(1,1),'kx')
plot(real_3t(1,1),imag_3t(1,1),'kx')
plot(real_4t(1,1),imag_4t(1,1),'kx')
plot(real_5t(1,1),imag_5t(1,1),'kx')
plot(real_6t(1,1),imag_6t(1,1),'kx')
plot(real_1(1,1),imag_1(1,1),'kx')
plot(real_2(1,1),imag_2(1,1),'kx')
plot(real_3(1,1),imag_3(1,1),'kx')
plot(real_4(1,1),imag_4(1,1),'kx')
plot(real_5(1,1),imag_5(1,1),'kx')
plot(real_6(1,1),imag_6(1,1),'kx')
legend('\tau_{p}=1.02 mm','\tau_{p}=1.52 mm')
xlabel('Re')
ylabel('Im')
grid on
title('Root Locus')




U_D1=P1(find(P1(:,6)>0)-1,3);
U_D=min(U_D1)
 
U_F1=P2(find(P2(:,6)>0),3);
U_F=min(U_F1)
k_F1=P2(find(P2(:,6)>0),1);
k_F=k_F1(1,1)

U_D1M=P1M(find(P1M(:,6)>0)-1,3);
U_DM=min(U_D1M)
 
U_F1M=P2M(find(P2M(:,6)>0),3);
U_FM=min(U_F1M)
k_F1M=P2M(find(P2M(:,6)>0),1);
k_FM=k_F1M(1,1)


U_F1t=P2t(find(P2t(:,6)>0),3);
U_Ft=min(U_F1t)
k_F1t=P2t(find(P2t(:,6)>0),1);
k_Ft=k_F1t(1,1)

