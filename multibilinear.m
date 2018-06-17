clear all
close all
clc 
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
%% initial setting
loop=100;
num=100; % number of moment to choose, need to be even number
%% example 1
% n=7;
% v=2;
% A=[-0.81 0.47 -0.43 1.6 0.26 -0.4 0.92; -0.61 -1.9 0.8 -1.6 2 0.98 -0.9; 0.5 -1.2 -2.1 -1.6 -1.1 0.14 -0.87; -1.3 2.1 0.47 -1.2 3.7 -1.2 -1.3; -0.24 -0.081 1.6 -3.6 -1.3 1.7 -2.6; 1.3 -0.96 -1.3 -0.57 -2.4 -2.4 -0.36; -0.16 15 -0.99 1.5 0.61 -2.2 -3.3];
% B=[0; 0;-0.196;1.42;0.292;0.198;1.59];
% C=[-0.804 0 0.835 -0.244 0.216 -1.17 -1.15];
% N=diag(-1*ones(1,n));

%% example 2 random generated
% n=100;
% v=4;
% eps=0.1;
% sys=rss(n);
% A=sys.a-eps*eye(n);
% B=sys.b;
% C=sys.c;
% N=rand(n)/100;
%% example 3
n=9;
v=4;
A=[-0.812 0.472 -0.433 1.58 0.258 -0.404 0.923 0.1 0; -0.612 -1.89 0.804 -1.6 2.01 0.983 -0.904 0.2 0;0.502 -1.18 -2.12 -1.57 -1.14 0.14 -0.868 0 0.4;-1.28 2.12 0.472 -1.16 3.65 -1.23 -1.26 0 0 ; -0.242 -0.0812 1.58 -3.58 -1.26 1.72 -2.59 0 0;1.29 -0.959 -1.28 -0.571 -2.43 -2.44 -0.365 0 0;-0.156 1.53 -0.992 1.47 0.611 -2.23 -3.28 0 0; 0 0 0 0 0 0 0 -0.2 0; 0 0 0 0 0 0 0 0 -1.3];
B=[0.2; 0 ;-0.31;1.52;0.292;0.198;1.59;0.01;0];
C=[-0.9 0 0.84 -0.33 0.216 -1.17 -1.15 0 0];
N=[-0.448 -0.329 -0.298 -0.028 -0.028 -0.0665 -0.028 -0.028 -0.028; -0.209 -0.28 -0.266 -0.028 -0.028 -0.014 -0.028 -0.028 -0.028; -0.19 -0.171 -0.308 -0.102 -0.028 -0.028 -0.0595 -0.028 -0.028; -0.028 -0.028 -0.15 -0.308 -0.028 -0.028 -0.0805 -0.028 -0.0245; -0.028 -0.028 -0.028 -0.028 -0.308 -0.028 -0.028 -0.028 -0.028; -0.0511 -0.0196 -0.028 -0.028 -0.028 -0.308 -0.0427 -0.028 -0.028;-0.028 -0.028 -0.0805 -0.115 -0.028 -0.0525 -0.308 -0.028 -0.028;-0.028 -0.028 -0.028 -0.028 -0.028 -0.028 -0.028 -0.308 -0.028;-0.028 -0.028 -0.028 -0.0259 -0.028 -0.028 -0.028 -0.028 -0.308];

%% example 4
% n=5;
% v=2;
% A=[0 0 0.024 0 0; 1 0 -0.26 0 0; 0 1 0.9 0 0; 0 0 0.2 0 -0.06; 0 0 0.15 1 0.5];
% B=[0.8;0.6;0.4;0.2;0.5];
% C=[0.2 0.4 0.6 0.8 1];
% N=diag([0.1,0.2,0.3,0.4,0.5]);
%% example 5 1006th linear system with random N
% n=1006;
% v=2;
% %select A, B ,C  , S
% A_1=[0 100; -100 0];
% A_2=2*A_1;
% A_3=3*A_1;
% gamma=(-1:-1:-1000);
% A_1000=diag(gamma);
% A=blkdiag(A_1000,A_1,A_2,A_3);
% B=[10*ones(6,1);ones(1000,1)];
% C=B';
% N=rand(n)/100;
%% signal and input

%L norm <1
%L=[0.15 0 0.15 0];
L=repmat([1 0],1,v/2)*0.15;

%signal generator
S_1=[0 1; -1 0];
S=blkdiag(0.15*S_1,30*S_1);

x0=[1 ; 1 ; 1; 1 ;1 ;1 ; 1 ; 1 ; 1];
%x0=rand(n,1);

%w0=random('Normal',0,0.5,v,1);
w0=[0.254;-0.431;0.256;-0.138];

z0=[w0;x0];

%% system
time=80;
[t,z]=ode45(@(t,z) bifunc(z,A,B,N,S,L,v),[0,time],z0,opts);

w=z(:,1:v);
x=z(:,v+1:end);

y=C*x';
%% PI
PI=zeros(n,loop);
PIiodd=zeros(n,v,loop+1);
PIieven=zeros(n,1,loop);
% PIi(n,v,1)=zeros(n,1);

SylA=A;
SylS= -1 * S;
SylBL= -1 *B * L;
PIiodd(:,:,1)=sylvester(SylA,SylS,SylBL);
%PI(:,1)=C*PIi(:,:,1)*w0;
for i=2:2:loop

    BigA=-kron(L',A);
    PIieven(:,:,i)=(BigA'*BigA)\BigA'*reshape(N*PIiodd(:,:,i-1)*norm(L'*L),[v*n,1]);
    
    %PI(:,i)=PI(:,i-1)+C*PIieven(:,:,i)*w0^i;

    SylA=A;
    SylS= S;
    SylBL= -1 * N * PIieven(:,:,i) * L;
    PIiodd(:,:,i+1)=sylvester(SylA,SylS,SylBL);
    %PI(:,i+1)=PI(:,i)+C*PIiodd(:,:,i+1)*w0^(i+1);
end

%% Moments

for j=1:size(t)
    omega=w(j,:)';
for i=2:2:num
   
    PI(:,i)=PI(:,i-1)+PIiodd(:,:,i-1) * omega;
    
    omega=w(j,:)*omega;
    
    PI(:,i+1)=PI(:,i)+PIieven(:,:,i)*omega;
    
    omega=w(j,:)'*omega;
end
%for i=1:1:loop+1
hpi(j,1)=C*PI(:,end);
end

%% plot graph orginal system
figure(1)
hold on
plot(t,hpi);
title('Moment of original system')
hold off
print('1','-depsc');
%plot y vs t and moments vs t
figure(2)
hold on
plot(t,y);
plot(t,hpi);
title('steady state and moment of original system');
xlabel('time');
ylabel('output');
legend('y','moment')
hold off
print('2','-depsc');

%% reduced system

x_r0=random('Normal',0,0.5,v,1);
%x_r0=zeros(v,1);
%x_r0=[1;1]
z0=[w0;x_r0];
%% design of G and M
eigenval=(-2:-1:-v-1);
G=place(S',L',eigenval);
G=G';

M=eye(v);

%% simulation
%time=100;
[t1,z]=ode45(@(t,z) reducedbifunc(z,G,M,S,L,v),[0,time],z0,opts);

w_r=z(:,1:v);
x_r=z(:,v+1:end);

%% steady state
for j=1:size(t1)
    omega=x_r(j,:)';
for i=2:2:num
   
    PI(:,i)=PI(:,i-1)+PIiodd(:,:,i-1) * omega;
    
    omega=x_r(j,:)*omega;
    
    PI(:,i+1)=PI(:,i)+PIieven(:,:,i)*omega;
    
    omega=x_r(j,:)'*omega;
end
%for i=1:1:loop+1
hpi_r(j,1)=C*PI(:,end);
end

y_r=hpi_r;

%% plot graph reduced system
figure(3)
hold on
plot(t1,hpi_r);
title('Steady state of reduced system')
hold off
print('3','-depsc');
%plot y vs t and moments vs t
figure(4)
hold on
plot(t,y);
plot(t1,y_r);
title('steady state(original system) vs steady state(reduced system)');
xlabel('time');
ylabel('output');
legend('y','y_r')
hold off
print('4','-depsc');



%% linearised system bode plot

originalsys=ss(A,B,C,0);

F=S-G*L;
H=C*PIiodd(:,:,1);

reducedsys=ss(F,G,H,0);

figure(5)
hold on
bode(originalsys);
bode(reducedsys,'r');
title(sprintf('Bode plot v= %d',v));
hold off
print('5','-depsc');







