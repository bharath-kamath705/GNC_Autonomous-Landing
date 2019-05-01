clc
clear all
close all

% %aircraft data:altitude, moment of inertia etc
 data=xlsread('boeing747_data.xlsx');

% %Dimensional derivatives case 1 Mach 0.2 
 dd=xlsread('dimensional_derivatives_case1');

%reference conditions [g theta_ref u_ref]
ref=[32.2 0 50]; %u_ref data(3,1)is the velocity mentioned in data file

%Aircraft system X_dot=AX+Bu
%states are {del_u w q del_theta]'
[A1,B1]=long_model(data(:,1),dd,ref);

%landing parameters
gsa=3*pi/180; %glide slope angle
s_ref=30000;
h_ref=1500;
u_ref=200 %data(3,1);



%adaptive r(t)
del_ws=10;
C=[tan(gsa) -1 0 0; 0 1 0 0];
A1_app=[A1 zeros(4,2); -C zeros(2,2)];
B1_app =[B1; zeros(2,2)];
D1=[zeros(4,1);-u_ref*tan(gsa); del_ws]; %-u_ref*tan(gsa)

B1_app(2,2)=0.2;
disp('controlability matrix rank')
P = [B1_app A1_app*B1_app A1_app^2*B1_app A1_app^3*B1_app A1_app^4*B1_app A1_app^5*B1_app ];
rank(P)


%eigen structure
%Two desired eigenvalues (dutch roll) are determined by desired damping ratio and
%natural frequency zeta and omega_n respectively
zeta1=0.1; 
omega_n1=3; 
lambda_d(1) = -2+j*0.5;
lambda_d(2) = conj(lambda_d(1));
lambda_d(3) = -0.05+ j*0;
lambda_d(4) = -2.9;
lambda_d(5) = -3.5;
lambda_d(6) = -4.5;



% Extracting desired eigenvectors directly from the null-space. 
for i=1:6
mat(:,:,i) = [(lambda_d(i)*eye(6)-A1_app) B1_app];
nullspace(:,:,i) = null(mat(:,:,i),'r');
vu(:,i) = 0.2*i*nullspace(:,1,i)+0.5*nullspace(:,2,i);
V(:,i)=vu(1:6,i);
U(:,i)=vu(7:8,i);
end

K = U*inv(V); %returns K as complex variable but with 0 imaginary part
K=real(K);



% x0=[0 0 0 0 0 0 s_ref h_ref];
% [t,x] = ode45('gsa_land_R',[0:0.02:100],x0,[],A1_app,B1_app,K,D1,u_ref);
% plot(t,x(:,8))


tol=0.1
ti=0
del_t=1
t_tot=185
n=t_tot/del_t
tf=del_t
X=zeros(1,8)
T=0
x0=[0 0 0 0 0 0 s_ref h_ref];
for i=1:n
    
    [t,x] = ode45('gsa_land_R',[ti tf],x0,[],A1_app,B1_app,K,D1,u_ref);
    
    if x(end,8)>(x(end,7)*tan(gsa)+tol)
        D1(end)=D1(end)-0.5;
    elseif x(end,8)<(x(end,7)*tan(gsa)-tol)
        D1(end)=D1(end)+0.5;
    else
        ;
    end
    
    X=[X;x];
    T=[T;t];
    ti=tf
    tf=tf+del_t;
    x0=x(end,:);
end

plot(T(2:end),X(2:end,2))

figure
plot(X(2:end,7),X(2:end,8))
hold on
plot(X(2:end,7),tan(gsa)*X(2:end,7),'r')
