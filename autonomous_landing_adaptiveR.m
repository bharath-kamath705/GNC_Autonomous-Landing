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
ws=10;
C=[tan(gsa) -1 0 0; 0 1 0 0];
A1_app=[A1 zeros(4,2); -C zeros(2,2)];
B1_app =[B1; zeros(2,2)];
D1=[zeros(4,1);-u_ref*tan(gsa); ws]; %-u_ref*tan(gsa)

B1_app(2,2)=0.2;
disp('controlability matrix rank')
P = [B1_app A1_app*B1_app A1_app^2*B1_app A1_app^3*B1_app A1_app^4*B1_app A1_app^5*B1_app ];
rank(P)


%------------eigen structure----------------------
%desired eigen values
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
%-----------------------------------------------------------------------

% x0=[0 0 0 0 0 0 s_ref h_ref];
% [t,x] = ode45('gsa_land_R',[0:0.02:100],x0,[],A1_app,B1_app,K,D1,u_ref);
% plot(t,x(:,8))

%--------------- Simulation----------------------------------
tol=0.1 %tolerance determines how close to gsa is acceptable
del_ws=0.5 %change in reference input
del_t=1 %time step size
t_tot=185 %total time of simulation
n=t_tot/del_t

ti=0
tf=del_t
X=zeros(1,8) %final states will be stored here
T=0 %final time will be stored here
x0=[0 0 0 0 0 0 s_ref h_ref];

for i=1:n
    
    [t,x] = ode45('gsa_land_R',[ti tf],x0,[],A1_app,B1_app,K,D1,u_ref);
    
    if x(end,8)>(x(end,7)*tan(gsa)+tol)
        D1(end)=D1(end)-del_ws;
    elseif x(end,8)<(x(end,7)*tan(gsa)-tol)
        D1(end)=D1(end)+del_ws;
    else
        ;
    end
    
    X=[X;x];
    T=[T;t];
    ti=tf
    tf=tf+del_t;
    x0=x(end,:);
end


%remember first row of X and T are zeros from initialization 
%first row should be ignored
figure(1)
hold on
plot(X(2:end,7),X(2:end,8)) %plots flight path
plot(X(2:end,7),tan(gsa)*X(2:end,7),'r') %plots desired glide slope
plot(X(2:end,7),zeros(size(X(2:end,7))),'k','LineWidth',2)
plot(0,0,'.','MarkerSize',25)
ylabel('h (ft)','FontSize',16),xlabel('s (ft)','FontSize',16)


%plot of primary states
figure(2),
subplot(4,1,1),plot(T(2:end),u_ref+X(2:end,1))
grid,ylabel('U_{ref}+\Delta u (ft/s)','FontSize',15),xlabel('t(s)')
subplot(4,1,2),plot(T(2:end),X(2:end,2))
grid,ylabel('w (ft/s)','FontSize',15),xlabel('t(s)')
subplot(4,1,3),plot(T(2:end),X(2:end,3))
grid,ylabel('q (rad/s)','FontSize',15),xlabel('t(s)')
subplot(4,1,4),plot(T(2:end),X(2:end,4))
grid,ylabel('\theta (rad)','FontSize',15),xlabel('t(s)')

%plots of control input
%subplot(2,3,5),plot(t2,u),grid,title('\delta_e','FontSize',16),xlabel('t(s)'),ylabel('rad','FontSize',16)
