 clc
clear all
close all

%aircraft data:altitude, moment of inertia etc
 data=xlsread('boeing747_data.xlsx');

%Dimensional derivatives case 1 Mach 0.2 
 dd=xlsread('dimensional_derivatives_case1');

%reference conditions [g theta_ref u_ref]
ref=[32.2 0 50]; %u_ref data(3,1)is the velocity mentioned in data file

%Aircraft system X_dot=AX+Bu
%states are {del_u w q del_theta]'
[A1,B1]=long_model(data(:,1),dd,ref);

%landing parameters
gsa=3*pi/180; %glide slope angle
s_ref=30000;
h_ref=2000;
u_ref=200 %data(3,1);



%adaptive r(t)
ws=10;
C=[tan(gsa) -1 0 0 0 0; 0 0 0 0 tan(gsa)  -1];
A1_app=[A1 zeros(4,4);
    -1 zeros(1,7);
    0 1 zeros(1,6);
    -C zeros(2,2)];
B1_app=[B1; zeros(4,2)];
D1=[zeros(4,1);-u_ref;0;-u_ref*tan(gsa); 0]; %-u_ref*tan(gsa)

%assume arbitrary values for derivative of Z wrt thrust
B1_app(2,2)=0.2;

%check controllability
disp('controlability matrix rank')
P = [B1_app A1_app*B1_app A1_app^2*B1_app A1_app^3*B1_app A1_app^4*B1_app A1_app^5*B1_app A1_app^6*B1_app A1_app^7*B1_app];
rank(P)


%------------eigen structure----------------------
%desired eigen values
lambda1_d=[-2+j*0.5 -2-j*0.5 -0.05+ j*0 -2.9 -3.5 -4.5 -5 -6]


% Extracting desired eigenvectors directly from the null-space. 
for i=1:length(lambda1_d)
mat(:,:,i) = [(lambda1_d(i)*eye(length(lambda1_d))-A1_app) B1_app];
nullspace(:,:,i) = null(mat(:,:,i),'r');
vu(:,i) = 0.2*i*nullspace(:,1,i)+0.5*nullspace(:,2,i);
V(:,i)=vu(1:length(lambda1_d),i);
U(:,i)=vu(length(lambda1_d)+1:end,i);
end

K1 = U*inv(V); %returns K as complex variable but with 0 imaginary part
K1=real(K1);
%-----------------------------------------------------------------------


%--------------- Simulation----------------------------------
tol=80   %tolerance determines how close to gsa is acceptable
del_ws=0.1 %change in reference input
del_t=1    %time step size
t_tot=220 %total time of simulation
n=t_tot/del_t

ti=0
tf=del_t
X=zeros(1,8) %final states will be stored here
T=0          %final time will be stored here
x0=[0 0 0 0  s_ref h_ref 0 0];

 [t,x] = ode45('gsa_land',[ti tf],x0,[],A1_app,B1_app,K1,D1);

% for i=1:n
%     
%    
%     
% %     if x(end,8)>(x(end,7)*tan(gsa)+tol)
% %         D1(end)=D1(end)-del_ws;
% %     elseif x(end,8)<(x(end,7)*tan(gsa)-tol)
% %         D1(end)=D1(end)+del_ws;
% %     else
% %         ;
% %     end
%     
%     X=[X;x];
%     T=[T;t];
%     ti=tf
%     tf=tf+del_t;
%     x0=x(end,:);
% end


%remember first row of X and T are zeros from initialization 
%first row should be ignored
figure(1)
hold on
%plot flight path,desired glide slope and runway start respectively
plot(x(:,5),x(:,6))              
plot(x(:,5),tan(gsa)*x(:,6),'r') 
plot(x(:,5),zeros(size(x(:,6))),'k','LineWidth',2)
plot(0,0,'.','MarkerSize',25)
ylabel('h (ft)','FontSize',16),xlabel('s (ft)','FontSize',16)
legend({'Flight path','Desired path','Runway','Runway start'},'FontSize',16,'Location','North')


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
figure(3)
u=-X(:,1:6)*K1'
subplot(2,1,1),plot(T(2:end),u(2:end,1))
grid,ylabel('\delta_e (rad)','FontSize',15),xlabel('t(s)')
subplot(2,1,2),plot(T(2:end),u(2:end,2))
grid,ylabel('\delta_p (lb)','FontSize',15),xlabel('t(s)')


