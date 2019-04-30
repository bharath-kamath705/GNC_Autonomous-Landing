
%aircraft data:altitude, moment of inertia etc
data=xlsread('boeing747_data.xlsx');

%Dimensional derivatives case 1 Mach 0.2 
dd=xlsread('dimensional_derivatives_case1');

%reference conditions [g theta_ref u_ref]
ref=[32.2 0 data(3,1)]; %u_ref is the velocity mentioned in data file

%Aircraft system X_dot=AX+Bu
%states are {del_u w q del_theta]'
[A1,B1]=long_model(data(:,1),dd,ref);

%landing parameters
gsa=3*pi/180; %glide slope angle
rw_d0=30000;
h_ref=1500;
u_ref=data(3,1);

%C matrices and d for reference tracking
C1=[zeros(1,5) 1];
C2=[zeros(1,4) -tan(gsa) 0];
d1=rw_d0*tan(gsa);
d2=h_ref;

% %state equations appended with equations for distance altitude and tracker
% A1_app=[A1 zeros(4,3);
%     1 zeros(1,6);
%     0 1 zeros(1,5);
%     C2-C1 0];
% B1_app=[B1; zeros(3,2)];
% D=[zeros(4,1);u_ref;0;d1-d2];

%adaptive r(t)
del_ws=-20;
C=[tan(gsa) -1 0 0; 0 1 0 0];
A1_app=[A1 zeros(4,2); -C zeros(2,2)];
B1_app =[B1; zeros(2,2)];
D1=[zeros(4,1);-u_ref*tan(gsa); del_ws];

B1_app(2,2)=0.1;
disp('controlability matrix rank')
P = [B1_app A1_app*B1_app A1_app^2*B1_app A1_app^3*B1_app A1_app^4*B1_app A1_app^5*B1_app ];
rank(P)


%eigen structure
%Two desired eigenvalues (dutch roll) are determined by desired damping ratio and
%natural frequency zeta and omega_n respectively
zeta1=0.2; 
omega_n1=3; 
lambda_d(1) = -zeta1*omega_n1 + j*omega_n1*sqrt(1-zeta1^2);
lambda_d(2) = conj(lambda1_d);
lambda_d(3) = -7 + j*0;
lambda_d(4) = 0.1*real(lambda1_d);
lambda_d(5) = 0.31*real(lambda1_d);
lambda_d(6) = 0.41*real(lambda1_d);


% Extracting desired eigenvectors directly from the null-space. 
for i=1:6
mat(:,:,i) = [(lambda_d(i)*eye(6)-A1_app) B1_app];
nullspace(:,:,i) = null(mat(:,:,i),'r');
vu(:,i) = (1+rand(1))*nullspace(:,1,i)+(1+rand(1))*nullspace(:,2,i);
V(:,i)=vu(1:6,i);
U(:,i)=vu(7:8,i);
end

K = U*inv(V); %returns K as complex variable but with 0 imaginary part




x0=[0 0 0 0 0.1 0.1];
[t,x] = ode45('gsa_land',[0:0.02:100],x0,[],A1_app,B1_app,K,D1);%,alpha_T)