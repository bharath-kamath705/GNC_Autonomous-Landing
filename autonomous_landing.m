
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
del_u_ref=-10
C=[1 0 0 0; 0 1 0 0];
A1_app=[A1 zeros(4,2); -C zeros(2,2)]
B1_app =[B1; zeros(2,2)]
D1=[zeros(4,1);del_u_ref; del_u_ref*tan(gsa)]


%eigen structure
%Two desired eigenvalues (dutch roll) are determined by desired damping ratio and
%natural frequency zeta and omega_n respectively
zeta1=0.2; 
omega_n1=2; 
lambda1_d = -zeta1*omega_n1 + j*omega_n1*sqrt(1-zeta1^2);
lambda2_d = conj(lambda1_d);
%roll eigen value is desired to be at -7
lambda3_d = -7 + j*0;
%spiral mode arbitrarily chosen 
lambda4_d = 0.1*real(lambda1_d);
lambda5_d = 0.1*real(lambda1_d);
lambda6_d = 0.1*real(lambda1_d);


% Extracting desired eigenvectors directly from the null-space. 
mat1 = [(lambda1_d*eye(6)-A1_app) B1_app];
nullspace1 = null(mat1,'r');
v1u1 = 0.22*nullspace1(:,1)+0.3*nullspace1(:,2);

mat2 = [(lambda2_d*eye(6)-A1_app) B1_app];
nullspace2 = null(mat2,'r');
v2u2 = 0.62*nullspace2(:,1)+0.3*nullspace2(:,2);

mat3 = [(lambda3_d*eye(6)-A1_app) B1_app];
nullspace3 = null(mat3,'r');
v3u3 = 0.12*nullspace3(:,1)+0.3*nullspace3(:,2);

mat4 = [(lambda4_d*eye(6)-A1_app) B1_app];
nullspace4 = null(mat4,'r');
v4u4 = 0.52*nullspace4(:,1)+0.3*nullspace4(:,2);

mat5 = [(lambda5_d*eye(6)-A1_app) B1_app];
nullspace5 = null(mat5,'r');
v5u5 = 0.77*nullspace5(:,1)+0.3*nullspace5(:,2);

mat6 = [(lambda6_d*eye(6)-A1_app) B1_app];
nullspace6 = null(mat6,'r');
v6u6 = 0.9*nullspace6(:,1)+0.3*nullspace6(:,2);

v1 = v1u1(1:6); u1 = v1u1(7:8);
v2 = v2u2(1:6); u2 = v2u2(7:8);
v3 = v3u3(1:6); u3 = v3u3(7:8);
v4 = v4u4(1:6); u4 = v4u4(7:8);
v5 = v5u5(1:6); u5 = v5u5(7:8);
v6 = v6u6(1:6); u6 = v6u6(7:8);

V = [v1 v2 v3 v4 v5 v6 ];
U = [u1 u2 u3 u4 u5 u6 ];

K = U*inv(V); %returns K as complex variable but with 0 imaginary part
K=real(K)

x0=[0 0 0 0 0.1 0.1];
[t,x] = ode45('gsa_land',[0:0.02:100],x0,[],A1_app,B1_app,K,D1);%,alpha_T)