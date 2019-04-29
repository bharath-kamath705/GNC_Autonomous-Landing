
g=32.2 %acceleration due to gravity

%aircraft data:altitude, moment of inertia etc
data=xlsread('boeing747_data.xlsx')
%case 1
h=data(1,1);
M_num=data(2,1);
V=data(3,1);
W=data(4,1);
I_x=data(5,1); I_y=data(6,1); I_z=data(7,1); I_zx=data(8,1);
xi=data(9,1);
C_d=data(10,1);
m=W;
%reference conditions
theta0=0
u0=V

%Dimensional derivatives case 1 Mach 0.2 
dd=xlsread('dimensional_derivatives_case1')
X_u=dd(1,1); Z_u=dd(1,2); M_u=dd(1,3);
X_w=dd(2,1); Z_w=dd(2,2); M_w=dd(2,3);
X_q=dd(3,1); Z_q=dd(3,2); M_q=dd(3,3);
X_wdot=dd(4,1); Z_wdot=dd(4,2); M_wdot=dd(4,3);
X_e=dd(5,1); Z_e=dd(5,2); M_e=dd(5,3);

%case 1 state matrix A
A1=zeros(4,4)
A1(1,:)=[X_u/m, X_w/m, 0, -g*cos(theta0)]
A1(2,:)=[Z_u, Z_w, Z_q+m*u0,-m*g*sin(theta0)]*inv(m-Z_wdot)
A1(3,1)=inv(I_y)*(M_u+(M_wdot*Z_u*inv(m-Z_wdot)));
A1(3,2)=inv(I_y)*(M_w+(M_wdot*Z_w*inv(m-Z_wdot)));
A1(3,3)=inv(I_y)*(M_q+(M_wdot*(Z_q+m*u0)*inv(m-Z_wdot)));
A1(3,4)=-M_wdot*m*g*sin(theta0)*inv(I_y*(m-Z_wdot));
A1(4,:)=[0, 0, 1, 0]