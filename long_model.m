function [A,B] = long_model(data,dd,ref)

%reference conditions
g=ref(1);
theta_ref=ref(2);
u_ref=ref(3);

%aircraft data
h=data(1);
M_num=data(2);
V=data(3);
W=data(4);
I_x=data(5); I_y=data(6); I_z=data(7); I_zx=data(8);
xi=data(9);
C_d=data(10);
m=W;

%dimensional derivatives
X_u=dd(1,1); Z_u=dd(1,2); M_u=dd(1,3);
X_w=dd(2,1); Z_w=dd(2,2); M_w=dd(2,3);
X_q=dd(3,1); Z_q=dd(3,2); M_q=dd(3,3);
X_wdot=dd(4,1); Z_wdot=dd(4,2); M_wdot=dd(4,3);
X_e=dd(5,1); Z_e=dd(5,2); M_e=dd(5,3);

X_p=-0.3*m*g; %(derivative w.r.t thrust input)


%compute state matrix A
A=zeros(4,4);
A(1,:)=[X_u/m, X_w/m, 0, -g*cos(theta_ref)];
A(2,:)=[Z_u, Z_w, Z_q+m*u_ref,-m*g*sin(theta_ref)]*inv(m-Z_wdot);
A(3,1)=inv(I_y)*(M_u+(M_wdot*Z_u*inv(m-Z_wdot)));
A(3,2)=inv(I_y)*(M_w+(M_wdot*Z_w*inv(m-Z_wdot)));
A(3,3)=inv(I_y)*(M_q+(M_wdot*(Z_q+m*u_ref)*inv(m-Z_wdot)));
A(3,4)=-M_wdot*m*g*sin(theta_ref)*inv(I_y*(m-Z_wdot));
A(4,:)=[0, 0, 1, 0];

%compute matrix B
B=zeros(4,2);
B(1,:)=inv(m)*[X_e, X_p];
B(2,:)=inv(m-Z_wdot)*[Z_e, 0];
B(3,:)=inv(I_y)*[M_e+M_wdot*Z_e*inv(m-Z_wdot), 0];
B(4,:)=[0,0];
end