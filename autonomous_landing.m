
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

%Dimensional derivatives case 1 Mach 0.2 
dd=xlsread('dimensional_derivatives_case1')
X_u=dd(1,1); Z_u=dd(1,2); M_u=dd(1,3);
X_w=dd(2,1); Z_w=dd(2,2); M_w=dd(2,3);
X_q=dd(3,1); Z_q=dd(3,2); M_q=dd(3,3);
X_wdot=dd(4,1); Z_wdot=dd(4,2); M_wdot=dd(4,3);
X_e=dd(5,1); Z_e=dd(5,2); M_e=dd(5,3);

