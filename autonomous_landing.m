
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

%C matrices for reference tracking
C1=[zeros(1,5) 1];
C2=[zeros(1,4) -tan(gsa) 0];
d1=rw_d0*tan(gsa);
d2=h_ref;

%A1 and B1 appended with equations for distance altitude and tracker
A1_app=[A1 zeros(4,3);
    1 zeros(1,6);
    0 1 zeros(1,5);
    C2-C1 0];
B1_app=[B1; zeros(3,2)];
D=[zeros(4,1);u_ref;0;d1-d2];





