
%aircraft data:altitude, moment of inertia etc
data=xlsread('boeing747_data.xlsx')

%Dimensional derivatives case 1 Mach 0.2 
dd=xlsread('dimensional_derivatives_case1')

%reference conditions [g theta_ref u_ref]
ref=[32.2 0 data(3,1)] %u_ref is the velocity mentioned in data file



A1=long_model(data(:,1),dd,ref)