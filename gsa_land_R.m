function dxdt = gsa_land_R(t,x,flag,A,B,K,D,u_ref);




dxdt=zeros(8,1);
dxdt(1:6) = (A-B*K)*x(1:6)+D;
dxdt(7)=-u_ref-x(1);
dxdt(8)=-x(2);
