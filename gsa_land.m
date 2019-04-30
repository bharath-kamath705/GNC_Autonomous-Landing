function dxdt = gsa_land(t,x,flag,A,B,K,D);


dxdt = (A-B*K)*x+D;

