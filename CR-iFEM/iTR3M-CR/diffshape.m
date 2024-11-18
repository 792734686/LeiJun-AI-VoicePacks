syms b1 b2 b3 c1 c2 c3 s t
L1 = 1/3-s-t;
L2 = s+1/3;
L3 = t+1/3;
Nu1 = 0.5*L1*(b3*L2-b2*L3);
Nu2 = 0.5*L2*(b1*L3-b3*L1);
Nu3 = 0.5*L3*(b2*L1-b1*L2);
Nv1 = 0.5*L1*(c3*L2-c2*L3);
Nv2 = 0.5*L2*(c1*L3-c3*L1);
Nv3 = 0.5*L3*(c2*L1-c1*L2);
dNs = diff([Nu1,Nu2,Nu3;Nv1,Nv2,Nv3],s);
dNt = diff([Nu1,Nu2,Nu3;Nv1,Nv2,Nv3],t);
