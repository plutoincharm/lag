function qdot = fun4barLagcODE(q,t)

global m L r MI T2a tf

m2 = m(1);
m3 = m(2);
m4 = m(3);
MI2 = MI(1);
MI3 = MI(2);
MI4 = MI(3);
L2 = L(1);
L3 = L(2);
L4 = L(3);
r2 = r(1);
r3 = r(2);
r4 = r(3);

T2 = T2a*t*(tf-t);
% T2 = T2a*cos(t);

omg2 = q(1);
omg3 = q(2);
omg4 = q(3);
tht2 = q(4);
tht3 = q(5);
tht4 = q(6);

Mmat=[L2.^2.*m3+MI2+m2.*r2.^2,L2.*m3.*r3.*cos(tht2+(-1).*tht3),0;L2.*...
      m3.*r3.*cos(tht2+(-1).*tht3),MI3+m3.*r3.^2,0;0,0,MI4+m4.*r4.^2];

Cq=[(-1).*L2.*sin(tht2),(-1).*L3.*sin(tht3),L4.*sin(tht4);L2.*cos( ...
     tht2),L3.*cos(tht3),(-1).*L4.*cos(tht4)];
 
rhs1=[T2+(-1).*L2.*m3.*omg3.^2.*r3.*sin(tht2+(-1).*tht3);L2.*m3.* ...
      omg2.^2.*r3.*sin(tht2+(-1).*tht3);0];  
  
rhs2=[L2.*omg2.^2.*cos(tht2)+L3.*omg3.^2.*cos(tht3)+(-1).*L4.*omg4.^2.* ...
      cos(tht4);L2.*omg2.^2.*sin(tht2)+L3.*omg3.^2.*sin(tht3)+(-1).*L4.* ...
      omg4.^2.*sin(tht4)];

invM = inv(Mmat); 
Gmat = Cq*invM*Cq';

% if rcond(Gmat)< 1e-4
%     rcond(Gmat)
%     det(Gmat)
%     t
%     pause
% end

lam = Gmat\(rhs2 - Cq*invM*rhs1);

alp = invM*(Cq'*lam+rhs1);

qdot = [alp;omg2;omg3;omg4];
  
  
  
  
  
  
  
  
  