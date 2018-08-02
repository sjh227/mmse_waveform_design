function [dVout] = dy(Xm,m,theta,S,globals)


noisesd = globals.noisesd;
NR = globals.NR;
L = globals.L;

Rninv = ((noisesd)^-2)*eye(2*NR*L);

x = [];
for l=1:L
    x = [x;real(Xm(:,l,m))];
end
for l=1:L
    x = [x;imag(Xm(:,l,m))];
end
s = [];
for l=1:L
    s = [s;real(S(:,l))];
end
for l=1:L
    s = [s;imag(S(:,l))];
end

H0 = H(theta,globals);
H2 = [];
for l=1:L
    H2 = blkdiag(H2,H0);
end
H1 = [real(H2),-imag(H2);imag(H2),real(H2)];


dVout = 2*transpose(H1)*Rninv*x - 2*transpose(H1)*Rninv*H1*s;



