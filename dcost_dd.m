function [dir_out] = dcost_dd(thetai,wi,S, Xm, thetam, PXm, globals)
%outputs a unit vector in direction of steepest descent

NXtheta = globals.NXtheta;
NT = globals.NT;
L = globals.L;
P = globals.P;
noisesd = globals.noisesd;
NR = globals.NR;

dir_out = zeros(NT,L);

sr = []; si = [];
for l=1:L
    sr = [sr;real(S(:,l))];
    si = [si;imag(S(:,l))];
end
s = [sr;si]; %express S as vector s'


PX = ones(NXtheta,1); %for later use
for m=1:NXtheta
    thetatemp.phi = thetam.phi(m,:);
    thetatemp.alpham = thetam.alpham(m,:);
    thetatemp.alphaa = thetam.alphaa(m,:);
    for l = 1:L
        mu = H(thetatemp,globals)*S(:,l);
        PX(m,1) = PX(m,1)*mvnpdf([real(Xm(:,l,m));imag(Xm(:,l,m))],[real(mu(:,1));imag(mu(:,1))],(noisesd^2)*eye(2*NR));
    end
end

PXnorm = sum(PX./PXm);


v2 = zeros(NR*L*2,1);
for m=1:NXtheta
    thetatemp.phi = thetam.phi(m,:);
    thetatemp.alpham = thetam.alpham(m,:);
    thetatemp.alphaa = thetam.alphaa(m,:);
    v2 = v2 + dy(Xm,m,thetatemp,S,globals)*(PX(m,1)/PXm(m,1))/PXnorm;
end

ds = zeros(NT*2*L,1);

for m=1:NXtheta
    thetatemp.phi = thetam.phi(m,:);
    thetatemp.alpham = thetam.alpham(m,:);
    thetatemp.alphaa = thetam.alphaa(m,:);
    
    [nom1, denom1] = terms_eval(Xm,m,S,wi,thetai,globals);
    u = ((nom1/denom1) - [thetam.phi(m,:),thetam.alpham(m,:),thetam.alphaa(m,:)])*transpose((nom1/denom1) - [thetam.phi(m,:),thetam.alpham(m,:),thetam.alphaa(m,:)]); 
    du = dUvec(Xm,nom1,denom1,thetam,m,S,wi,thetai,globals);
    v = (PX(m,1)/(PXm(m,1)))/PXnorm;
    dv = dy(Xm,m,thetatemp,S,globals)*v - v*v2;

    ds = ds + u*dv + v*du; %(duv + U*dvv)*(PX(m,1)/(PXm(m,1)*PXnorm));

end


%find component of ds tangental to s
ds = ds - ((ds'*s)/(s'*s))*s;

%normalise
ds = ds/(sqrt(ds'*ds));


%turn back into matrix form
for n=1:NT*L
    nt = 1+mod(n-1,NT);
    l = 1+floor((n-1)/NT);
    dir_out(nt,l) = ds(n,1) + ds(n+NT*L,1)*sqrt(-1);
end




