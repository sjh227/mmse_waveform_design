function [Xm, thetam, PXm] = Xtheta_sample(globals,S,theta_in,wi,thetam_in)

NXtheta = globals.NXtheta;
noisesd = globals.noisesd;
Nparts = size(wi,1);
NR = globals.NR;
L = globals.L;
Ntargets = globals.Ntargets;
MCfull = globals.MCfull;

if MCfull
    PXm = 1;
else
    PXm = ones(NXtheta,1);
end

Xm = zeros(NR, L, NXtheta); 
if nargin == 4 
    thetam.phi = zeros(NXtheta, Ntargets);
    thetam.alpham = zeros(NXtheta, Ntargets);
    thetam.alphaa = zeros(NXtheta, Ntargets);
    wicdf = wi;
    for i=2:Nparts
	wicdf(i) = wicdf(i-1) + wicdf(i);
    end
else
    thetam = thetam_in;
end
    

for m=1:NXtheta
    %pick theta
    if nargin == 4
        r = unifrnd(0,1);
        ii = 1;
        while r> wicdf(ii);
            ii = ii + 1;
        end
        thetam.phi(m,:) = theta_in.phi(ii,:);
        thetam.alpham(m,:) = theta_in.alpham(ii,:);
        thetam.alphaa(m,:) = theta_in.alphaa(ii,:);
    end
    thetatemp.phi = thetam.phi(m,:);
    thetatemp.alpham = thetam.alpham(m,:);
    thetatemp.alphaa = thetam.alphaa(m,:);
    for l=1:L
        mu = H(thetatemp,globals)*S(:,l);
        Xkr = zeros(NR,1);
        Xki = zeros(NR,1);
        for nr=1:NR
            Xkr(nr,1) = normrnd(real(mu(nr,1)),noisesd);
            Xki(nr,1) = normrnd(imag(mu(nr,1)),noisesd);
        end
        Xm(:,l,m) = Xkr+sqrt(-1)*Xki;
        if ~MCfull
            PXm(m,1) = PXm(m,1)*mvnpdf([real(Xm(:,l,m));imag(Xm(:,l,m))],[real(mu(:,1));imag(mu(:,1))],(noisesd^2)*eye(2*NR));
        end
    end


end
