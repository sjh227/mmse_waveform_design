function [cost] = waveform_cost(Si, thetai, wi, Xm, thetam, PXm, globals)

NXtheta = globals.NXtheta;
L = globals.L;
noisesd = globals.noisesd;
NR = globals.NR;
MCfull = globals.MCfull;
PXloc = ones(NXtheta,1);

cost = 0;

for m=1:NXtheta
    [nom1, denom1] = terms_eval(Xm,m,Si,wi,thetai,globals);
    thetatemp.phi = thetam.phi(m,:);
    thetatemp.alpham = thetam.alpham(m,:);
    thetatemp.alphaa = thetam.alphaa(m,:);
    if ~MCfull
        for l = 1:L
            mu = H(thetatemp,globals)*Si(:,l);
            PXloc(m,1) = PXloc(m,1)*mvnpdf([real(Xm(:,l,m));imag(Xm(:,l,m))],[real(mu(:,1));imag(mu(:,1))],(noisesd^2)*eye(2*NR));
        end
    end
    cost = cost + (PXloc(m,1)/PXm(m,1))*((nom1/denom1) - [thetam.phi(m,:),thetam.alpham(m,:),thetam.alphaa(m,:)])*transpose((nom1/denom1) - [thetam.phi(m,:),thetam.alpham(m,:),thetam.alphaa(m,:)]); %row vector
end


cost = cost/sum(PXloc./PXm);



