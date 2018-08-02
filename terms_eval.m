function [nom, denom] = terms_eval(Xm,m,S,wi,thetai,globals)

NR = globals.NR;
L = globals.L;
Nparts = size(wi,1);
noisesd = globals.noisesd;
Ntargets = globals.Ntargets;

nom = zeros(1,3*Ntargets);
denom = 0;

for np=1:Nparts
    PX = 1;
    thetatemp.phi = thetai.phi(np,:);
    thetatemp.alpham = thetai.alpham(np,:);
    thetatemp.alphaa = thetai.alphaa(np,:);
    for l = 1:L
        mu = H(thetatemp,globals)*S(:,l);
        PX = PX*mvnpdf([real(Xm(:,l,m));imag(Xm(:,l,m))],[real(mu(:,1));imag(mu(:,1))],(noisesd^2)*eye(2*NR));
    end
    nom = nom + wi(np)*[thetai.phi(np,:),thetai.alpham(np,:),thetai.alphaa(np,:)]*PX;
    denom = denom + wi(np)*PX;
end
    