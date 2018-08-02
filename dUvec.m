function [dUout] = dUvec(Xm,nom1,denom1,thetam,m,S,wi,thetai,globals)

%global Nparts noisesd L NR Ntargets alphrand

Nparts = size(wi,1);
noisesd = globals.noisesd;
L = globals.L;
NR = globals.NR;
Ntargets = globals.Ntargets;

%dUout = zeros(NR*L*2,1);

%nom2 = zeros(NR*L*2,1);
nom3 = zeros(NR*L*2,1);
nom2 = zeros(NR*L*2,3*Ntargets);


for i=1:Nparts
    PX = 1;
    thetatemp.phi = thetai.phi(i,:);
    thetatemp.alpham = thetai.alpham(i,:);
    thetatemp.alphaa = thetai.alphaa(i,:);
    for l = 1:L
        mu = H(thetatemp,globals)*S(:,l);
        PX = PX*mvnpdf([real(Xm(:,l,m));imag(Xm(:,l,m))],[real(mu(:,1));imag(mu(:,1))],(noisesd^2)*eye(NR*2));
    end   
    dVloc = dy(Xm,m,thetatemp,S,globals);
    nom2 = nom2 + wi(i)*PX*dVloc*[thetai.phi(i,:),thetai.alpham(i,:),thetai.alphaa(i,:)];
    nom3 = nom3 + wi(i)*PX*dVloc;
end
dUout = 2*((nom2/denom1) - (nom3*nom1/(denom1^2)))*transpose(((nom1/denom1)-[thetam.phi(m,:),thetam.alpham(m,:),thetam.alphaa(m,:)]));
