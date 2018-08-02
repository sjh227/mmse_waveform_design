function [thetai_out,wi_out,thetahat, globals] = particle_filter_simple(X,thetai,wi,S,globals)

%get globals
noisesd = globals.noisesd;
L = globals.L;
NR = globals.NR;
Ntargets = globals.Ntargets;
index = globals.index; %the index of this run

npnt = [180; 1770];
Nthr = 0.01;

if size(thetai.phi,1) == 0
    Nparts = npnt(Ntargets,1);
    %load from file
    [theta_out_full] = param_model([],globals);
    thetai_out.alpham = theta_out_full.alpham*ones(180,1);
    thetai_out.alphaa = theta_out_full.alphaa*ones(180,1);
    thetai_out.phi = [-90:89]';
    load('orthovec.mat')
    index2 = orthovec(index);
    strtemp = ['../waveform_design_ntargets/trial15/ortho/wik_'  num2str(index2) 'it30.mat'];
    load(strtemp)
    wi_out = wik(:,1);
    if index == 421
        wi_out = [wi_out(1:76,1); 0;0; wi_out(77:178,1)];
    end
else
    thetai_out = thetai;
    wi_out = wi;
    Nparts = size(wi,1);
    for np=1:Nparts
        thetatemp.phi = thetai_out.phi(np,:);
        thetatemp.alpham = thetai_out.alpham(np,:);
        thetatemp.alphaa = thetai_out.alphaa(np,:);
        mu = H(thetatemp,globals)*S;
        for l=1:L
            wi_out(np,1) = wi_out(np,1) * mvnpdf(real(X(:,l)),real(mu(:,l)),(noisesd^2)*eye(NR)) * mvnpdf(imag(X(:,l)),imag(mu(:,l)),(noisesd^2)*eye(NR));
        end
    end
    wi_out = wi_out/ sum(wi_out); 
    %remove particles with low weight
    %i=1;
    for i=1:180 % i<=size(wi_out,1)
        if wi_out(i,1) < ((1/npnt(Ntargets,1))*Nthr)
            wi_out(i,1) = 0;
            %wi_out = [wi_out(1:i-1,1);wi_out(i+1:size(wi_out,1))];
            %thetai_out.phi = [thetai_out.phi(1:i-1,:);thetai_out.phi(i+1:size(thetai_out.phi,1),:)];
            %thetai_out.alpham = [thetai_out.alpham(1:i-1,:);thetai_out.alpham(i+1:size(thetai_out.alpham,1),:)];
            %thetai_out.alphaa = [thetai_out.alphaa(1:i-1,:);thetai_out.alphaa(i+1:size(thetai_out.alphaa,1),:)];
        %else
        %    i=i+1;
        end
    end
wi_out = wi_out/ sum(wi_out); %normalise weightings
end
wi_out = wi_out/ sum(wi_out); 



thetahat = wi_out'*thetai_out.phi;
%globals.Nparts = Nparts;

    