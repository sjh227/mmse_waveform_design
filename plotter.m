function [voidvar] = plotter(wik,thetaik,Sk,binwidth,ind)
%currently only for 1 target
close all

Ncols = 3; %set here

%other terms, no need to set
Ktrue = size(wik,2);
K = Ncols*floor(Ktrue/Ncols);
Nparts = size(wik,1);
NT = size(Sk,1);
L = size(Sk,2)/Ktrue;
Ntargets = size(thetaik,2)/Ktrue;
Nbbw = [180;0;60;0;0;0;0;0;20];
binbase1 = [-90:1:89]';
binbase3 = [-89:3:88]';
binbase9 = [-86:9:85]';

%binints = globals.binints;
%binplot = globals.phivec;
%Nbins = size(binplot,1);

voidvar = 1;

figure
maxwave = -inf;
minwave = inf;
maxpdf = -inf;

for k=1:K
    %waveform plot
    wavevec = zeros(181,1);
    for thetaplot=-90:1:90
        at = zeros(NT,1);
        for nt=1:NT
            at(nt,1) = exp(-(nt-1)*pi*sqrt(-1)*cosd(90+thetaplot));
        end
        wavevec(thetaplot+91,1) = 10*log10(abs((1/L)*at'*conj(Sk(:,L*(k-1)+1:L*k)*Sk(:,L*(k-1)+1:L*k)')*at));
    end
    maxwave = max([maxwave;wavevec]);
    minwave = min([minwave;wavevec]);
    wavepos = 2*floor((k-1)/Ncols)*Ncols + 1 + mod(k-1,Ncols);
    subplot(2*(K/Ncols),Ncols,wavepos)
    plot([-90:1:90]',wavevec,'k')
    
    %pdf plot
    colourphi1 = ['k','r','b','g']';
    thetaikq = thetaik(:,(k-1)*Ntargets+1:k*Ntargets);
    pdfpos = 2*floor((k-1)/Ncols)*Ncols + Ncols + 1 + mod(k-1,Ncols);
    subplot(2*(K/Ncols),Ncols,pdfpos)
    hold on
    for q=1:Ntargets
        Nbins = Nbbw(binwidth(k,1),1);
        pdfvec = zeros(Nbins,1);
        for i=1:Nparts
            if binwidth(k,1) == 1
                j = thetaikq(i,q) + 91;
                pdfvec(j,1) = pdfvec(j,1) + wik(i,k)/binwidth(k,1);
            elseif binwidth(k,1) == 3
                if thetaikq(i,q) ~= 0
                    j = (thetaikq(i,q) + 92)/3;
                    pdfvec(j,1) = pdfvec(j,1) + wik(i,k)/binwidth(k,1);
                end
            elseif binwidth(k,1) == 9
                if thetaikq(i,q) ~= 0
                    j = (thetaikq(i,q) + 95)/9;
                    pdfvec(j,1) = pdfvec(j,1) + wik(i,k)/binwidth(k,1);
                end
            end
        end
        maxpdf = max([maxpdf;pdfvec]);
        if binwidth(k,1) == 1
            plot(binbase1,pdfvec,colourphi1(q,1))
        elseif binwidth(k,1) == 3
            if q == 1
                %plot(binbase3,pdfvec,colourphi1(q,1))
                plot(binbase3,pdfvec,'b')
            elseif q==2
                %plot(binbase3,pdfvec,'Color',[0.1,0.1,0.1])
                plot(binbase3,pdfvec,'m')
            end
        elseif binwidth(k,1) == 9
            plot(binbase9,pdfvec,colourphi1(q,1))
        end
    end
end

ymaxwave = maxwave + abs(0.1*maxwave);
yminwave = minwave - abs(0.1*maxwave);

for k=1:K
    wavepos = 2*floor((k-1)/Ncols)*Ncols + 1 + mod(k-1,Ncols);
    pdfpos = 2*floor((k-1)/Ncols)*Ncols + Ncols + 1 + mod(k-1,Ncols);
    subplot(2*(K/Ncols),Ncols,wavepos)
    hold on
    axis([-90 90 yminwave ymaxwave])
    targetglobals.Ntargets = Ntargets;
    targettheta = param_model([],targetglobals);
    for i=1:Ntargets
        if i==1
            plot([targettheta.phi(1,i), targettheta.phi(1,i)],[yminwave, ymaxwave],'b:')
        elseif i==2
            %plot([targettheta.phi(1,i), targettheta.phi(1,i)],[yminwave, ymaxwave],'Color',[0.1,0.1,0.1],'LineStyle','--')
            plot([targettheta.phi(1,i), targettheta.phi(1,i)],[yminwave, ymaxwave],'m:')
        elseif i==3
            plot([targettheta.phi(1,i), targettheta.phi(1,i)],[yminwave, ymaxwave],'b:')
        end
    end
    if nargin == 5
        title(['Iteration ',  num2str(ind)],'FontWeight','Normal')
    else
        title(['Pulse ',  num2str(k)],'FontWeight','Normal')
    end
    ylabel('Power/ dB')
    %xlabel('\phi/ ^O')
    subplot(2*(K/Ncols),Ncols,pdfpos)
    hold on
    axis([-90 90 0 1.1*maxpdf])
    targetglobals.Ntargets = Ntargets;
    targettheta = param_model([],targetglobals);
    for i=1:Ntargets
        if i==1
            plot([targettheta.phi(1,i), targettheta.phi(1,i)],[0, 1.1*maxpdf],'b:')
        elseif i==2
            plot([targettheta.phi(1,i), targettheta.phi(1,i)],[0, 1.1*maxpdf],'m:')
        elseif i==3
            plot([targettheta.phi(1,i), targettheta.phi(1,i)],[0, 1.1*maxpdf],'b:')
        end
    end
    xlabel('\phi/ ^o')
    ylabel('p(\phi)')
    box on
end

end
