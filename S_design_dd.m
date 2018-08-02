function [S_out] = S_design_dd(globals, thetai,wi,Si)
%adaptively design waveform using directional derivative to handle
%constraint

P = globals.P;
L = globals.L;
NT = globals.NT;
tolerance = globals.tolerance;
stepsize = globals.stepsize0;
kglob = globals.kglob;
MCfull = globals.MCfull;
maxits = globals.maxits;
printops = globals.printops;
%ds = zeros(2*NT*L,1);

if nargin == 1
    if NT == L
        T = randn(NT); 
        [S_out , ~, ~] = svd(T);
        S_out = (sqrt(P*L)/sqrt(trace(S_out*S_out')))*S_out;
    else
        S_out = eye(min(NT, L));
        S_out = padarray(S_out,[NT - min(NT,L), L - min(NT,L)],'post');
        S_out = (sqrt(P*L)/sqrt(trace(S_out*S_out')))*S_out;
    end
else
    if nargin == 3
        S = eye(min(NT, L));
        S = padarray(S,[NT - min(NT,L), L - min(NT,L)],'post');
        S = (sqrt(P*L)/sqrt(trace(S*S')))*S;      
    else
        S = Si;
    end
    [Xm, thetam, PXm] = Xtheta_sample(globals,(sqrt(P*L)/sqrt(trace(S*S')))*S,thetai,wi);
    [cost] =  waveform_cost((sqrt(P*L)/sqrt(trace(S*S')))*S, thetai, wi, Xm, thetam, PXm, globals);
    if printops
        disp(['cost_',num2str(kglob),'_t00: ',num2str(cost)])
    end
    needdir = true;
    its = 0;
    while stepsize > tolerance && its < maxits
        if needdir %only recalculate dir if new location
            [dir] = dcost_dd(thetai,wi,S, Xm, thetam,PXm,globals); %dir is normalised
        end
        S_temp = S - stepsize*dir;
        S_temp = (sqrt(P*L)/sqrt(trace(S_temp*S_temp')))*S_temp; %normalise to return to surface of hyper-sphere
        [cost_temp] = waveform_cost((sqrt(P*L)/sqrt(trace(S_temp*S_temp')))*S_temp, thetai, wi,Xm,thetam,PXm, globals);
        if printops
            if its < 9
                disp(['cost_', num2str(kglob),'_t0', num2str(its+1) ,': ' , num2str(cost_temp)])
            else
                disp(['cost_', num2str(kglob),'_t', num2str(its+1) ,': ' , num2str(cost_temp)])
            end
        end
        if cost_temp < cost
            S = S_temp; 
            cost = cost_temp;
            needdir = true;
        else
            stepsize = (1/2)*stepsize;
            needdir = false; 
        end
        its = its + 1;
    end
    S_out = S;
end

