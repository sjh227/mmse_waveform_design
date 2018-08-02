function [S_out] = S_design_rmb(globals, thetai,wi,Si)
%adaptively design waveform using directional derivative to handle
%constraint

P = globals.P;
L = globals.L;
NT = globals.NT;
%tolerance =  (1/L)*(globals.tolerance)^2; %squared because RMB operates on Rsk
stepsize = (globals.stepsize0)^2;
tolerance = stepsize*(1/2)^9;
kglob = globals.kglob;
maxits = globals.maxits;
printops = globals.printops;
%ds = zeros(2*NT*L,1);

if nargin == 1
    %if NT == L
    %    T = randn(NT); 
    %    [S_out , ~, ~] = svd(T);
    %    S_out = (sqrt(P*L)/sqrt(trace(S_out*S_out')))*S_out;
    %else
        S_out = eye(min(NT, L));
        S_out = padarray(S_out,[NT - min(NT,L), L - min(NT,L)],'post');
        S_out = (sqrt(P*L)/sqrt(trace(S_out*S_out')))*S_out;
    %end
else
    if nargin == 3
        S = eye(min(NT, L));
        S = padarray(S,[NT - min(NT,L), L - min(NT,L)],'post');
        S = (sqrt(P*L)/sqrt(trace(S*S')))*S;      
    else
        S = Si;
    end
    %[Xm, thetam, PXm] = Xtheta_sample(globals,(sqrt(P*L)/sqrt(trace(S*S')))*S,thetai,wi);
    Rs = (1/L)*S*S';
    globals.Rs = Rs; %save as global to pass to optimiser
    [T, Dinv] = DandT(thetai,wi,globals);
    cost = trace(T*Dinv*transpose(T));
    %[cost] =  waveform_cost((sqrt(P*L)/sqrt(trace(S*S')))*S, thetai, wi, Xm, thetam, PXm, globals);
    if printops
        disp(['cost_',num2str(kglob),'_t00: ',num2str(cost)])
    end
    needdir = true;
    its = 0;
    while stepsize > tolerance && its < maxits
        Rs = (1/L)*(S*S');
        globals.Rs = Rs; %save as global to pass to optimiser
        if needdir %only recalculate dir if new location
            [dir,~] = dcost_rmb(thetai,wi,T,Dinv,globals); %dir is normalised
%             if printops
%                 dirop = dir
%             end
        end
        Rs_temp = Rs + stepsize*dir;
        [~,pd] = chol(L*Rs_temp);
        if pd ~= 0 
            stepsize = (1/2)*stepsize;
            needdir = false; 
            its = its + 1;
            if printops
                if its < 9
                    disp(['inval', num2str(kglob),'_t0', num2str(its) ])
                else
                    disp(['inval', num2str(kglob),'_t', num2str(its) ])
                end
            end
        else
            S_temp = chol(L*Rs_temp);
            S_temp = (sqrt(P*L)/sqrt(trace(S_temp*S_temp')))*S_temp; %renormalise to be on the safe side
            cost_temp = trace(T*Dinv*transpose(T));
            %[cost_temp] = waveform_cost(S_temp, thetai, wi,Xm,thetam,PXm, globals);
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
    end
    S_out = S;
end

