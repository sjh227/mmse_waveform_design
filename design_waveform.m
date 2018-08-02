function [root_squared_error] = design_waveform(design_type)
%function to adaptively design waveform
close all
design_type = str2num(design_type); %condor passes in a string
%tout = []; %initialise incase not set later

%initialise all implicit 'inputs' as local variables (global variables don't work in compiled matlab, ie for Condor
ASNR = -3;
folder_index = 15;
globals.noisesd = 1/sqrt(2); %fixed: this is for each of the real and imaginary components so noise var is one
globals.NXtheta = 250; %Number of samples of the pair (X, theta) for Monte-Carlo integration
globals.L = 1; %ie number of snapshots to be designed before next waveform update
globals.NT = 5; %number of transmit elements
globals.NR = 5; %number of receive elements
%if design_type == 5
%   globals.L = max(globals.NT,globals.NR); %needed to determine Sk from Rsk in RMB method
%end
globals.P = (10^(ASNR/10))/(globals.NR*globals.L); %note that |\alpha|^2 is defined as equal to 1
if design_type == 5
    globals.L = max(globals.NT,globals.NR); %needed to determine Sk from Rsk in RMB method
    globals.P = (10^(ASNR/10))/(globals.NR*globals.L); %note that |\alpha|^2 is defined as equal to 1
end
globals.stepsize0 = (1/2)*sqrt(globals.P*globals.L);
globals.tolerance = globals.stepsize0*(1/2)^9;%((1/2)^12)*sqrt(globals.P/globals.NT);
globals.K = 30; %number of iterations of waveform design process
globals.MCfull = false; %if true, resample each point of descent, if false sample once
globals.maxits = 40; %in no line search maximum iterations
globals.printops = false; %whether to print outputs inc figs - set false for condor
globals.Ntargets = 1;
globals.linesearch = false; %for true use line search, for false simple adaptive stepsize
globals.J = 10; %for RMB only
%*****PF parameters set in PF******

%initialise other variables
theta = []; %initialise theta as empty
thetai.phi = []; thetai.alphaa = []; thetai.alpham = []; wi = []; %initialise particle filter particles and weights
globals.kglob = 0; %make k part of globals for displaying results in other functions
[S] = S_design_dd(globals); %initialise S - same for all design types
S0 = S; %store initialised S, so that gradient descent can occur from S_0 rather than S_{k-1} if desired

if ~globals.printops
    i=1; search_ind = true;
    while search_ind
    strtemp = ['./trial15/rmb_ind/ind_'  num2str(i) '.mat'];
        if ~(exist(fullfile(cd,strtemp),'file') ==2)
            num_save = i;
            globals.index = i;
            save_dummy = [1;0];
            save(['./trial' num2str(folder_index) '/rmb_ind/ind_' num2str(i) '.mat'],'save_dummy')
            search_ind = false;
            rng(i);
        else
            i=i+1;
        end
    end
else
    globals.index = 1;
end

% num_save = 421;
% globals.index = 421;
% save_dummy = [1;0];
% save(['./trial' num2str(folder_index) '/rmb_ind/ind_' num2str(421) '.mat'],'save_dummy')
% rng(421);

%num_save= randi(5000);


%vectors to hold all values of estimates and waveform for later plotting
%wik = [];
%thetaik = [];
%Sk = [];
squared_error = [];
%thetahatv = [];
%binwidth = [];

%check for folders and create if not existing
%if ~globals.printops
%    if design_type == 0
%        if exist(['ortho' num2str(folder_index)]) ~= 7
%            mkdir(['ortho' num2str(folder_index)]);
%        end
%    else
%        if exist(['adapt' num2str(folder_index)]) ~= 7
%            mkdir(['adapt' num2str(folder_index)]);
%        end
%    end
%end
    
K = globals.K;
for k=1:K
    %tic
    globals.kglob = k;
    theta = param_model(theta,globals); %evolve or initialise theta
    [X] = meas_model(theta,S,globals); %generate random observation X
    [thetai,wi,thetahat,globals] = particle_filter_simple(X,thetai,wi,S,globals); %estimate the parameter using a particle filter
    %binwidth = [binwidth;globals.binwidth];
%     if k == 1 && design_type == 5
%     %    sizewik = size(wi,1);
%         globals.L = max(globals.NT,globals.NR); %needed to determine Sk from Rsk in RMB method
%         globals.P = (10^(ASNR/10))/(globals.NR*globals.L); %note that |\alpha|^2 is defined as equal to 1
%         globals.stepsize0 = 2*pi*(1/8)*sqrt(globals.P*globals.L);
%         globals.tolerance = globals.stepsize0*(1/2)^9;
%     end
    %disp(['np=', num2str(size(wi,1))])
    %wik = [wik,[wi;zeros(sizewik-size(wi,1),1)]];
    %thetaik = [thetaik, [thetai.phi;zeros(sizewik-size(wi,1),globals.Ntargets)]];
    %Sk = [Sk,S];
    squared_error = [squared_error; (thetahat - theta.phi)*(thetahat - theta.phi)']; %currently only works for one target
    if ~globals.printops && design_type == 5
            save(['./trial' num2str(folder_index) '/rmb/se_' num2str(num_save) 'it' num2str(k) '.mat'],'squared_error')
    end
    %thetahatv = [thetahatv; thetahat];
    if k<globals.K %no need to design the final waveform
        if design_type == 1
            [S] = S_design(globals,thetai,wi); %design the next waveform our method descending from identity
        elseif design_type == 2
            [S] = S_design(globals,thetai,wi,S); %design the next waveform - our method descending from previous S
        elseif design_type == 3
            [S] = S_design_mmse(globals,thetahat); %mmse design - only for 1 target
        elseif design_type == 4
            [S] = S_design_inversion(globals,X); %inversion waveform approach
        elseif design_type == 5
            [S] = S_design_rmb(globals,thetai,wi); %bounded mmse design method
        elseif design_type == 6
            [S] = S_design_dd(globals,thetai,wi); %gradient descent using directional derivative over hyper-sphere of max power
        elseif design_type == 7
            [S] = S_design_rmb_rand(globals,thetai,wi); %random method based on rmb format
        else
            [S] = S_design(globals); %not adaptive design type = 0
        end
    end
    if ~globals.printops && k == K
        if design_type == 0
            save(['./trial' num2str(folder_index) '/ortho/se_' num2str(num_save) 'it' num2str(k) '.mat'],'squared_error')
            save(['./trial' num2str(folder_index) '/ortho/Sk_' num2str(num_save) 'it' num2str(k) '.mat'],'Sk')
            save(['./trial' num2str(folder_index) '/ortho/wik_' num2str(num_save) 'it' num2str(k) '.mat'],'wik')
            save(['./trial' num2str(folder_index) '/ortho/thetaik_' num2str(num_save) 'it' num2str(k) '.mat'],'thetaik')
            save(['./trial' num2str(folder_index) '/ortho/thetahatv_' num2str(num_save) 'it' num2str(k) '.mat'],'thetahatv')
            save(['./trial' num2str(folder_index) '/ortho/binwidth_' num2str(num_save) 'it' num2str(k) '.mat'],'binwidth')
        elseif design_type == 1
            save(['./trial' num2str(folder_index) '/adapt/se_' num2str(num_save) 'it' num2str(k) '.mat'],'squared_error')
            save(['./trial' num2str(folder_index) '/adapt/Sk_' num2str(num_save) 'it' num2str(k) '.mat'],'Sk')
            save(['./trial' num2str(folder_index) '/adapt/wik_' num2str(num_save) 'it' num2str(k) '.mat'],'wik')
            save(['./trial' num2str(folder_index) '/adapt/thetaik_' num2str(num_save) 'it' num2str(k) '.mat'],'thetaik')
            save(['./trial' num2str(folder_index) '/adapt/thetahatv_' num2str(num_save) 'it' num2str(k) '.mat'],'thetahatv')
            save(['./trial' num2str(folder_index) '/adapt/binwidth_' num2str(num_save) 'it' num2str(k) '.mat'],'binwidth')
        elseif design_type == 3
            save(['./trial' num2str(folder_index) '/mmse/se_' num2str(num_save) 'it' num2str(k) '.mat'],'squared_error')
            save(['./trial' num2str(folder_index) '/mmse/Sk_' num2str(num_save) 'it' num2str(k) '.mat'],'Sk')
            save(['./trial' num2str(folder_index) '/mmse/wik_' num2str(num_save) 'it' num2str(k) '.mat'],'wik')
            save(['./trial' num2str(folder_index) '/mmse/thetaik_' num2str(num_save) 'it' num2str(k) '.mat'],'thetaik')
            save(['./trial' num2str(folder_index) '/mmse/thetahatv_' num2str(num_save) 'it' num2str(k) '.mat'],'thetahatv')
            save(['./trial' num2str(folder_index) '/mmse/binwidth_' num2str(num_save) 'it' num2str(k) '.mat'],'binwidth')
        elseif design_type == 4
            save(['./trial' num2str(folder_index) '/inverse/se_' num2str(num_save) 'it' num2str(k) '.mat'],'squared_error')
            save(['./trial' num2str(folder_index) '/inverse/Sk_' num2str(num_save) 'it' num2str(k) '.mat'],'Sk')
            save(['./trial' num2str(folder_index) '/inverse/wik_' num2str(num_save) 'it' num2str(k) '.mat'],'wik')
            save(['./trial' num2str(folder_index) '/inverse/thetaik_' num2str(num_save) 'it' num2str(k) '.mat'],'thetaik')
            save(['./trial' num2str(folder_index) '/inverse/thetahatv_' num2str(num_save) 'it' num2str(k) '.mat'],'thetahatv')
            save(['./trial' num2str(folder_index) '/inverse/binwidth_' num2str(num_save) 'it' num2str(k) '.mat'],'binwidth')
        %elseif design_type == 5
        %    save(['./trial' num2str(folder_index) '/rmb/se_' num2str(num_save) 'it' num2str(k) '.mat'],'squared_error')
        %    save(['./trial' num2str(folder_index) '/rmb/Sk_' num2str(num_save) 'it' num2str(k) '.mat'],'Sk')
        %    save(['./trial' num2str(folder_index) '/rmb/wik_' num2str(num_save) 'it' num2str(k) '.mat'],'wik')
        %    save(['./trial' num2str(folder_index) '/rmb/thetaik_' num2str(num_save) 'it' num2str(k) '.mat'],'thetaik')
        %    save(['./trial' num2str(folder_index) '/rmb/thetahatv_' num2str(num_save) 'it' num2str(k) '.mat'],'thetahatv')
        %    save(['./trial' num2str(folder_index) '/rmb/binwidth_' num2str(num_save) 'it' num2str(k) '.mat'],'binwidth')
        elseif design_type == 6
            save(['./trial' num2str(folder_index) '/adapt_dd/se_' num2str(num_save) 'it' num2str(k) '.mat'],'squared_error')
            save(['./trial' num2str(folder_index) '/adapt_dd/Sk_' num2str(num_save) 'it' num2str(k) '.mat'],'Sk')
            save(['./trial' num2str(folder_index) '/adapt_dd/wik_' num2str(num_save) 'it' num2str(k) '.mat'],'wik')
            save(['./trial' num2str(folder_index) '/adapt_dd/thetaik_' num2str(num_save) 'it' num2str(k) '.mat'],'thetaik')
            save(['./trial' num2str(folder_index) '/adapt_dd/thetahatv_' num2str(num_save) 'it' num2str(k) '.mat'],'thetahatv')
            save(['./trial' num2str(folder_index) '/adapt_dd/binwidth_' num2str(num_save) 'it' num2str(k) '.mat'],'binwidth')
        elseif design_type == 7
            save(['./trial' num2str(folder_index) '/rmbrand/se_' num2str(num_save) 'it' num2str(k) '.mat'],'squared_error')
            save(['./trial' num2str(folder_index) '/rmbrand/Sk_' num2str(num_save) 'it' num2str(k) '.mat'],'Sk')
            save(['./trial' num2str(folder_index) '/rmbrand/wik_' num2str(num_save) 'it' num2str(k) '.mat'],'wik')
            save(['./trial' num2str(folder_index) '/rmbrand/thetaik_' num2str(num_save) 'it' num2str(k) '.mat'],'thetaik')
            save(['./trial' num2str(folder_index) '/rmbrand/thetahatv_' num2str(num_save) 'it' num2str(k) '.mat'],'thetahatv')
            save(['./trial' num2str(folder_index) '/rmbrand/binwidth_' num2str(num_save) 'it' num2str(k) '.mat'],'binwidth')
        end
    end
    %toc
    %tout = [tout;toc];
end

if globals.printops
    [voidvar] = plotter(wik,thetaik,Sk,binwidth);
end

root_squared_error = squared_error.^0.5;

