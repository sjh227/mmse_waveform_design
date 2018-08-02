function [theta_out] = param_model(theta_in,globals)

Ntargets = globals.Ntargets; %for up to 4 targets

if size(theta_in,1) == 0
    theta_out.phi = [-40,20,62,80]; %OLD [-20,...]
    theta_out.alphaa = [46,329,228,35]; %randomly selected
    theta_out.alpham = [1,1,1,1];
    theta_out.phi = theta_out.phi(1,1:Ntargets);
    theta_out.alphaa = theta_out.alphaa(1,1:Ntargets);
    theta_out.alpham = theta_out.alpham(1,1:Ntargets);
    theta_out.alpham = theta_out.alpham/(sqrt(theta_out.alpham*theta_out.alpham'));
else
    theta_out = theta_in; %can use for moving targets
end
