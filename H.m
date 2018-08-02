function [Hop] = H(theta,globals)
%By convention consider left-most element of the array to be element 1 for
%transmit and receive

NT = globals.NT;
NR = globals.NR;
Ntargets = globals.Ntargets;

at = zeros(NT,1);
ar = zeros(NR,1);
Hop = zeros(NR,NT);

for i = 1:Ntargets
    phi = theta.phi(1,i);
    for nt=1:NT
        at(nt,1) = exp(-(nt-1)*pi*sqrt(-1)*cosd(90+phi));
    end
    for nr=1:NR
        ar(nr,1) = exp(-(nr-1)*pi*sqrt(-1)*cosd(90+phi));
    end
    Hop = Hop + theta.alpham(1,i)*(cosd(theta.alphaa(1,i)) + sqrt(-1)*sind(theta.alphaa(1,i)))*ar*transpose(at);
end

end

 