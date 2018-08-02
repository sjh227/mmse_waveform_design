function[T, Dinv] = DandT(thetai,wi,globals)
%normalised direction according to BRMB design method

Rs = globals.Rs;
J = globals.J;
NP = size(wi,1);
L = globals.L;

%generate test points: matrix T of test points size: (1,J) for one target
sdtheta = sqrt(wi'*(thetai.phi.*thetai.phi) - (wi'*thetai.phi)^2);

T = [];
j=1;
while j <= J
    T_temp = min(89,max(-90,round(normrnd(0,sdtheta))));
    newT = ismember(T_temp,T);
    if ~newT
        T = [T,T_temp];
        j = j+1;
    end
end



%generate D
% D = zeros(J);
% for i=1:J
%     for j=1:J
%         pm = 1;
%         for np = 1:NP
%             theta = thetai.phi(np);
%             thetaii.phi = thetai.phi(np) + T(i);
%             thetaii.alphaa = thetai.alphaa(np);
%             thetaii.alpham = thetai.alpham(np);
%             thetaj.phi = thetai.phi(np) + T(j);
%             thetaj.alphaa = thetai.alphaa(np);
%             thetaj.alpham = thetai.alpham(np);
%             theta0.phi = thetai.phi(np);
%             theta0.alphaa = thetai.alphaa(np);
%             theta0.alpham = thetai.alpham(np);
%             if -91 < theta + T(i) && -91 < theta + T(j) && 90 > theta + T(i) && 90 > theta + T(j) 
%                 A = (H(thetaii,globals) - H(theta0,globals))'*(H(thetaj ,globals) - H(theta0,globals)); %noting that noise covariance is defined as the identity
%                 D(i,j) = D(i,j) + (exp(2*L*trace(real(A*Rs)))*wi(91+theta + T(i))*wi(91+theta + T(j))/wi(91+theta));       
%             else
%                 pm = pm - wi(np);
%             end
%             if pm > 0
%                 D(i,j) = D(i,j)/pm;
%             else
%                 D(i,j) = 0;
%             end
%         end
%     end
% end
% 
% 
% %generate G
% Dinv = (D-ones(J))^-1;

D = zeros(J);
for i=1:J
    for j=1:J
        pm = 1;
        for np = 1:NP
            theta = thetai.phi(np);
            thetaii.phi = thetai.phi(np) + T(i);
            thetaii.alphaa = thetai.alphaa(np);
            thetaii.alpham = thetai.alpham(np);
            thetaj.phi = thetai.phi(np) + T(j);
            thetaj.alphaa = thetai.alphaa(np);
            thetaj.alpham = thetai.alpham(np);
            theta0.phi = thetai.phi(np);
            theta0.alphaa = thetai.alphaa(np);
            theta0.alpham = thetai.alpham(np);
            if -91 < theta + T(i) && -91 < theta + T(j) && 90 > theta + T(i) && 90 > theta + T(j) 
                A = (H(thetaii,globals) - H(theta0,globals))'*(H(thetaj ,globals) - H(theta0,globals)); %noting that noise covariance is defined as the identity
                D(i,j) = D(i,j) + (exp(2*L*trace(real(A*Rs)))*wi(91+theta + T(i))*wi(91+theta + T(j))/wi(91+theta));       
            else
                pm = pm - wi(np);
            end
        end
        if pm > 0
            D(i,j) = D(i,j)/pm;
        else
            D(i,j) = 0;
        end
    end
end


%generate G
Dinv = (D-ones(J))^-1;


