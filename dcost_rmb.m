function[dir,Dinv] = dcost_rmb(thetai,wi,T,Dinv,globals)
%normalised direction according to BRMB design method

Rs = globals.Rs;
NT = globals.NT;
J = globals.J;
NP = size(wi,1);
L = globals.L;

Rsvec = [];
for i=1:size(Rs,2)
    Rsvec = [Rsvec;Rs(:,i)];
end

Aeq = zeros(1,(L*(L+1)/2)+2*(L*(L-1)/2));
for i=1:L
    Aeq(1,i) = 1; %encode trace constraint
    %Aeq(2,L^2+i) = 1; %encode trace constraint
end
beq = [0];
%generate test points: matrix T of test points size: (1,J) for one target
% sdtheta = sqrt(wi'*(thetai.phi.*thetai.phi) - (wi'*thetai.phi)^2);
% 
% T = [];
% j=1;
% while j <= J
%     T_temp = min(89,max(-90,round(normrnd(0,sdtheta))));
%     newT = ismember(T_temp,T);
%     if ~newT
%         T = [T,T_temp];
%         j = j+1;
%     end
% end


%generate L
Lk = zeros(NT*J*J,NT);
for i=1:NT
    for j=1:NT
        Lij = zeros(J);
        for l=1:J
            for m =1:J
                ind1 = l;
                ind2 =  m;
                pm = 1;
                for np = 1:NP
                    theta = thetai.phi(np);
                    thetal.phi = thetai.phi(np) + T(l);
                    thetal.alphaa = thetai.alphaa(np);
                    thetal.alpham = thetai.alpham(np);
                    thetam.phi = thetai.phi(np) + T(m);
                    thetam.alphaa = thetai.alphaa(np);
                    thetam.alpham = thetai.alpham(np);
                    theta0.phi = thetai.phi(np);
                    theta0.alphaa = thetai.alphaa(np);
                    theta0.alpham = thetai.alpham(np);
                    if -91 < theta + T(l) && -91 < theta + T(m) && 90 > theta + T(l) && 90 > theta + T(m) 
                        A = (H(thetal,globals) - H(theta0,globals))'*(H(thetam ,globals) - H(theta0,globals)); %noting that noise covariance is defined as the identity
                        Ael = A(j,i); %i,j swapped as transpose
                        Lij(ind1,ind2) = Lij(ind1,ind2) + L*(exp(2*L*trace(real(A*Rs)))*Ael*wi(91+theta + T(l))*wi(91+theta + T(m))/wi(91+theta));       
                    else
                        pm = pm - wi(np);
                    end
                end
                if pm > 0
                    Lij(ind1,ind2) = Lij(ind1,ind2)/pm;
                else
                    Lij(ind1,ind2) = 0;
                end
            end
        end
        vecL = [];
        for vi=1:size(Lij,2)
            vecL = [vecL;Lij(:,vi)];
        end
        Lk((i-1)*(J^2)+1:i*(J^2),j) = vecL;
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
%         end
%         if pm > 0
%             D(i,j) = D(i,j)/pm;
%         else
%             D(i,j) = 0;
%         end
%     end
% end
% 
% 
% %generate G
% Dinv = (D-ones(J))^-1;
G = Dinv*transpose(T)*T*Dinv;

%evaluate gradR
vecG = [];
for i=1:size(G,2)
    vecG = [vecG;G(:,i)];
end
gradR = -kron(eye(NT),transpose(vecG))*Lk;

%code in the optimisation problem
dR0 = zeros((L*(L+1)/2)+2*(L*(L-1)/2),1); %correct size for dR0
validinit = optim_fun(dR0,gradR);
if ~isnan(validinit) && validinit ~= inf && validinit ~= -inf
    dR = fmincon(@(dR)optim_fun(dR,gradR),dR0,[],[],Aeq,beq); %,[],[],@(dR)pd_fun(dR,Rs));
else
    dR = dR0;
end

% dRm = [];
% 
% 
% for i=1:L
%     dRtemp = dR((i-1)*L+1:i*L,1) + sqrt(-1)*dR(L^2+(i-1)*L+1:L^2+i*L,1);
%     dRm = [dRm,dRtemp];
% end

dRm = zeros(L,L);

relocmat = [1,6,7,8,9;0,2,10,11,12;0,0,3,13,14;0,0,0,4,15;0,0,0,0,5];
imlocmat = [0,16,17,18,19;0,0,20,21,22;0,0,0,23,24;0,0,0,0,25;0,0,0,0,0];

for i=1:L
    for j=i:L
        if i==j
            dRm(i,j) = dR(i,1); %leading diagonal
        else
            repart = dR(relocmat(i,j),1);
            impart = dR(imlocmat(i,j),1);
            dRm(i,j) = repart + sqrt(-1)*impart;
            dRm(j,i) = repart - sqrt(-1)*impart;
        end
    end
end

dir = dRm/(sqrt(trace(dRm'*dRm)));
% dir = [];
% dirsize = sqrt(size(dR,1));
% for i=1:dirsize
%     dir = [dir,dR((i-1)*dirsize+1:i*dirsize)];
% end

