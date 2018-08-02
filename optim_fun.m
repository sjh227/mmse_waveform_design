function[cost] = optim_fun(dR,gradR)
%cost function

L=5;

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


cost = real(trace(gradR*dRm'));

end

