function[c,ceq] = pd_fun(dR,Rs)



dRm = [];
L = sqrt(size(dR,1)/2);


for i=1:L
    dRtemp = dR((i-1)*L+1:i*L,1) + sqrt(-1)*dR(L^2+(i-1)*L+1:L^2+i*L,1);
    dRm = [dRm,dRtemp];
end



[~,ceq] = chol(Rs+dRm);


c=[];