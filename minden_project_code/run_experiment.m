function result = run_experiment(FID, n, nMC,v_list, Meas, Mat)
% result = run_experiment(FID, n, nMC,v_list, Meas, Mat)
% Input:
%   FID - the matlab file identifier of output file
%   n - the dimension of the array being recovered
%   nMC - number of Monte Carlo trials
%   v_list - n_list (m,r) pairs to try
%   Meas - string naming the measurement ensemble
%   Mat - string naming the matrix ensemble
% Output:
%   result - n_list by 6 array (n,r,m,nMC,err0,err1)

tol = 1e-2;
n_list = size(v_list,1);
result = zeros(n_list,6);
nErr = 2;

for i_list = 1:n_list
    m = v_list(i_list,1);
    r = v_list(i_list,2);
    res = zeros(nMC,nErr);
    for mc = 1:nMC
        [y,A,X0] = problem_instance(n,r,m,Meas,Mat);
        [X1,val] = solveNuc_Square_CVX(A,y,n,Mat);
        assert(~isnan(val));
        obj0 = norm(X0,'fro');
        err1 = norm(X1-X0,'fro') ./ obj0;
        err0 = err1 < tol;
        res(mc,:) = [err0, err1];
        my_str = sprintf('MC Sample(%i) =[ %i %i %i %i %0.3f %0.3f]\n',mc,n,r,m,mc,...
            res(mc,1),...
            res(mc,2));
        disp(my_str);
        fprintf(FID,my_str);
    end
    result(i_list,:) = [n r m nMC mean(res) ];
    fprintf(FID,'result(%i,:)= [%i %i %i %i ', i_list,n,r,m,nMC);
    fprintf(FID,'%6.3f %6.3f]\n',...
       result(i_list,5),...
       result(i_list,6));
end
end