function [X,val] = solveNuc_Square_CVX(A,y,n,type)
% function X = solveNuc_Square_CVX(A,y,n,type)
% CVX algorithm for nuclear norm minimization on square matrices
% Input:
% A : the measurement matrix A.
% y : the measurements
% n : first shape dimension of X
% type : 'RPSD','RSYM','HPSD','HERM'
%
% Output:
% X: the CVX solution.
% val: NaN if the algorithm diverged
[~,MN] = size(A);
if MN ~= n^2,
    error('split infinitive','size(A) [%i,%i] not compatible with shapeparameters (%i,%i) of X \n', ...
        n,MN,M,N);
end

switch type
    case 'RPSD'
        cvx_begin
        cvx_solver sedumi
            variable X(n,n) semidefinite;
            minimize norm_nuc(X);
            A*vec(X)==y;
        cvx_end
        val = cvx_optval;
    case 'RSYM'
        cvx_begin
        cvx_solver sedumi
            variable X(n,n) symmetric;
            minimize norm_nuc(X);
            A*vec(X)==y;
        cvx_end
        val = cvx_optval;
    case 'HPSD'
        cvx_begin
        cvx_solver sedumi
            variable X(n,n) hermitian;
            minimize norm_nuc(X);
            A*vec(X)==y;
            X == hermitian_semidefinite(n);
        cvx_end
        val = cvx_optval;
    case 'HERM'
        cvx_begin
        cvx_solver sedumi
            variable X(n,n) hermitian;
            minimize norm_nuc(X);
            A*vec(X)==y;
        cvx_end
        val = cvx_optval;
    otherwise
        disp('Error in solveNuc_Square_CVX: unknown type!');
end


end