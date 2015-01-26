function [y,A,X0] = problem_instance(n,r,m,Meas,Mat)
% [y,A,X0] = problem_instance(n,r,Meas,Mat,m)
% Generate a random problem instance

%tol = 1e-10;
A = matrix_sample(m,n,Meas,Mat);
X0 = random_rank_r(n,r,Mat);
y = A*vec(X0);
%y(abs(y) < tol) = 0;


end