%random_rank_r testing
clear;clc;
n = 100;
r = 25;


%RPSD
A = random_rank_r(n,r,'RPSD');
lambda = eig(A);
assert(min(lambda) > -1e-3);
assert(norm(A-A','fro') == 0);

%RSYM
A = random_rank_r(n,r,'RSYM');
lambda = eig(A);
lambda = lambda(abs(lambda) > 1e-3);
assert(max(abs(abs(lambda)-1)) < 1e-3);
assert(norm(A-A','fro') == 0);

%HPSD
A = random_rank_r(n,r,'HPSD');
lambda = eig(A);
assert(min(lambda) > -1e-3);
assert(norm(A-A','fro') == 0);

%HERM
A = random_rank_r(n,r,'HERM');
lambda = eig(A);
lambda = lambda(abs(lambda) > 1e-3);
assert(max(abs(abs(lambda)-1)) < 1e-3);
assert(norm(A-A','fro') == 0);

