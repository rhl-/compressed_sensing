%Unit tests for matrix_sample
clear;clc;

m = 128;
n = 128;
disp('Testing Entry...');
tic();
A = matrix_sample(m,n,'Entry','HPSD');
toc
% Test: assert that there is only one nonzero
for i = 1:m
    Ai = A(i,:);
    assert(nnz(Ai) == 1);
end


disp('Testing Perm...');
tic;
A = matrix_sample(m,n,'Perm','HPSD');
toc
% Test: assert that there is only one nonzero per row/column
for i = 1:m
    Ai = reshape(A(i,:),n,n);
    for k = 1:n
        assert(nnz(Ai(k,:)) == 1);
        assert(nnz(Ai(:,k)) == 1);
    end
end


disp('Testing RSPerm...');
tic();
A = matrix_sample(m,n,'RSPerm','HPSD');
toc()
% Test: since perm test worked, just need to check signs
B = n*A.*A;

assert(abs(sum(sum(B)) - nnz(B)) < 1e-4);


disp('Testing CSPerm...');
tic();
A = matrix_sample(m,n,'CSPerm','HPSD');
toc
% Test: since perm test worked, just need to check complex signs
B = n^2*A.*A.*A.*A;
assert(abs(sum(sum(B)) - nnz(B)) < 1e-4);


disp('Testing RGPerm...');
tic();
A = matrix_sample(m,n,'RGPerm','HPSD');
toc
% Test: since perm test worked, just need to check that statistics of nonzeros are correct;
B = A(:);
B=B(B~=0);
disp('This should be close to zero');
disp(mean(B));
disp('This should be close to zero');
disp(abs(mean(B.^2) - 1/n))


disp('Testing CGPerm...');
tic();
A = matrix_sample(m,n,'CGPerm','HPSD');
toc
% Test: since perm test worked, just need to check that statistics of nonzeros are correct;
B = A(:);
B=B(B~=0);
disp('This should be close to zero');
disp(abs(mean(B)));
disp('This should be close to zero');
disp(abs(mean(abs(B).^2) - 1/n))


disp('Testing RDirac...');
tic();
A = matrix_sample(m,n,'RDirac','HPSD');
toc
% Test: check that fro norm is right
disp('This should be 1');
disp(norm(A(1,:),2))

disp('Testing CDirac...');
tic();
A = matrix_sample(m,n,'CDirac','HPSD');
toc
% Test: check that fro norm is right
disp('This should be 1');
disp(norm(A(1,:),2))

Adir = A;




disp('Testing RGauss...');
tic();
A = matrix_sample(m,n,'RGauss','HPSD');
toc
% Test: Just need to check that statistics of nonzeros are correct;
B = A(:);
disp('This should be close to zero');
disp(abs(mean(B)));
disp('This should be close to zero');
disp(abs(mean(abs(B).^2) - 1/n^2))


disp('Testing CGauss...');
tic();
A = matrix_sample(m,n,'CGauss','HPSD');
toc
B = A(:);
disp('This should be close to zero');
disp(abs(mean(B)));
disp('This should be close to zero');
disp(abs(mean(abs(B).^2) - 1/n^2))
