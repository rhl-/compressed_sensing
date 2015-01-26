function A = matrix_sample(m,n,meas,mat)
% A = matrix_sample(m,n,type)
% Matrix sampling
% Inputs
% m number of measurements
% n dimension of n times n matrices
% meas string: ?Entry?,?Perm?,?RSPerm?,?CSPerm?,?RGPerm?,?CGPerm?,?RDirac?,?CDirac?,?RGauss?,?CGauss?
% mat string: 'RPSD', 'RSYM', 'HPSD', 'HERM'
% Outputs
% A an array m by n^2, containing the sampling matrices as row vectors

switch meas
    case 'Entry'
        entryMat();
    case 'Perm'
        permMat();
    case 'RSPerm'
        rspermMat();
    case 'CSPerm'
        cspermMat();
    case 'RGPerm'
        rgpermMat();
    case 'CGPerm'
        cgpermMat();
    case 'RDirac'
        I = eye(2);
        wx = [0,1;1,0];
        wy = [0,-1;1,0];
        wz = [1,0;0,-1];
        mySet ={I,wx,wy,wz};
        diracMat();
    case 'CDirac'
        I = eye(2);
        sx = [0,1;1,0];
        sy = [0,-sqrt(-1);sqrt(-1),0];
        sz = [1,0;0,-1];
        mySet ={I,sx,sy,sz};
        diracMat();
    case 'RGauss'
        rgaussMat();
    case 'CGauss'
        cgaussMat();
    otherwise
        disp('Error in matrix_sample: unrecognized type!');
        
    if size(A) ~= size(unique(A,'rows'))
        A = matrix_sample(A,m,n,meas,mat);
    end
end


    function entryMat()
        %Sample a random coordinate for each matrix from the upper triangle
        %randints = randi(n^2,m,1);
        %randintsI = randi(n,m,1);
        %randintsJ = randi(n,m,1);
        %randints = sub2ind([n,n],randintsI,randintsJ);
        population = find(triu(ones(n)));
        dests = randsample(population,m,'false');
        A = sparse(1:m,dests,1,m,n^2);
    end

    function permMat()
        %Each matrix is a random permutation scaled by 1/sqrt(n)
        A = zeros(m,n^2);
        for i = 1:m
            perm = randperm(n);
            Ai = sparse(1:n,perm,1,n,n);
            A(i,:) = Ai(:);
        end
        A = A / sqrt(n);
    end
    function rspermMat()
        %Add a sign to the permutation
        permMat();
        signs = (rand(m,n^2) < 0.5);
        A = A.*signs;
    end
    function cspermMat()
        %Add a "complex sign" to the permutation
        permMat();
        randomMat = rand(m,n^2);
        % a quarter get the identity, the rest get i
        signs = sqrt(-1) * (randomMat > 0.25);
        A = A.*signs;
        % a quarter have the identity, another quarter are i, the rest get i again
        signs = signs.*(randomMat > 0.5);
        A = A.*signs;
        % a quarter have the identity, a quarter have i, a quarter have -1,
        % the rest get i again
        signs = signs.*(randomMat > 0.75);
        A = A.*signs;
    end
    function rgpermMat()
        %Modulate the permutation with a normal
        permMat();
        A = A.*randn(m,n^2);
    end
    function cgpermMat()
        %Modulate the permutation with a normal complex
        permMat();
        A = 1/sqrt(2)*A.*(randn(m,n^2) + sqrt(-1)*randn(m,n^2));
    end
    function diracMat()
        %Ughhh kronecker
        J = log2(n);
        assert(J == round(J));

        A = zeros(m,n^2);
        % Every number from 1 to n^2 is identified with a matrix in the
        % ensemble
        choices = randperm(n^2);
        ctr = 1;
        for i = 1:m
            parity = 1;
            while parity == 1
                %Keep looking for 1 with even parity
                crnt = choices(ctr);
                [seq,parity] = int2sequence(crnt,J);
                ctr = ctr + 1;
                if mat(1) == 'H'
                    break;
                end
            end
            Ai = 1;
            for j = 1:J
                Ai = kron(Ai,mySet{seq(j)+1});
            end
            
            if mat(1) ~= 'H'
                assert (norm(Ai +Ai.') > 1e-8)
            end
            %assert (norm(Ai + Ai') > 1e-8);
            A(i,:) = Ai(:);
        end
%         for j = 1:m
%             flag = false;
%             while(~flag)
%                     Ai = 1;
%                 for i = 1:J
%                     b = mySet{randi(4,1,1)};
%                     Ai = kron(Ai,b);
%                 end
% %                 if norm(Ai - Ai.') < 1e-5 || mat(1) == 'H'
% %                     % Matrix is symmetric we are good
% %                     % Or we are complex
% %                     flag = true;
% %                 end
%             end
%             A(j,:) = Ai(:);
%         end
        A = A/sqrt(n);
        A = sparse(A);
    end
%     function cdiracMat()
%         %Ughhh kronecker
%         J = log2(n);
%         assert(J == round(J));
%         I = eye(2);
%         sx = [0,1;1,0];
%         sy = [0,-sqrt(-1);sqrt(-1),0];
%         sz = [1,0;0,-1];
%         mySet ={I,sx,sy,sz};
% 
%         A = zeros(m,n^2);
%         for j = 1:m
%             flag = false;
%             while(~flag)
%                 Ai = 1;
%                 for i = 1:J
%                    b = mySet{randi(4,1,1)};
%                    Ai = kron(Ai,b);
%                 end
%                 if norm(Ai - Ai.') < 1e-5 || mat(1) == 'H'
%                     % Matrix is symmetric we are good
%                     % Or we are complex
%                     flag = true;
%                 end
%             end
%             A(j,:) = Ai(:);
%         end
%         A = A/sqrt(n);
%         A = sparse(A);
%         
%     end
    function rgaussMat()
        A = 1/n * randn(m,n^2);
    end
    function cgaussMat()
        A = 1/n * (randn(m,n^2) + sqrt(-1) * randn(m,n^2))/sqrt(2);
    end
    
    function [seq,parity] = int2sequence(int,J)
        %Parity will count the number of times the antisymmetric matrix
        %appears in this measurement
        parity = 0;
        seq = zeros(J,1);
        for i = 1:J
            % pop off the "tail" of the number, which is the integer
            % encoding of the index into the 4-character alphabet for the
            % current "b"
            tail = mod(int,4);
            seq(i) = tail;
            if tail == 2 %0,1,2,3 == I,x,y,z so 2 is y
                % If we chose y, flip the parity bit
                parity = 1-parity;
            end
            % Get ready to grab the next letter encoding
            int = (int - tail) / 4;            
        end
    end

end