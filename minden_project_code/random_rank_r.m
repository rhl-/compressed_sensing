function X0 = random_rank_r(n,r,type)
% X0 = random_rank_r(n,r,type)
% Matrix sampling
% Inputs
% n dimension of n times n matrices
% r rank of matrix 1 <= r <= n
% type string: ?RPSD?,?RSYM?,?HPSD?,?HERM?
% Outputs
% X0 an array n by n of rank r

    switch type
        
        case 'RPSD'
            rpsdMat();
        case 'RSYM'
            rsymMat();
        case 'HPSD'
            hpsdMat();
        case 'HERM'
            hermMat();
        otherwise
            disp('Error in random_rank_r: unknown type');
    end

    function rpsdMat()
         raw = randn(n,r);
        [U,~,~] = qr(raw,0);
        if( size(U,2) ~= r ),
            error('TOO MUCH HAAR','Expected %i cols but found U (%i) ', ...
                r,size(U,2));
        end
        X0 = U * U';
    end

    function rsymMat()
         raw = randn(n,r);
        [U,~,~] = qr(raw,0);
        if( size(U,2) ~= r ),
            error('TOO MUCH HAAR','Expected %i cols but found U (%i) ', ...
                r,size(U,2));
        end
        X0 = U *diag(sign(randn(r,1)))* U';
    end

    function hpsdMat()
        raw = randn(n,r) + sqrt(-1.).* randn(n,r);
        [U,~,~] = qr(raw,0);
        if( size(U,2) ~= r ),
            error('TOO MUCH HAAR','Expected %i cols but found U (%i) ', ...
                r,size(U,2));
        end
        X0 = U * U';
    end

    function hermMat()
        raw = randn(n,r) + sqrt(-1.).* randn(n,r);
        [U,~,~] = qr(raw,0);
        if( size(U,2) ~= r ),
            disp(size(U))
            disp(r)
            error('TOO MUCH HAAR','Expected %i cols but found U (%i) ', ...
                r,size(U,2));
        end
        X0 = U *diag(sign(randn(r,1)))* U';
    end

end