function [lambda,w] = mineig(A)
    % [w,lambda] = eig(A);
    % [lambda,jmin] = min(diag(lambda));
    % w = w(:,jmin);
    % % eig is slower than eigs
    % % run eig_time.m
    % https://www.mathworks.com/matlabcentral/answers/20522-find-max-min-eigenvalue-of-a-symmetric-matrix
    [w,lambda] = eigs(A,1,'smallestreal');
end 