function v = svec3(A)
    m = size(A,1);
    v = zeros(m*(m+1)/2,1);
    k = 1;
    for j = 1:m
        for i = 1:j
            if j > i
                v(k) = sqrt(2)*A(i,j);
                k = k + 1;
            else % j == i
                v(k) = A(i,i);
                k = k + 1;
            end
        end
    end
end