function elem = gen_elem(coord,rem_type)
    
    nn = size(coord,1);

    % all [a b] with a and b in 1:nn 
    [n1,n2] = meshgrid(1:nn,1:nn);
    elem = [n1(:) n2(:)];
    
    % Remove [a a]
    elem(elem(:,1)==elem(:,2),:) = [];
    
    % Remove [a b] if a>b
    elem = sort(elem,2);
    elem = unique(elem,'rows');
    
    % Remove bars longer than distmax
    n1 = elem(:,1);
    n2 = elem(:,2);
    c1x = coord(n1,1);
    c1y = coord(n1,2);
    c1z = coord(n1,3);
    c2x = coord(n2,1);
    c2y = coord(n2,2);
    c2z = coord(n2,3);
    dx = c2x-c1x;
    dy = c2y-c1y;
    dz = c2z-c1z;
    lengths = sqrt(dx.^2 + dy.^2 + dz.^2);
    if rem_type == 1
        distmax = max(lengths); % this way nothing is removed
    elseif rem_type == 2
        distmax = sqrt(2);
    end
    elem(lengths>distmax,:) = [];
    
    % Remove elements longer and in the same line than smaller ones
    nel = size(elem,1);
    elem2rem = false(nel,1);
    for i = 1:nel
        element = elem(i,:);
        n1 = element(1);
        n2 = element(2);
        v = [coord(n2,1)-coord(n1,1)
             coord(n2,2)-coord(n1,2)
             coord(n2,3)-coord(n1,3)];
        for n3 = 1:nn
            if n3 ~= n1 && n3 ~= n2
                w = [coord(n3,1)-coord(n1,1)
                     coord(n3,2)-coord(n1,2)
                     coord(n3,3)-coord(n1,3)];
                if norm(cross(v,w)) < 1e-6 % v and w are parallel
                    [~,j] = max(abs(v));
                    if 0 < w(j)/v(j) && w(j)/v(j) < 1
                        elem2rem(i) = true;
                    end
                end
            end
        end
    end
    elem(elem2rem,:) = [];
end