function coord = gen_coord(nx,ny,nz)
    index = 1;
    coord = zeros(nx*ny*nz,3);
    for k = 1:nz
        for j = 1:ny
            for i = 1:nx
                coord(index,:) = [i-1 j-1 k-1];
                index = index + 1;
            end
        end
    end
end