% Initialization
nx = 5;
ny = nx;

ux0 = 0;
ux1 = 0;
uy0 = 0;
uy1 = 1;

% Initialize Au = b
u = zeros(nx * ny, 1);
b = zeros(nx * ny, 1);
%A = zeros(nx * ny);
ind = 0;

is = 0;
for j = 1 : ny
    is = is + 2;
end

for i = 1 : nx
    is = is + 2;
end
for ri = 2 : nx-1
    for rj = 2 : nx-1
        is = is + 5;
    end
end

iv = zeros(is, 1);
jv = zeros(is, 1);
val = zeros(is, 1);

% Assign Boundary
for j = 1 : ny
    i = 1;
    ij = getIJ(i, j, nx);
    u(ij) = ux0;
    b(ij) = ux0;
    ind = ind + 1;
    iv(ind) = ij;
    jv(ind) = ij;
    val(ind) = 1;
    
    i = nx;
    ij = getIJ(i, j, nx);
    u(ij) = ux1;
    b(ij) = ux1;
    ind = ind + 1;
    iv(ind) = ij;
    jv(ind) = ij;
    val(ind) = 1;
end

for i = 1 : nx
    j = 1;
    ij = getIJ(i, j, nx);
    u(ij) = uy0;
    b(ij) = uy0;
    ind = ind + 1;
    iv(ind) = ij;
    jv(ind) = ij;
    val(ind) = 1;
    
    j = ny;
    ij = getIJ(i, j, nx);
    u(ij) = uy1;
    b(ij) = uy1;
    ind = ind + 1;
    iv(ind) = ij;
    jv(ind) = ij; 
    val(ind) = 1;
end


% Construct A matrix
for ri = 2 : nx-1
    for rj = 2 : nx-1
        rij = getIJ(ri, rj, nx);
        
        ci = ri + 1;
        cj = rj;
        cij = getIJ(ci, cj, nx);
        ind = ind + 1;
        iv(ind) = rij;
        jv(ind) = cij;
        val(ind) = 1;

        ci = ri - 1;
        cj = rj;
        cij = getIJ(ci, cj, nx);
        ind = ind + 1;
        iv(ind) = rij;
        jv(ind) = cij;       
        val(ind) = 1;
        
        ci = ri;
        cj = rj + 1;
        cij = getIJ(ci, cj, nx);
        ind = ind + 1;
        iv(ind) = rij;
        jv(ind) = cij;
        val(ind) = 1;
        
        ci = ri;
        cj = rj - 1;
        cij = getIJ(ci, cj, nx);
        ind = ind + 1;
        iv(ind) = rij;
        jv(ind) = cij;
        val(ind) = 1;

        ci = ri;
        cj = rj;
        cij = getIJ(ci, cj, nx);
        ind = ind + 1;
        iv(ind) = rij;
        jv(ind) = cij;
        val(ind) = -4;
    end
end

As = sparse(iv, jv, val, nx * ny, nx * ny);
usol = As \ b;
usol = reshape(usol, [nx, ny]);
close all;
figure
contourf(0 : 1 / (nx - 1) : 1, 0 : 1 / (ny - 1) : 1, usol)
figure
spy(As)