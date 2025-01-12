% elasticity parameters
E = 1e5;       % Young's modulus
vv = 0.3;      % Poisson's ratio
lamda = vv*E/((1+vv)*(1-2*vv));
miu = E/(2*(1+vv));
D = [lamda+2*miu, lamda, 0; lamda, lamda+2*miu, 0; 0, 0, miu]; % isotropic elastic material

% force function
f_x = @(x, y) 0; % body force in x-direction
f_y = @(x, y) 0; % body force in y-direction

% traction on boundary (Neumann BCs)
traction = @(x, y) [10, 0]; % for example, 10 kPa in x-direction, 0 in y-direction

% quadrature rule
n_int_xi = 3;
n_int_eta = 3;
n_int = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% mesh generation
nodes = msh.POS(:, 1: 2);
triangular = msh.TRIANGLES(:, 1: 3);
% imported from gmsh

n_np_x = length(nodes(:, 1)) ^ 0.5;  % number of nodal points in x-direction
n_np_y = length(nodes(:, 2)) ^ 0.5;  % number of nodal points in y-direction
n_np = n_np_x * n_np_y;  % total number of nodal points

n_en = 4; % number of nodes in each element
n_el_x = n_np_x - 1; % number of elements in x-direction
n_el_y = n_np_y - 1; % number of elements in y-direction
n_el = n_el_x * n_el_y; % total number of elements

% generate the nodal coordinates
xx_coor = zeros(n_np, 1);
for i = 1: n_np
    xx_coor(i) = x_coor(nodes, i);
end
yy_coor = zeros(n_np, 1);
for i = 1: n_np
    yy_coor(i) = y_coor(nodes, i);
end

% construct data structure
% IEN array
IEN = zeros(n_el, n_en);
for ee = 1: n_el
    for aa = 1: n_en
        IEN(ee, aa) = IEN_quadrilateral(aa, ee, triangular);
    end
end

% convert the physical position into global node number
% because node number distribution is irregular
index = zeros(n_np, 1);
for i=1:n_el_y
    for j=1:n_el_x
        k1 = (i-1)*n_np_x+j;
        k2 = (i-1)*n_el_x+j;
        index(k1) = triangular(2*k2-1, 1);
    end
    index((i-1)*n_np_x+n_np_x) = triangular(2*k2-1, 3);
end
aa=(n_np_y-1)*n_np_x;
bb=(n_el_y-1)*n_el_x;
for i=1:n_el_x
    index(aa + i) = triangular(2*(bb+i)-1, 2);
end
index(n_np) = triangular(2*n_el, 3);

% ID array (include all DOFs for x = 1)
ID = zeros(n_np, 2); % each node has 2 DOFs (x and y)
counter = 0;
for ny = 2:n_np_y-1
    for nx = 2:n_np_x
        index0 = index((ny-1)*n_np_x + nx);
        counter = counter + 1;
        ID(index0, :) = [2*counter-1, 2*counter]; 
    end
end

n_eq = 2 * counter; % total number of equations

% LM array
LM = zeros(n_el, 2 * n_en); 
for ee = 1:n_el
    for aa = 1:n_en
        node = IEN(ee, aa);
        LM(ee, 2*aa-1:2*aa) = ID(node, :);
    end
end

% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 36 * n_eq);
% K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);

% loop over elements to assemble stiffness matrix and load vector
for ee = 1:n_el
    x_ele = xx_coor(IEN(ee, :));
    y_ele = yy_coor(IEN(ee, :));
    % coordinates of nodal points in each element
    
    k_ele = zeros(8, 8); % element stiffness matrix
    % f_ele = zeros(8, 1); % element load vector

    for ll = 1: n_int
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        for aa = 1 : n_en   
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end
        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi; % determinant of Jacobi matrix

        % compute B matrix 
        B = zeros(3, 8);
        for aa = 1: n_en
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            Na_x1 = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            Na_x2 = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
            B(:, (2*aa-1):(2*aa)) = [Na_x1, 0; 0, Na_x2; Na_x2, Na_x1];
        end
        
        % compute k_ele and f_ele
        k_ele = k_ele + weight(ll) * detJ * (B' * D * B);
        % for aa = 1: n_en
        %     f_ele(2*aa-1) = f_ele(2*aa-1) + weight(ll) * detJ * f_x(x_ele(aa), y_ele(aa));
        %     f_ele(2*aa) = f_ele(2*aa) + weight(ll) * detJ * f_y(x_ele(aa), y_ele(aa));
        % end
    end
    % assemble global stiffness matrix and load vector
    for aa = 1:2*n_en
        for bb = 1:2*n_en
            P = LM(ee, aa);
            Q = LM(ee, bb);
            if P > 0 && Q > 0
                K(P, Q) = K(P, Q) + k_ele(aa, bb);
            end
        end
        % if LM(ee, aa) > 0
        %     F(LM(ee, aa)) = F(LM(ee, aa)) + f_ele(aa);
        % end
    end
end

% apply Neumann BCs on x = 1
boundaryNodes = find(xx_coor == 1); % right boundary (x = 1)
for i = 1:length(boundaryNodes)-1
    n1 = boundaryNodes(i);
    n2 = boundaryNodes(i+1);
    x1 = xx_coor(n1); y1 = yy_coor(n1);
    x2 = xx_coor(n2); y2 = yy_coor(n2);
    L = sqrt((x2-x1)^2 + (y2-y1)^2);

    for ll = 1:2 % 2-point Gauss quadrature for 1D
        xi = [-sqrt(1/3), sqrt(1/3)]; % Gauss points
        w = [1, 1]; % weights
        N1 = (1-xi(ll))/2; N2 = (1+xi(ll))/2;
        x = N1*x1 + N2*x2; y = N1*y1 + N2*y2;
        t = traction(x, y);
        
        % add contributions to global force vector
        if ID(n1, 1)
            F(ID(n1, 1)) = F(ID(n1, 1)) + w(ll) * L/2 * t(1);
        end
        if ID(n1, 2)
            F(ID(n1, 2)) = F(ID(n1, 2)) + w(ll) * L/2 * t(2);
        end
        if ID(n2, 1)
            F(ID(n2, 1)) = F(ID(n2, 1)) + w(ll) * L/2 * t(1);
        end
        if ID(n2, 2)
            F(ID(n2, 2)) = F(ID(n2, 2)) + w(ll) * L/2 * t(2);
        end
    end
end

% solve the linear system
d = K \ F;

disp_xx = zeros(n_np, 1);
disp_yy = zeros(n_np, 1);
for ii = 1 : n_np
  index1 = ID(ii, 1); 
  index2 = ID(ii, 2);
  % ii is global node number and index is global equation number
  if index1 > 0
      disp_xx(ii) = d(index1);
  end
  if index2 > 0
      disp_yy(ii) = d(index2);
  end
end

plot3(xx_coor, yy_coor, disp_xx)
%plot3(xx_coor, yy_coor, disp_yy)