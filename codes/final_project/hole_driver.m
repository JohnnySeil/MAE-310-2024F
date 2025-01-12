% elasticity parameters
E = 1e5;       % Young's modulus
vv = 0.3;      % Poisson's ratio

% Plane stress material stiffness matrix
D = E / (1 - vv^2) * [1, vv, 0; vv, 1, 0; 0, 0, (1 - vv) / 2];

% Load mesh generated from Gmsh
load('mesh.mat'); % Replace 'mesh.mat' with your actual mesh file name
nodes = msh.POS(:, 1:2);        % Nodal coordinates
triangular = msh.TRIANGLES(:, 1:3); % Element connectivity
tol = 1e-6; % Tolerance for boundary condition checks

% Generate nodal coordinates
n_np = size(nodes, 1); % Total number of nodes
n_el = size(triangular, 1); % Total number of elements
xx_coor = nodes(:, 1);
yy_coor = nodes(:, 2);

% Initialize arrays
n_en = 4; % Number of nodes per element (quadrilateral elements)
disp_xx = zeros(n_np, 1); % Displacement in x
disp_yy = zeros(n_np, 1); % Displacement in y
F = zeros(2 * n_np, 1); % Global force vector

% ID Array: Map global DOFs
ID = (1:2 * n_np)';
for i = 1:n_np
    if abs(xx_coor(i)) < tol % Symmetry boundary at x = 0
        ID(2 * i - 1) = 0; % Fix u_x
    end
    if abs(yy_coor(i)) < tol % Symmetry boundary at y = 0
        ID(2 * i) = 0; % Fix u_y
    end
end

% Apply traction boundary condition (σ_xx = 10 kPa) on x = 4
traction_force = 10; % Traction in kPa
for i = 1:n_np
    if abs(xx_coor(i) - 4) < tol % Right boundary
        F(2 * i - 1) = F(2 * i - 1) + traction_force; % Traction in x-direction
    end
end

% Stiffness matrix assembly
K = sparse(2 * n_np, 2 * n_np); % Sparse global stiffness matrix
for ee = 1:n_el
    % Element stiffness matrix assembly
    element_nodes = triangular(ee, :);
    x_ele = xx_coor(element_nodes);
    y_ele = yy_coor(element_nodes);
    [k_ele, f_ele] = computeElementStiffness(x_ele, y_ele, D); % Element stiffness and force
    for aa = 1:n_en
        for bb = 1:n_en
            P = 2 * element_nodes(aa) - 1; % x DOF
            Q = 2 * element_nodes(bb) - 1; % x DOF
            K(P, Q) = K(P, Q) + k_ele(2 * aa - 1, 2 * bb - 1); % xx term
            K(P + 1, Q + 1) = K(P + 1, Q + 1) + k_ele(2 * aa, 2 * bb); % yy term
        end
    end
end

% Apply boundary conditions by removing constrained DOFs
free_dofs = ID > 0; % Identify free DOFs
K_reduced = K(free_dofs, free_dofs);
F_reduced = F(free_dofs);

% Solve the reduced system
d_reduced = K_reduced \ F_reduced;

% Map results back to full displacement vector
disp = zeros(2 * n_np, 1);
disp(free_dofs) = d_reduced;

% Separate displacements
disp_xx = disp(1:2:end); % Displacement in x
disp_yy = disp(2:2:end); % Displacement in y

% Compute stresses
stress_xx = zeros(n_el, 1);
stress_yy = zeros(n_el, 1);
stress_xy = zeros(n_el, 1);
for ee = 1:n_el
    element_nodes = triangular(ee, :);
    x_ele = xx_coor(element_nodes);
    y_ele = yy_coor(element_nodes);
    ux_ele = disp_xx(element_nodes);
    uy_ele = disp_yy(element_nodes);
    [strain, stress] = computeStress(x_ele, y_ele, ux_ele, uy_ele, D); % Compute stress
    stress_xx(ee) = stress(1);
    stress_yy(ee) = stress(2);
    stress_xy(ee) = stress(3);
end

% Visualization
figure;
trisurf(triangular, xx_coor, yy_coor, disp_xx, 'EdgeColor', 'none');
colorbar;
title('Displacement in X');
xlabel('X [m]');
ylabel('Y [m]');

figure;
trisurf(triangular, xx_coor, yy_coor, stress_xx, 'EdgeColor', 'none');
colorbar;
title('\sigma_{xx} Stress Distribution');
xlabel('X [m]');
ylabel('Y [m]');

function [k_ele, f_ele] = computeElementStiffness(x_ele, y_ele, D)
    % Quadrature points and weights
    [xi, eta, w] = Gauss2D(3, 3);
    n_int = length(w);
    n_en = length(x_ele); % Number of nodes per element
    k_ele = zeros(2 * n_en, 2 * n_en);
    f_ele = zeros(2 * n_en, 1);

    for i = 1:n_int
        [N, dN_dxi, dN_deta] = shapeFunctions(xi(i), eta(i));
        [J, detJ, dN_dx, dN_dy] = computeJacobian(x_ele, y_ele, dN_dxi, dN_deta);
        B = strainDisplacementMatrix(dN_dx, dN_dy);
        k_ele = k_ele + B' * D * B * detJ * w(i);
    end
end

function [strain, stress] = computeStress(x_ele, y_ele, ux_ele, uy_ele, D)
    % Compute strain and stress at the element center
    [N, dN_dxi, dN_deta] = shapeFunctions(0, 0);
    [J, detJ, dN_dx, dN_dy] = computeJacobian(x_ele, y_ele, dN_dxi, dN_deta);
    B = strainDisplacementMatrix(dN_dx, dN_dy);
    u_ele = [ux_ele; uy_ele];
    strain = B * u_ele;
    stress = D * strain;
end