function [val_xi, val_eta] = Quad_grad1(aa, xi, eta) %% for triangular elements with three nodal points

if aa == 1
    val_xi  = 0;
    val_eta = 1;
elseif aa == 2
    val_xi  = -1;
    val_eta = -1;
elseif aa == 3
    val_xi  = 1;
    val_eta = 0;
else
    error('Error: value of a should be 1,2,or 3.');
end

% EOF