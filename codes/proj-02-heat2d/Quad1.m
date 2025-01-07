function val = Quad1(aa, xi, eta) %% for triangular elements with three nodal points

if aa == 1
    val = eta;
elseif aa == 2
    val = 1-xi-eta;
elseif aa == 3
    val = xi;
else
    error('Error: value of a should be 1,2,or 3.');
end

% EOF