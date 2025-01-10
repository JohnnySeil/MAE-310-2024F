function A = IEN_quadrilateral(a, e, triangular)
if a == 1
    A = triangular(e * 2 - 1, 1);
elseif a == 2
    A = triangular(e * 2 - 1, 3);
elseif a == 3
    A = triangular(e * 2, 3);
elseif a == 4
    A = triangular(e * 2 - 1, 2);
end
end