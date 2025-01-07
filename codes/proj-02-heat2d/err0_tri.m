n=[30,40,50,60,70,80,90,100];
for i=1:8
    [err0(i),err1(i),h(i)]=error1(n(i));
    err0(i)=log(err0(i));
    err1(i)=log(err1(i));
    h(i)=log(h(i));
end
polyfit(h,err0,1)
plot(h,err0,'r-o')
hold on
polyfit(h,err1,1)
plot(h,err1,'b-o')
grid on


function [errrr0, errrr1, hx]=error1(nn)

kappa = 1.0; % conductivity

% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% mesh generation
n_en   = 4;               % number of nodes in an element
n_el_x = nn;               % number of elements in x-dir
n_el_y = nn;               % number of elements in y-dir
n_el   = n_el_x * n_el_y; % total number of elements

n_np_x = n_el_x + 1;      % number of nodal points in x-dir
n_np_y = n_el_y + 1;      % number of nodal points in y-dir
n_np   = n_np_x * n_np_y; % total number of nodal points

x_coor = zeros(n_np, 1);
y_coor = x_coor;

hx = 1.0 / n_el_x;        % mesh size in x-dir
hy = 1.0 / n_el_y;        % mesh size in y-dir

% generate the nodal coordinates
for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx; % nodal index
    x_coor(index) = (nx-1) * hx;
    y_coor(index) = (ny-1) * hy;
  end
end

% IEN array  
IEN = zeros(n_el*2, n_en-1);
for ex = 1 : n_el_x
  for ey = 1 : n_el_y
    ee = 2*((ey-1) * n_el_x + ex)-1; % element index
    IEN(ee  , 1) = (ey-1) * n_np_x + ex;
    IEN(ee  , 2) =  ey    * n_np_x + ex;
    IEN(ee  , 3) =  ey    * n_np_x + ex + 1;
    IEN(ee+1, 1) =  ey    * n_np_x + ex + 1;
    IEN(ee+1, 2) = (ey-1) * n_np_x + ex + 1; 
    IEN(ee+1, 3) = (ey-1) * n_np_x + ex;
  end
end
% reconstruct the IEN array

% ID array
ID = zeros(n_np,1);
counter = 0;
for ny = 2 : n_np_y - 1
  for nx = 2 : n_np_x - 1
    index = (ny-1)*n_np_x + nx;
    counter = counter + 1;
    ID(index) = counter;  
  end
end
% ID array remains the same

n_eq = counter;

LM = ID(IEN);

% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 9 * n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
for ee = 1 : n_el*2
  x_ele = x_coor( IEN(ee, 1:n_en-1) );
  y_ele = y_coor( IEN(ee, 1:n_en-1) );
  
  k_ele = zeros(n_en-1, n_en-1); % element stiffness matrix
  f_ele = zeros(n_en-1, 1);    % element load vector
  
  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en-1
      x_l = x_l + x_ele(aa) * Quad1(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad1(aa, xi(ll), eta(ll));    
      [Na_xi, Na_eta] = Quad_grad1(aa, xi(ll), eta(ll));
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
    end
    
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
    
    for aa = 1 : n_en-1
      Na = Quad1(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Quad_grad1(aa, xi(ll), eta(ll));
      Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
      Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
      
      f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;
      
      for bb = 1 : n_en-1
        Nb = Quad1(bb, xi(ll), eta(ll));
        [Nb_xi, Nb_eta] = Quad_grad1(bb, xi(ll), eta(ll));
        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
        
        k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
      end % end of bb loop
    end % end of aa loop
  end % end of quadrature loop
 
  for aa = 1 : n_en-1
    PP = LM(ee, aa);
    if PP > 0
      F(PP) = F(PP) + f_ele(aa);
      
      for bb = 1 : n_en-1
        QQ = LM(ee, bb);
        if QQ > 0
          K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
        else
          % modify F with the boundary data
          % here we do nothing because the boundary data g is zero or
          % homogeneous
        end
      end  
    end
  end
end

% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, 1);

for ii = 1 : n_np
  index = ID(ii);
  if index > 0
    disp(ii) = dn(index);
  else
    % modify disp with the g data. Here it does nothing because g is zero
  end
end

%plot3(x_coor,y_coor,disp)
% save the solution vector and number of elements to disp with name
% HEAT.mat
% save("HEAT", "disp", "n_el_x", "n_el_y");

% EOF
npq=10; % select quadrature points in one element
[xii, etaa, w] = Gauss2D(npq, npq);
errrr0=0; errrr1=0;

for ee = 1:n_el*2
    xx_ele = x_coor(IEN(ee, :));
    yy_ele = y_coor(IEN(ee, :));
    u_ele = disp(IEN(ee, :));
    
    for ll = 1 : npq*npq % totsl number of quadrature points per element
        xx_l = 0; yy_l = 0; uh = 0; duh_dxi = 0; duh_deta = 0;
        dx_dxi = 0; dx_deta = 0; dy_dxi = 0; dy_deta = 0;
        
        for aa = 1:n_en-1
            xx_l = xx_l + xx_ele(aa)*Quad1(aa,xii(ll),etaa(ll));
            yy_l = yy_l + yy_ele(aa)*Quad1(aa,xii(ll),etaa(ll));
            uh = uh + u_ele(aa)*Quad1(aa,xii(ll),etaa(ll));
            [Na_xii, Na_etaa] = Quad_grad1(aa, xii(ll), etaa(ll));
            dx_dxi = dx_dxi + xx_ele(aa)*Na_xii;
            dy_dxi = dy_dxi + yy_ele(aa)*Na_xii;
            dx_deta = dx_deta + xx_ele(aa)*Na_etaa;
            dy_deta = dy_deta + yy_ele(aa)*Na_etaa;
            duh_dxi = duh_dxi + u_ele(aa)*Na_xii;
            duh_deta = duh_deta + u_ele(aa)*Na_etaa;
        end
        detj = dx_dxi * dy_deta - dx_deta * dy_dxi; % determinant of Jacobi matrix
        duh_dx = 1/detj * (duh_dxi*dy_deta - duh_deta*dy_dxi);
        duh_dy = 1/detj * (-duh_dxi*dx_deta + duh_deta*dx_dxi);

        errrr0 = errrr0 + w(ll) * detj * (uh - exact(xx_l, yy_l))^2;
        errrr1 = errrr1 + w(ll) * detj * ( (uh - exact(xx_l, yy_l))^2 + (duh_dx - exact_x(xx_l, yy_l))^2 + (duh_dy - exact_y(xx_l, yy_l))^2 );
    end
end

errrr0 = errrr0^0.5;
errrr1 = errrr1^0.5;



end
