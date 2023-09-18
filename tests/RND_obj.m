function val = RND_obj(x)
%RND_OBJ Reactor network design
%   1e-5 <= x5, x6 <= 16
%   It possesses two local solutions at x = (16, 0) with f = -0.37461
%   and at x = (0, 16) with f = -0.38808.
%   The global optimum is (x5, x6; f) = (3.036504, 5.096052; -0.388812).
%
%   Source:
%     Babu, B. V., and Angira, Rakesh (2006).
%     Modified differential evolution (MDE) for optimization of nonlinear
%     chemical processes.
%     Computers and Chemical Engineering 30, 989-1002.

x5 = x(1); x6 = x(2);
k1 = 0.09755988; k2 = 0.99*k1; k3 = 0.0391908; k4 = 0.9*k3;
val = -( k2*x6*(1 + k3*x5) + k1*x5*(1 + k2*x6) )/( (1 + k1*x5)*(1 + k2*x6)* ...
                                                   (1 + k3*x5)*(1 + k4*x6) );

end

