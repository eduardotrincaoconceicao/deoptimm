function val = HEND_obj(x)
%HEND_OBJ Heat exchanger network design
%   100 <= x1 <= 10000,   1000 <= x2, x3 <= 10000,
%   10 <= x4, x5 <= 1000
%   The global optimum is (x1, x2, x3, x4, x5: f) =
%   (579.19, 1360.13, 5109.92, 182.01, 295.60; 7049.25).
%
%   Source:
%     Babu, B. V., and Angira, Rakesh (2006).
%     Modified differential evolution (MDE) for optimization of nonlinear
%     chemical processes.
%     Computers and Chemical Engineering 30, 989-1002.

val = x(1) + x(2) + x(3);

end

