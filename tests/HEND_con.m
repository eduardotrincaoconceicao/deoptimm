function val = HEND_con(x)
%HEND_CON Heat exchanger network design
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

x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4); x5 = x(5);
val = [100*x1 - x1*(400 - x4) + 833.33252*x4 - 83333.333;
       x2*x4 - x2*(400 - x5 + x4) - 1250*x4 + 1250*x5;
       x3*x5 - x3*(100 + x5) - 2500*x5 + 1250000];

end

