function val = alkylation_obj(x)
%ALKYLATION_OBJ Optimal operation of alkylation unit
%
%   Variable   Lower Bound   Upper Bound
%   ------------------------------------
%   x1                1500          2000
%   x2                   1           120
%   x3                3000          3500
%   x4                  85            93
%   x5                  90            95
%   x6                   3            12
%   x7                 145           162
%   ------------------------------------
%   The maximum profit is $1766.36 per day, and the optimal
%   variable values are x1 = 1698.256922, x2 = 54.274463,
%   x3 = 3031.357313, x4 = 90.190233, x5 = 95.0,
%   x6 = 10.504119, x7 = 153.535355.
%
%   Source:
%     Babu, B. V., and Angira, Rakesh (2006).
%     Modified differential evolution (MDE) for optimization of nonlinear
%     chemical processes.
%     Computers and Chemical Engineering 30, 989-1002.

x1 = x(1); x3 = x(3);
val = 1.715*x1 + 0.035*x1*x(6) + 4.0565*x3 +10.0*x(2) - 0.063*x3*x(5);

end

