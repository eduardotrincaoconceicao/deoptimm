function val = alkylation_con(x)
%ALKYLATION_CON Optimal operation of alkylation unit
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

x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);
x5 = x(5); x6 = x(6); x7 = x(7);
val = [0.0059553571*x6^2*x1 + 0.88392857*x3 - 0.1175625*x6*x1 - x1;
       1.1088*x1 + 0.1303533*x1*x6 - 0.0066033*x1*x6^2 - x3;
       6.66173269*x6^2 + 172.39878*x5 - 56.596669*x4 - 191.20592*x6 - 10000;
       1.08702*x6 + 0.32175*x4 - 0.03762*x6^2 - x5 + 56.85075;
       0.006198*x7*x4*x3 + 2462.3121*x2 - 25.125634*x2*x4 - x3*x4;
       161.18996*x3*x4 + 5000.0*x2*x4 - 489510.0*x2 - x3*x4*x7;
       0.33*x7 - x5 + 44.333333;
       0.022556*x5 - 0.007595*x7 - 1.0;
       0.00061*x3 - 0.0005*x1 - 1.0;
       0.819672*x1 - x3 + 0.819672;
       24500.0*x2 - 250.0*x2*x4 - x3*x4;
       1020.4082*x4*x2 - 1.2244898*x3*x4 - 100000*x2;
       6.25*x1*x6 + 6.25*x1 - 7.625*x3 - 100000;
       1.22*x3 - x6*x1 - x1 + 1.0];

end

