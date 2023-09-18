function val = RND_con(x)
%RND_CON Reactor network design
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

val = sqrt(x(1)) + sqrt(x(2)) - 4;

end

