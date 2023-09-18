function val = g11_obj(x)
%G11_OBJ CEC2006 benchmarks
%   -1 <= xi <= 1 (i = 1, 2)
%   The optimal solution is x* = (+-1/sqrt(2), 1/2)
%   and the optimal value f(x*) = 0.75.
%
%   Source:
%     Runarsson, Thomas P., and Yao, Xin (2000).
%     Stochastic ranking for constrained evolutionary optimization.
%     IEEE Transactions on Evolutionary Computing 4, 284-294.

val = x(1)^2 + (x(2) - 1)^2;

end

