function val = sf1(x)
%SF1 Schaffer 1 problem
%   -100 <= x1, x2 <= 100
%   The number of local minima is not known but the global minimum is
%   located at x* = (0, 0) with f(x*) = 0.
%
%   Source:
%     Ali, M. Montaz, Khompatraporn, Charoenchai, and
%     Zabinsky, Zelda B. (2005).
%     A numerical evaluation of several stochastic algorithms on selected
%     continuous global optimization test problems.
%     Journal of Global Optimization 31, 635-672.

temp = x(1)^2 + x(2)^2;
val = 0.5 + (sin(sqrt(temp))^2 - 0.5)/(1 + 0.001*temp)^2;

end

