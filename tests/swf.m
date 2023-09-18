function val = swf(x)
%SWF Schwefel problem
%   -500 <= xi <= 500, i = {1, 2, ..., n}
%   The number of local minima for a given n is not known, but the global
%   minimum f(x*) = -418.9829n is located at x* = (s, s, ..., s),
%   s = 420.97.
%
%   Source:
%     Ali, M. Montaz, Khompatraporn, Charoenchai, and
%     Zabinsky, Zelda B. (2005).
%     A numerical evaluation of several stochastic algorithms on selected
%     continuous global optimization test problems.
%     Journal of Global Optimization 31, 635-672.

% val = -dot( x, sin(sqrt(abs(x))) );
val = -sum( x .* sin(sqrt(abs(x))) );

end

