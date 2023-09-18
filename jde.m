function res = jde(lower, upper, fn, constr, meq, eps, varargin)
%JDE Nonlinear constrained minimization via Differential Evolution.
%   This implementation of the jDE algorithm [1] uses the
%   DE/rand/1/either-or mutation strategy [2]. It also immediately replaces
%   each worse parent vector, thereby utilizing only one population [3]
%   instead of two as in the original DE algorithm [4].
%
%   Examples:
%{
      jde(-[100 100], [100 100], @sf1, [], [], [], ...
          'NP', 50, 'tol', 1e-7, 'maxiter', 800, 'trace', true)
      jde(repmat(-500, 10, 1), repmat(500, 10, 1), @swf, [], [], [], ...
          'tol', 1e-7, 'trace', true)

      jde([1e-5 1e-5], [16 16], @RND_obj, @RND_con, [], [], ...
          'NP', 40, 'tol', 1e-7, 'trace', true)
      jde([100 1000 1000 10 10], [10000 10000 10000 1000 1000], ...
          @HEND_obj, @HEND_con, [], [], ...
          'tol', 1e-7, 'trace', true)
      jde([1500 1 3000 85 90 3 145], [2000 120 3500 93 95 12 162], ...
          @alkylation_obj, @alkylation_con, [], [], ...
          'tol', 1e-7, 'trace', true)
%}

%   References:
%     [1] Brest, Janez, Greiner, Saso, Boskovic, Borko, Mernik, Marjan, and
%         Zumer, Viljem (2006).
%         Self-adapting control parameters in differential evolution: a
%         comparative study on numerical benchmark problems.
%         IEEE Transactions on Evolutionary Computation 10, 646-657.
%
%     [2] Price, Kenneth V., Storn, Rainer M., and
%         Lampinen, Jouni A. (2005).
%         Differential Evolution: a practical approach to global
%         optimization.
%         Springer, Berlin, pp. 117-118.
%
%     [3] Babu, B. V., and Angira, Rakesh (2006).
%         Modified differential evolution (MDE) for optimization of
%         non-linear chemical processes.
%         Computers and Chemical Engineering 30, 989-1002.
%
%     [4] Storn, Rainer, and Price, Kenneth (1997).
%         Differential evolution - a simple and efficient heuristic for
%         global optimization over continuous spaces.
%         Journal of Global Optimization 11, 341-359.

%   Copyright 2014, 2019, Eduardo L. T. Conceicao
%   Available under the GPL-3


    function parseInVarargin(argsToParse)
    % Parse input arguments
        for ii = 1:2:length(argsToParse)
            validateattributes(argsToParse{ii}, {'char'}, {'nonempty'}, ...
                               mfilename)
            validatestring(argsToParse{ii}, ...
                          {'NP', 'Fl', 'Fu', ...
                           'tau_F', 'tau_CR', 'tau_pF', ...
                           'jitter_factor', ...
                           'tol', 'maxiter', 'fnscale', ...
                           'compare_to', ...
                           'add_to_init_pop', 'trace', 'triter', ...
                           'details'}, ...
                           mfilename, argsToParse{ii});
            switch argsToParse{ii}
                case 'NP', NP = argsToParse{ii+1};
                case 'Fl', Fl = argsToParse{ii+1};
                case 'Fu', Fu = argsToParse{ii+1};
                case 'tau_F', tau_F = argsToParse{ii+1};
                case 'tau_CR', tau_CR = argsToParse{ii+1};
                case 'tau_pF', tau_pF = argsToParse{ii+1};
                case 'jitter_factor', jitter_factor = argsToParse{ii+1};
                case 'tol', tol = argsToParse{ii+1};
                case 'maxiter', maxiter = argsToParse{ii+1};
                case 'fnscale', fnscale = argsToParse{ii+1};
                case 'compare_to', compare_to = argsToParse{ii+1};
                case 'add_to_init_pop', add_to_init_pop = argsToParse{ii+1};
                case 'trace', trace = argsToParse{ii+1};
                case 'triter', triter = argsToParse{ii+1};
                case 'details', details = argsToParse{ii+1};
            end
        end
    end

    function x = handle_bounds(x, u)
    % Check feasibility of bounds and enforce parameters limits
    % by a deterministic variant of bounce-back resetting
    % (also known as midpoint base)
    % Price, KV, Storn, RM, and Lampinen, JA (2005)
    % Differential Evolution: A Practical Approach to Global Optimization.
    % Springer, p 206
        bad = x > upper;
        x(bad) = 0.5*(upper(bad) + u(bad));
        bad = x < lower;
        x(bad) = 0.5*(lower(bad) + u(bad));
    end

    % Mutate/recombine
    function trial = performReproduction()
        ignore = rand(d, 1) > CRtrial;
        if all(ignore)                      % ensure that trial gets at least
            ignore(randperm(d, 1)) = false; % one mutant parameter
        end
        % Source for trial is the base vector plus weighted differential
        if rand(1) <= pFtrial
            trial = X_base + Ftrial.*(X_r1 - X_r2);
        else
            trial = X_base + 0.5*(Ftrial + 1).*(X_r1 + X_r2 - 2*X_base);
        end
        % or trial parameter comes from target vector X_i itself.
        trial(ignore) = X_i(ignore);
    end

    % Evaluate/select ... -------------------------------------------------
    function child_()
        ftrial = fn(trial); % Evaluate trial with your function
        if ftrial <= fpop(i)
            pop(:, i) = trial;
            fpop(i) = ftrial;
            F(:, i) = Ftrial;
            CR(i) = CRtrial;
            pF(i) = pFtrial;
        end
    end

    % Equality constraints are present alongside the inequalities
    function child_constr_all()
    % Zhang, Haibo, and Rangaiah, G. P. (2012).
    % An efficient constraint handling method with integrated differential
    % evolution for numerical and engineering optimization.
    % Computers and Chemical Engineering 37, 74-88.
        htrial = constr1(trial);
        TAVtrial = sum( max(htrial, 0) );
        if TAVtrial > mu
            if TAVtrial <= TAVpop(i) % trial and target are both
                pop(:, i) = trial;   % unfeasible, the one with smaller
                hpop(:, i) = htrial; % constraint violation is chosen
                F(:, i) = Ftrial;    % or trial vector when both are
                CR(i) = CRtrial;     % solutions of equal quality
                pF(i) = pFtrial;
                TAVpop(i) = TAVtrial;
            end
        elseif TAVpop(i) > mu        % trial is feasible and target is not
            pop(:, i) = trial;
            fpop(i) = fn(trial);
            hpop(:, i) = htrial;
            F(:, i) = Ftrial;
            CR(i) = CRtrial;
            pF(i) = pFtrial;
            TAVpop(i) = TAVtrial;
        else                         % between two feasible solutions, the
            ftrial = fn(trial);      % one with better objective function
            if ftrial <= fpop(i)     % value is chosen
                pop(:, i) = trial;   % or trial vector when both are
                fpop(i) = ftrial;    % solutions of equal quality
                hpop(:, i) = htrial;
                F(:, i) = Ftrial;
                CR(i) = CRtrial;
                pF(i) = pFtrial;
                TAVpop(i) = TAVtrial;
                FF = sum(TAVpop <= mu)/NP;
                mu = mu*(1 - FF/NP);
            end
        end
    end

    % Only inequality constraints are present
    function child_constr_ineq()
        htrial = constr1(trial);
        TAVtrial = sum( max(htrial, 0) );
        if TAVtrial > mu
            if TAVtrial <= TAVpop(i) % trial and target are both unfeasible
                pop(:, i) = trial;
                hpop(:, i) = htrial;
                F(:, i) = Ftrial;
                CR(i) = CRtrial;
                pF(i) = pFtrial;
                TAVpop(i) = TAVtrial;
            end
        elseif TAVpop(i) > mu        % trial is feasible and target is not
            pop(:, i) = trial;
            fpop(i) = fn(trial);
            hpop(:, i) = htrial;
            F(:, i) = Ftrial;
            CR(i) = CRtrial;
            pF(i) = pFtrial;
            TAVpop(i) = TAVtrial;
            FF = sum(TAVpop <= mu)/NP;
            mu = mu*(1 - FF/NP);
        else                         % two feasible solutions
            ftrial = fn(trial);
            if ftrial <= fpop(i)
                pop(:, i) = trial;
                fpop(i) = ftrial;
                hpop(:, i) = htrial;
                F(:, i) = Ftrial;
                CR(i) = CRtrial;
                pF(i) = pFtrial;
                TAVpop(i) = TAVtrial;
                FF = sum(TAVpop <= mu)/NP;
                mu = mu*(1 - FF/NP);
            end
        end
    end
    % ... -----------------------------------------------------------------

    function ind = which_best_(x)
        [~, ind] = min(x);
    end

    function ind = which_best_constr(x)
        idx = TAVpop <= mu;
        if all(idx)
            ind = which_best_(x);
        elseif any(idx)
            ind = find(idx);
            ind = ind(which_best_(x(idx)));
        else
            ind = which_best_(TAVpop);
        end
    end

    function h = constr_eq(par)
        h = constr(par);
        h(equalIndex) = abs(h(equalIndex)) - eps;
    end

    function conv = stop_crit()
    % Zielinski, Karin, and Laur, Rainer (2008).
    % Stopping criteria for differential evolution in
    % constrained single-objective optimization.
    % In: U. K. Chakraborty (Ed.), Advances in Differential Evolution,
    % SCI 143, Springer-Verlag, pp 111-138
        conv = ( feval(compare_to, fpop) - fpop(x_best_ind) ) / fnscale;
    end

    function conv = rule_()
        conv = converge >= tol;
    end

    function conv = rule_constr()
        conv = (converge >= tol) || any(hpop(:, x_best_ind) > 0);
    end

narginchk(3, 36) % Check number of input arguments
% Check input parameters
validateattributes(lower, ...
                   {'numeric'}, {'vector', 'finite'}, ...
                   mfilename, 'lower', 1)
validateattributes(upper, ...
                   {'numeric'}, {'vector', 'finite'}, ...
                   mfilename, 'upper', 2)
d = length(lower);
if length(upper) ~= d
    error('lower must have same length as upper')
end
lower = lower(:);
upper = upper(:);
if any(lower > upper)
    error('lower <= upper is not fulfilled')
end
validateattributes(fn, {'function_handle'}, {'nonempty'}, mfilename, 'fn', 3)

% Deal with missing arguments
% Set default values
if nargin < 4 || isempty(constr), constr = []; end
if nargin < 5 || isempty(meq), meq = 0; end
if nargin < 6 || isempty(eps), eps = 1e-5; end
NP = 10*d;
Fl = 0.1;
Fu = 1;
tau_F = 0.1;
tau_CR = 0.1;
tau_pF = 0.1;
jitter_factor = 0.001;
tol = 1e-15;
maxiter = 200*d;
fnscale = 1;
compare_to = 'median';
add_to_init_pop = [];
trace = false;
triter = 1;
details = false;
parseInVarargin(varargin)

% Check input parameters
if ~isempty(constr)
    validateattributes(constr, ...
                       {'function_handle'}, {'nonempty'}, ...
                       mfilename, 'constr', 4)
    validateattributes(meq, ...
                       {'numeric'}, {'scalar', 'integer', 'nonnegative'}, ...
                       mfilename, 'meq', 5)
    validateattributes(eps, ...
                       {'numeric'}, {'positive', 'finite'}, ...
                       mfilename, 'eps', 6)
    if isscalar(eps), eps = repmat(eps, meq, 1); end
    if length(eps) ~= meq
        error('eps must be either of length meq, or length 1')
    end
end
validateattributes(NP, ...
                   {'numeric'}, {'scalar', 'integer'}, ...
                   mfilename, 'NP')
validateattributes(Fl, ...
                   {'numeric'}, {'scalar','nonnegative'}, ...
                   mfilename, 'Fl')
validateattributes(Fu, ...
                   {'numeric'}, {'scalar','nonnegative'}, ...
                   mfilename, 'Fu')
if Fl > Fu
    error('Fl <= Fu is not fulfilled')
end
validateattributes(tau_F, ...
                   {'numeric'}, {'scalar', '>=', 0, '<=', 1}, ...
                   mfilename, 'tau_F')
validateattributes(tau_CR, ...
                   {'numeric'}, {'scalar', '>=', 0, '<=', 1}, ...
                   mfilename, 'tau_CR')
validateattributes(tau_pF, ...
                   {'numeric'}, {'scalar', '>=', 0, '<=', 1}, ...
                   mfilename, 'tau_pF')
validateattributes(jitter_factor, ...
                   {'numeric'}, {'scalar','nonnegative'}, ...
                   mfilename, 'jitter factor')
validateattributes(tol, ...
                   {'numeric'}, {'scalar'}, ...
                   mfilename, 'tol')
validateattributes(maxiter, ...
                   {'numeric'}, {'scalar', 'integer','nonnegative'}, ...
                   mfilename, 'maxiter')
validateattributes(fnscale, ...
                   {'numeric'}, {'scalar', 'positive'}, ...
                   mfilename, 'fnscale')
validateattributes(compare_to, ...
                   {'char'}, {'nonempty'}, mfilename,'compare_to')
validatestring(compare_to, {'median', 'mean'}, mfilename,'compare_to');
if ~isempty(add_to_init_pop)
    validateattributes(add_to_init_pop, ...
                       {'numeric'}, {'2d', 'nrows', d, 'nonnan'}, ...
                       mfilename, 'add_to_init_pop')
    if any( (add_to_init_pop < repmat(lower, 1, size(add_to_init_pop, 2))) | ...
            (add_to_init_pop > repmat(upper, 1, size(add_to_init_pop, 2))) )
        error('add_to_init_pop must be between lower and upper')
    end
end
validateattributes(trace, {'logical'}, {'scalar'}, mfilename, 'trace')
validateattributes(triter, ...
                   {'numeric'}, {'scalar', 'integer', '>=', 1}, ...
                   mfilename, 'triter')
validateattributes(details, {'logical'}, {'scalar'}, mfilename, 'details')

% Initialization:
if ~isempty(constr)
    if meq > 0
        equalIndex = 1:meq;
        constr1 = @(par) constr_eq(par);
        child_constr = @() child_constr_all();
    else
        constr1 = @(par) constr(par);
        child_constr = @() child_constr_ineq();
    end
end

use_jitter = ~(jitter_factor == 0);

pop = unifrnd( repmat(lower, 1, NP), repmat(upper, 1, NP) );
if ~isempty(add_to_init_pop)
    pop = horzcat(pop, add_to_init_pop);
    NP = size(pop, 2);
end
if NP < 4
    error('NP must be at least 4')
end
% Combine jitter with dither
% Storn, Rainer (2008).
% Differential evolution research - trends and open questions.
% In: U. K. Chakraborty (Ed.), Advances in Differential Evolution,
% SCI 143, Springer-Verlag, pp 11-12
if use_jitter
    F = unifrnd(Fl, Fu, 1, NP) .* (1 + jitter_factor*unifrnd(-0.5, 0.5, d, 1));
else
    F = unifrnd(Fl, Fu, 1, NP);
end
CR = rand(NP, 1);
pF = rand(NP, 1);
popIndex = 1:NP;
fpop = zeros(1, NP);
for k = popIndex
    fpop(k) = fn(pop(:, k));
end
if ~isempty(constr)
    hpop = zeros(length( constr1(pop(:, 1)) ), NP);
    TAVpop = zeros(1, NP);
    for k = popIndex
        hpop(:, k) = constr1(pop(:, k));
        TAVpop(k) = sum( max(hpop(:, k), 0) );
    end
    mu = median(TAVpop);
end

if ~isempty(constr)
    child = @() child_constr();
    which_best = @(x) which_best_constr(x);
    rule = @() rule_constr();
else
    child = @() child_();
    which_best = @(x) which_best_(x);
    rule = @() rule_();
end
x_best_ind = which_best(fpop);
converge = stop_crit();
convergence = 0;
iteration = 0;

while rule() % generation loop
    if iteration >= maxiter
        warning('maximum number of iterations reached without convergence')
        convergence = 1;
        break
    end
    iteration = iteration + 1;

    for i = popIndex % Start loop through population

        % Equalize the mean lifetime of all vectors
        % Price, KV, Storn, RM, and Lampinen, JA (2005)
        % Differential Evolution: A Practical Approach to Global Optimization.
        % Springer, p 284
        i = mod(iteration + i, NP) + 1;

        % Self-adjusting parameter control scheme
        if rand(1) <= tau_F
            % Combine jitter with dither
            if use_jitter
                Ftrial = unifrnd(Fl, Fu) * ...
                         (1 + jitter_factor*unifrnd(-0.5, 0.5, d, 1));
            else
                Ftrial = unifrnd(Fl, Fu);
            end
        else
            Ftrial = F(:, i);
        end

        if rand(1) <= tau_CR
            CRtrial = rand(1);
        else
            CRtrial = CR(i);
        end

        if rand(1) <= tau_pF
            pFtrial = rand(1);
        else
            pFtrial = pF(i);
        end

        % DE/rand/1/either-or/bin
        X_i = pop(:, i);
        % Randomly pick 3 vectors all diferent from target vector
%         r = datasample( [1:(i-1) (i+1):NP], 3, 'Replace', false);
        r = [1:(i-1) (i+1):NP];
        r = r( randperm(NP-1, 3) );
        X_base = pop(:, r(1));
        X_r1 = pop(:, r(2));
        X_r2 = pop(:, r(3));

        trial = handle_bounds(performReproduction(), X_base);

        child()

        x_best_ind = which_best(fpop);
    end

    converge = stop_crit();
    if trace && (mod(iteration, triter) == 0)
        fprintf('%d: <%g> (%g)', ...
                iteration, converge, fpop(x_best_ind))
        fprintf(' %g', pop(:, x_best_ind))
        if ~isempty(constr)
            fprintf(' {')
            fprintf(' %d', find(hpop(:, x_best_ind) > 0))
            fprintf(' }')
        end
        fprintf('\n')
    end
end

res = struct('par', pop(:, x_best_ind), ...
             'value', fpop(x_best_ind), ...
             'iter', iteration, ...
             'convergence', convergence);
if details
    res.poppar = pop;
    res.popcost = fpop;
end
             
end