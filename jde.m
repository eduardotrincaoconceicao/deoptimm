function res = jde(lower, upper, fn, options)
%JDE Nonlinear constrained minimization via Differential Evolution.
%   This implementation of the jDE algorithm [1] uses the
%   DE/rand/1/either-or mutation strategy [2]. It also immediately replaces
%   each worse parent vector, thereby utilizing only one population [3]
%   instead of two as in the original DE algorithm [4].
%
%   Examples:
%{
      jde(-[100 100], [100 100], @sf1, ...
          NP = 50, tol = 1e-7, maxiter = 800, trace = true)
      jde(repmat(-500, 10, 1), repmat(500, 10, 1), @swf, ...
          tol = 1e-7, trace = true)

      jde([1e-5 1e-5], [16 16], @RND_obj, constr = @RND_con, ...
          NP = 40, tol = 1e-7, trace = true)
      jde([100 1000 1000 10 10], [10000 10000 10000 1000 1000], ...
          @HEND_obj, constr = @HEND_con, ...
          tol = 1e-7, trace = true)
      jde([1500 1 3000 85 90 3 145], [2000 120 3500 93 95 12 162], ...
          @alkylation_obj, constr = @alkylation_con, ...
          tol = 1e-7, trace = true)
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

%   Copyright 2014, 2019, 2023, Eduardo L. T. Conceicao
%   Available under the GPL-3


    arguments
        lower (:, 1) {mustBeNonempty, mustBeNumeric, mustBeFinite}
        upper (:, 1) {mustBeNonempty, mustBeNumeric, mustBeFinite}
        fn function_handle {mustBeNonempty}
        options.constr = []
        options.meq (1, 1) ...
            {mustBeNonempty, mustBeInteger, mustBeNonnegative} = 0
        options.eps (:, 1) ...
            {mustBeNonempty, mustBePositive, mustBeFinite} = 1e-5
        options.NP (1, 1) ...
            {mustBeNonempty, mustBeInteger, mustBeNonnegative} ...
            = 10*length(lower)
        options.Fl (1, 1) ...
            {mustBeNonempty, mustBeNonnegative, mustBeFinite} = 0.1
        options.Fu (1, 1) ...
            {mustBeNonempty, mustBeNonnegative, mustBeFinite} = 1
        options.tau_F (1, 1) ...
            {mustBeNonempty, mustBeInRange(options.tau_F, 0, 1)} = 0.1
        options.tau_CR (1, 1) ...
            {mustBeNonempty, mustBeInRange(options.tau_CR, 0, 1)} = 0.1
        options.tau_pF (1, 1) ...
            {mustBeNonempty, mustBeInRange(options.tau_pF, 0, 1)} = 0.1
        options.jitter_factor (1, 1) ...
            {mustBeNonempty, mustBeNonnegative, mustBeFinite} = 0.001
        options.tol (1, 1) ...
            {mustBeNonempty, mustBeNumeric, mustBeFinite} = 1e-15
        options.maxiter (1, 1) ...
            {mustBeNonempty, mustBeInteger, mustBeNonnegative} ...
            = 200*length(lower)
        options.fnscale (1, 1) ...
            {mustBeNonempty, mustBePositive, mustBeFinite} = 1
        options.compare_to ...
            {mustBeTextScalar, ...
             mustBeMember(options.compare_to, ["median", "max"])} = "median"
        options.add_to_init_pop = []
        options.trace (1, 1) ...
            {mustBeA(options.trace, "logical")} = false
        options.triter (1, 1) ...
            {mustBeNonempty, mustBeInteger, ...
             mustBeGreaterThanOrEqual(options.triter, 1)} = 1
        options.details (1, 1) ...
            {mustBeA(options.details, "logical")} = false
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
        ftrial = fn(trial);
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
                pop(:, i) = trial;   % infeasible, the one with smaller
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
            if TAVtrial <= TAVpop(i) % trial and target are both infeasible
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
        h = h(:);
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

% Extract arguments
constr = options.constr;
meq = options.meq;
eps = options.eps;
NP = options.NP;
Fl = options.Fl;
Fu = options.Fu;
tau_F = options.tau_F;
tau_CR = options.tau_CR;
tau_pF = options.tau_pF;
jitter_factor = options.jitter_factor;
tol = options.tol;
maxiter = options.maxiter;
fnscale = options.fnscale;
compare_to = options.compare_to;
add_to_init_pop = options.add_to_init_pop;
trace = options.trace;
triter = options.triter;
details = options.details;

% Check input parameters
d = length(lower);
if length(upper) ~= d
    error('lower must have same length as upper')
end
if any(lower > upper)
    error('lower <= upper is not fulfilled')
end
if ~isempty(constr)
    validateattributes(constr, ...
                       {'function_handle'}, {'nonempty'}, ...
                       mfilename, 'constr')
    if ~(isscalar(eps) || length(eps) == meq)
        error('eps must be either of length meq, or length 1')
    end
end
if Fl > Fu
    error('Fl <= Fu is not fulfilled')
end
if ~isempty(add_to_init_pop)
    validateattributes(add_to_init_pop, ...
                       {'numeric'}, {'2d', 'nrows', d, 'finite', 'nonnan'}, ...
                       mfilename, 'add_to_init_pop')
    if any( (add_to_init_pop < repmat(lower, 1, size(add_to_init_pop, 2))) | ...
            (add_to_init_pop > repmat(upper, 1, size(add_to_init_pop, 2))) )
        error('add_to_init_pop must be between lower and upper')
    end
end

% Initialization
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
validateattributes(fpop, ...
                   {'numeric'}, {'row', 'ncols', NP, 'nonnan'}, ...
                   mfilename, 'fpop')
if ~isempty(constr)
    hpop = zeros(length( constr1(pop(:, 1)) ), NP);
    TAVpop = zeros(1, NP);
    for k = popIndex
        hpop(:, k) = constr1(pop(:, k));
        TAVpop(k) = sum( max(hpop(:, k), 0) );
    end
    validateattributes(hpop, ...
                       {'numeric'}, {'2d', 'ncols', NP, 'nonnan'}, ...
                       mfilename, 'hpop')
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

while rule() % Generation loop
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