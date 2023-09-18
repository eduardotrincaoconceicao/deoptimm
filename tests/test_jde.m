% Copyright 2019, Eduardo L. T. Conceicao
% Available under the GPL-3


rng('default');

% Unconstrained test problems ---------------------------------------------

%% Test 1: Schaffer 1 problem
out = jde(-[100 100], [100 100], @sf1, [], [], [], ...
          'NP', 50, 'tol', 1e-7, 'maxiter', 800);
assert(abs(out.value) <= 1e-4)

%% Test 2: Schwefel problem
out = jde(repelem(-500, 10), repelem(500, 10), @swf, [], [], [], ...
          'tol', 1e-7);
expVal = -418.9829*10;
assert(abs(expVal - out.value) <= 1e-4*abs(expVal));

% Only inequality constraints ---------------------------------------------

%% Test 3: Reactor network design
out = jde([1e-5 1e-5], [16 16], @RND_obj, @RND_con, [], [], ...
          'NP', 40, 'tol', 1e-7);
expVal = -0.388812;
assert(abs(expVal - out.value) <= 1e-2*abs(expVal));

%% Test 4: Heat exchanger network design
out = jde([100 1000 1000 10 10], [10000 10000 10000 1000 1000], ...
          @HEND_obj, @HEND_con, [], [], ...
          'tol', 1e-4);
expVal = 7049.25;
assert(abs(expVal - out.value) <= 1e-3*abs(expVal));

%% Test 5: Optimal operation of alkylation unit
out = jde([1500 1 3000 85 90 3 145], [2000 120 3500 93 95 12 162], ...
          @alkylation_obj, @alkylation_con, [], [], ...
          'tol', 0.1);
expVal = -1766.36;
assert(abs(expVal - out.value) <= 1e-2*abs(expVal));

% Only equality constraints -----------------------------------------------

%% Test 6: Test problem g11 from CEC2006 benchmarks
out = jde(-[1 1], [1 1], @g11_obj, @g11_con, 1, [], ...
          'tol', 1e-7);
expVal = 0.75;
assert(abs(expVal - out.value) <= 1e-2*abs(expVal));
