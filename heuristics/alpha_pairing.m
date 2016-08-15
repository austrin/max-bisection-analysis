% Compute ratio alpha(c, f) from the pairing algorithm.
%
% Searches for worst configuration by using Matlab's Optimization
% Toolbox, decrease risk of error by running the search many times
% with different starting points.
%
% Three optional arguments:
% - trials (default 50): number of searches to run.
% - eps (default 1e-13): tolerance (TolFun and TolX) option for fmincon
% - display (default 'none'): display option for fmincon ('iter' sometimes useful)
function a = alpha_pairing(c, f, trials, eps, display);

if nargin < 3
  trials = 50;
end;
if nargin < 4
  eps = 1e-13;
end;
if nargin < 5
  display = 'none';
end;

options = optimset('Display', display, ...
		   'Diagnostics', 'off', ...
		   'Algorithm', 'active-set', ...
		   'LargeScale', 'off', ...
		   'FunValCheck', 'on', ...
		   'MaxIter', 1000, ...
		   'TolFun', eps, ...
		   'TolX', eps);

% Triangle inequality constraints
A =  [-1    -1    -1
      -1     1    1
       1    -1    1
       1     1    -1];
B = [ 1-eps
      1-eps
      1-eps
      1-eps ];
% Upper and lower bounds on the variables
hi = [1-eps 1-eps 1-eps];
lo = -hi;

a = 2;

% Some known trouble-some configurations that we want to make sure to
% use as starting points.
tweak_configurations = [0.9999, -0.9999,  -0.9998;
			0.176945, 0.176945, -0.646110;
			0.55  -0.5 -0.687;
			0.995 -0.70  -0.70;]';

% Remember what the worst case configuration was for the last search:
% we will use it as first starting point this time
global alpha_pairing_last_worst; 
if length(alpha_pairing_last_worst) == 3
    tweak_configurations = [alpha_pairing_last_worst, tweak_configurations];
end;

for i = 1:trials
    if i <= size(tweak_configurations, 2)
      % First go through known trouble-some configurations
      v_start = tweak_configurations(:, i);
    else
      % After that just go with random starting points
      v_start = random_configuration;
    end;
    [v, aa] = fmincon(@alpha_pairing_wrapper, v_start, A, B, [], [], lo, hi, [], options, c, f);
    if (aa < a)
        a = aa;
        alpha_pairing_last_worst = v;
    end;
end;


function a=alpha_pairing_wrapper(v, c, f)
  a = alpha_cf(v(1), v(2), v(3), c, f);
%  disp(sprintf('%.6f ', [v' a]));
