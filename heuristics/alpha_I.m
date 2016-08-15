% Compute ratio alpha(mu1, mu2, rho, I1, I2) from the pairing algorithm.
%
% mu1, mu2, rho form a configuration
% I1, I2 are intervals (in the form of 2x1 or 1x2-dimensional vectors)
%
% Note: this function does NOT accept matrix arguments.
%
function ratio = alpha_I(mu1, mu2, rho, I1, I2)
  tr = trho(mu1, mu2, rho);
  
  % Find maximum lambda (i.e., minimum ratio)
  % on the four combinations of endpoints
  L = max(Lambda([I1(1) I1(1) I1(2) I1(2)], ...
		 [I2(1) I2(2) I2(1) I2(2)], ...
		 [tr tr tr tr]));

  if tr > 1e-15
    % trho > 0: find minimum ratio on the five other possibilities
    % For these we must be a bit more careful and check that they are in range.
    if inival(0, I1) && inival(0, I2)
      L = max(L, Lambda(0, 0, tr));
    end;

    g = 1 - 2*Phi(normalinv((1-I1(1))/2)/tr);
    if inival(g, I2)
      L = max(L, Lambda(I1(1), g, tr));
    end;
    
    g = 1 - 2*Phi(normalinv((1-I1(2))/2)/tr);
    if inival(g, I2)
      L = max(L, Lambda(I1(2), g, tr));
    end;
    
    g = 1 - 2*Phi(normalinv((1-I2(1))/2)/tr);
    if inival(g, I1)
      L = max(L, Lambda(g, I2(1), tr));
    end;
    
    g = 1 - 2*Phi(normalinv((1-I2(2))/2)/tr);
    if inival(g, I1)
      L = max(L, Lambda(g, I2(2), tr));
    end;
  end;

  if rho >= 1
    ratio = 1;
  else
    ratio = 2*(1-L) ./ (1-rho);
    % Truncate excessive values (can happen if rho is close to 1)
    if ratio >= 100000
      ratio = 100000;
    end;
  end;
 
%  disp(sprintf('%.4f %.4f %.4f (trho %.4f) [%.4f, %.4f] [%.4f, %.4f] -> L %.5f, rat %.5f', mu1, mu2, rho, tr, I1, I2, L, ratio));


% Utility function to check if x is in the interval I.
function ret = inival(x, I)
  ret = (x-I(1))*(x-I(2)) < 0;
