% Compute ratio alpha_{c,f}(mu1, mu2, rho) from the pairing algorithm.
%
% mu1, mu2, rho form a configuration
% 0 <= c <= 1 is the linear scaling factor
% f: [0,1] -> [0,1] is the "bias boosting" function
%     f must be able to take matrices as inputs.
%
% mu1, mu2, rho may be matrices, in which case they must all have the same dimension.
% The result in this case is a matrix of the same dimension giving the ratios on all
% the indicated configurations.
%
function ratio = alpha_cf(mu1, mu2, rho, c, f)
  mu1 = trunc(mu1, -1, 1);
  mu2 = trunc(mu2, -1, 1);
  rho = trunc(rho, -1, 1);

  % The "weak bounds"
  A1 = c*mu1;
  B1 = A1 + (1-c)*sign(mu1).*f(abs(mu1));
  A2 = c*mu2;
  B2 = A2 + (1-c)*sign(mu2).*f(abs(mu2));

%  [mu1, mu2, rho, A1, B1, A2, B2],
  % The boost that is given to one of the vertices in case of biases of opposite sign
  boost = (1-c)*f(min(abs(mu1), abs(mu2)));

  ratio = ones(size(mu1));
  for i = 1:size(mu1, 1)
    for j = 1:size(mu1, 2)
      I1 = [A1(i,j), B1(i,j)];
      I2 = [A2(i,j), B2(i,j)];
      if sign(mu1) == sign(mu2)
	ratio(i,j) = alpha_I(mu1(i,j), mu2(i,j), rho(i,j), I1, I2);
      else
	% Two cases depending on which vertex gets the boost
	J1 = [I1(1) + sign(mu1(i,j))*boost(i,j), I1(2)];
	J2 = [I2(1) + sign(mu2(i,j))*boost(i,j), I2(2)];
	ratio(i,j) = min(alpha_I(mu1(i, j), mu2(i, j), rho(i, j), J1, I2), ...
			 alpha_I(mu1(i, j), mu2(i, j), rho(i, j), I1, J2));
      end;
    end;
    if mod(i, 10) == 0
      disp(sprintf('%.3f%%', 100*i/size(mu1,1)));
    end;
  end;

function X = trunc(X, lo, hi)
   X(find(X < lo)) = lo;
   X(find(X > hi)) = hi;
