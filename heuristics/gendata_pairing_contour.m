% Computes the ratio alpha_{c,f}(mu1, mu2, rho) from the pairing algorithm along the
% boundary of the polytope
%
% Parameters:
% 0 <= c <= 1 is the linear scaling factor
% f: [0,1] -> [0,1] is the "bias boosting" function
%     f must be able to take matrices as inputs.
%  n integer giving the granularity (step size will be 1/n so we will have 2*n+1 intervals)
%
function gendata_pairing_contour(c, f, n)
  N = 2*n+1;
  mu = -1:1/n:1;

  mu1 = mu'*ones(1, N);
  mu2 = ones(N, 1)*mu;

  rhomin = -1 + abs(mu1+mu2);
  rhomax = 1 - abs(mu1-mu2);

  % disp('Computing alpha along rho_min');
  % alpha1 = alpha_cf(mu1, mu2, rhomin, c, f);

  % fid = fopen('pairing-contour-rhomin.txt', 'w');
  % for i = 1:N
  %   fprintf(fid, '%.14f %.14f %.14f\n', [mu1(i,:); mu2(i,:); alpha1(i,:)]);
  %   fprintf(fid, '\n');
  % end;
  % fclose(fid);

  % disp('Computing alpha along rho_max');
  % alpha2 = alpha_cf(mu1, mu2, rhomax, c, f);

  % fid = fopen('pairing-contour-rhomax.txt', 'w');
  % for i = 1:N
  %   fprintf(fid, '%.14f %.14f %.14f\n', [mu1(i,:); mu2(i,:); alpha2(i,:)]);
  %   fprintf(fid, '\n');
  % end;
  % fclose(fid);

  rho = mu2;
  mu2min = max(rhomin, (-1+1e-5)*ones(size(rho)));
  mu2max = rhomax;

  disp('Computing alpha along mu2_min');
  alpha3 = alpha_cf(mu1, mu2min, rho, c, f);

  fid = fopen('pairing-contour-mu2min.txt', 'w');
  for i = 1:N
    fprintf(fid, '%.14f %.14f %.14f\n', [mu1(i,:); rho(i,:); alpha3(i,:)]);
    fprintf(fid, '\n');
  end;
  fclose(fid);

  % disp('Computing alpha along mu2_min');
  % alpha4 = alpha_cf(mu1, mu2max, rho, c, f);

  % fid = fopen('pairing-contour-mu2max.txt', 'w');
  % for i = 1:N
  %   fprintf(fid, '%.14f %.14f %.14f\n', [mu1(i,:); rho(i,:); alpha4(i,:)]);
  %   fprintf(fid, '\n');
  % end;
  % fclose(fid);

  % Note first level should not appear.
  levels = [0.877, 0.878, 0.88, 0.89, 0.9, 0.95, 1];

  figure;
  [c1, h1] = contour(mu1, rho, alpha3, levels);
  clabel(c1, h1);

%  figure;
%  [c2, h2] = contour(mu1, rho, alpha4, levels);
%  clabel(c2, h2);

