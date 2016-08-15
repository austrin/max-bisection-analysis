function gendata_linear_limit()
  phi1 = [0.176945, 0.176945, -0.646110];
  phi2 = [1, -1, -1];
  c = (0:0.001:1)';

  Phi1 = ones(size(c))*phi1;
  Phi2 = ones(size(c))*phi2;

  alpha1 = alpha(Phi1(:,1), Phi1(:,2), Phi1(:,3), c*phi1(1), c*phi1(2));
  alpha2 = alpha(Phi2(:,1), Phi2(:,2), Phi2(:,3), c*phi2(1), c*phi2(2));
  sdp1 = (1-phi1(3))/2;
  sdp2 = (1-phi2(3))/2;
  val1 = alpha1*sdp1;
  val2 = alpha2*sdp2;
  
  mixprob = 0.931935;

  fid = fopen('linear-bad-mixed-config.txt', 'w');

  alphamix = (mixprob*val1 + (1-mixprob)*val2)/(mixprob*sdp1 + (1-mixprob)*sdp2);

  fprintf(fid, '%.10f %.10f\n', [c'; alpha1']);
  fprintf(fid, '\n\n');
  fprintf(fid, '%.10f %.10f\n', [c'; alpha2']);
  fprintf(fid, '\n\n');
  fprintf(fid, '%.10f %.10f\n', [c'; alphamix']);

  max(alphamix),

  plot(c, alpha1, c, alpha2, c, alphamix);
  axis([0 1 0.873 0.874]);