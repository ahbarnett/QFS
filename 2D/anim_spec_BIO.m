% script to animate spectrum of QFS BIO vs N.
% Barnett 8/29/19

tol = 1e-10;
curvemeth = 'n';   % 'n' slightly worse spec than 'i'
interior = true;
figure;
for N=50:20:600
  clf;
  lam = spec_BIO(tol,[],curvemeth,N,interior);
  title(sprintf('spec(A_{QFS}) on-surf BIO: tol=%5.3g, N=%d',tol,N))
  drawnow;
end
% we see a little jiggling but v stable & well-cond, for curvemeht='i'.
% 1e-10: curvemeth='n' shows cond(A) grow weakly to <1e2 once N > N_tol.
%  It then stabilizes back down to 1e1.

