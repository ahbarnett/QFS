function fun=interpfunfromgrid(y,ncomp)  % gives a func handle (scalar or vec)
% this is for 1d grids, periodic only. Wrapper to perispecinterparb from BIE2D.
if ncomp==1
  [~,fun] = perispecinterparb(y,nan);   % scalar spectral interpolant
elseif ncomp==2
  N = numel(y)/ncomp;
  [~,fun1] = perispecinterparb(y(1:N),nan);   % vector spectral interpolant
  [~,fun2] = perispecinterparb(y(N+1:end),nan);   % vector spectral interpolant
  fun = @(t) [fun1(t); fun2(t)];                 % stack cmpts
else, error('interpfunfromgrid only for ncomp=1 or 2!');
end
