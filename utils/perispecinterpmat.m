function I = perispecinterpmat(Nout,Nin)
% PERISPECINTERPMAT.  Return dense matrix interpolating between PTR grids
%
% See BIE2D/utils/perispecinterp for convention. Uses N calls to that for now.

I = nan(Nout,Nin);
for i=1:Nin
  v = zeros(Nin,1); v(i) = 1;  % build interp matrix
  I(:,i) = perispecinterp(v,Nout);
end
