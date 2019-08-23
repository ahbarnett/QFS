function I = perispecinterpmat(Nout,Nin)
% PERISPECINTERPMAT.  Return dense matrix interpolating between PTR grids
%
% I = perispecinterpmat(Nout,Nin) returns matrix Nout-by-Nin in size.
%
% See BIE2D/utils/perispecinterp for convention.
% Uses N calls to that for now, so is O(N^2 ln N), not the best.

I = nan(Nout,Nin);
for i=1:Nin
  v = zeros(Nin,1); v(i) = 1;       % unit vector to send in
  I(:,i) = perispecinterp(v,Nout);
end
