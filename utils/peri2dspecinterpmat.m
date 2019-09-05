function I = peri2dspecinterpmat(Nout,Nin)
% PERI2DSPECINTERPMAT.  Return dense matrix interp between double-PTR grids
%
% I = peri2dspecinterpmat(Nout,Nin) returns matrix prod(Nout)-by-prod(Nin)
%  in size.
%
% See BIE3D/utils/peri2dspecinterp for convention.
% Uses N=prod(Nin) calls to that for now, so is O(N^2 ln N), not the best.
% Note the ordering of the 2D arrays in the rows, or cols, of matrix I, is
%  the natural Fortran ordering of matlab.
%
% Without arguments, does self-test.

if nargin==0, test_peri2dspecinterpmat; return; end
I = kron(perispecinterpmat(Nout(2),Nin(2)),perispecinterpmat(Nout(1),Nin(1)));

function I = peri2dspecinterpmat_old(Nout,Nin)
nout=prod(Nout); nin=prod(Nin);    % matrix size
I = nan(nout,nin);
for i=1:nin
  v = zeros(nin,1); v(i) = 1;      % unit vector to send in
  vout = peri2dspecinterp(reshape(v,Nin),Nout);
  I(:,i) = vout(:);
end
% timing: (acceptable)
% tic; I = peri2dspecinterpmat([200 200],[100 100]); toc  % 10 sec on i7

%%%%%%%%%%%%
function test_peri2dspecinterpmat   % compares new (kron) and old. that's it.
Nin = [30 50];
Nout = [20 70];
tic;
I = peri2dspecinterpmat(Nout,Nin);
toc
tic;
I2 = peri2dspecinterpmat_old(Nout,Nin);
toc
fprintf('I size %d-by-%d\n',size(I,1),size(I,2))
norm(I-I2,'fro')    % notice this is *not* small if the symmetrizing differs!

