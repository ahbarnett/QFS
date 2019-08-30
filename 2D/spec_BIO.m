function lam = spec_BIO(tol,verb,curvemeth,N,interior)
% SPEC_BIO.  Plot spectrum of discrete bdry integral operator given by QFS
%
% lam = spec_BIO(tol,verb,curvemeth,N,interior)
%  outputs and plots spectrum of A (on-surf operator) appoximated by QFS at
%  tolerance tol (which should be no smaller than say 1e-14).
%  ext Laplace DLP only for now. All args are optional (first 2 can be empty).
%  Also prints cond(A).
%
% See also: anim_spec_BIO
%
% Barnett 8/29/19
if nargin<1 || isempty(tol), tol = 1e-10; end   % don't make too close to emach!
if nargin<2 || isempty(verb), verb = 0; end
o = []; o.factor='l';             % QFS opts
o.verb = (verb>1);
if nargin>=3, o.curvemeth=curvemeth; end
if nargin<5, interior = true; end

k = 0;                 % wavenumber (0 for now)
a = .3; w = 5;         % smooth wobbly radial shape params

srcker = @LapSLP;                 % choose QFS src rep
sgn = sign_from_side(interior);

if nargin<4
  N = ceil(40*log10(1/tol));        % N getting QFS tol for starfish
end
b = wobblycurve(1,a,w,N);
q = qfs_create(b,interior,@LapDLP,srcker,tol,o);       % QFS approx to DLP
Q = q.qfsco(eye(N));              % send in all poss densities: Q maps dens to co
B = srcker(b,q.s);                % maps QFS src to bdry vals.
A = B*Q;                          % the on-surf DLP eval operator
fprintf('N=%5d: cond(A_{QFS}) = %.3g\n',N,cond(A))
%svd(A) % note: for exterior case, is a single 0 eigval; interior not.
lam = eig(A);
plot(lam,'+'); axis equal
hold on; plot(-0.5*sgn,0,'r*'); % supposed accumulation pt of spec
% since A = +-1/2 + D is what QFS approximates


