function lam = spec_BIO(tol,verb,curvemeth,N,interior,lp)
% SPEC_BIO.  Plot spectrum of discrete 2nd-kind bdry integral op found by 2D QFS
%
% lam = spec_BIO(tol,verb,curvemeth,N,interior)
%  outputs and plots spectrum of A (on-surf operator) appoximated by QFS with
%  tolerance tol (which should be no smaller than say 1e-14), for a simple
%  curve in 2D.
%
%  Args interior (true=int, false=exterior) and lp (='D', 'S') allow 4 cases:
%    int D:  op should approximate -1/2 + D, well cond
%    int S:  op should approximate 1/2 + D^T, nullity 1
%    ext D:  op should approximate 1/2 + D, nullity 1
%    ext S:  op should approximate -1/2 + D^T, well cond
%
%  All args are optional, and args 1,2,4 can be empty.
%  Also prints cond(A).
%  Returns: all eigenvalues of A.
%
% Purpose: to check that 1-sided QFS gives a good spectrum accumulating at +-1/2
%
% See also: anim_spec_BIO
%
% Examples of four 2nd-kind operators above:
%  figure; spec_BIO(1e-10,0,'i',[],true,'D');
%  figure; spec_BIO(1e-10,0,'i',[],true,'S');
%  figure; spec_BIO(1e-10,0,'i',[],false,'D');
%  figure; spec_BIO(1e-10,0,'i',[],false,'S');
%
% Barnett 8/29/19
if nargin<1 || isempty(tol), tol = 1e-10; end   % don't make too close to emach!
if nargin<2 || isempty(verb), verb = 0; end
o = []; o.factor='l';             % QFS opts
o.verb = (verb>1);
if nargin>=3, o.curvemeth=curvemeth; end
if nargin<5, interior = true; end
if nargin<6, lp='D'; end

k = 0;                 % wavenumber (0 for now)
a = .3; w = 5;         % smooth wobbly radial shape params

srcker = @LapSLP;                 % choose QFS src rep
JRsgn = 1; if lp=='S', JRsgn=-1; end  % jump rel sign

if nargin<4 || isempty(N)
  N = ceil(40*log10(1/tol));        % N getting QFS tol for starfish
end
b = wobblycurve(1,a,w,N);
if lp=='D'
  q = qfs_create(b,interior,@LapDLP,srcker,tol,o);       % QFS approx to DLP
  B = srcker(b,q.s);              % maps QFS src to bdry vals.
elseif lp=='S'
  q = qfs_create(b,interior,@LapSLP,srcker,tol,o);       % QFS approx to SLP
  [~,B] = srcker(b,q.s);          % maps QFS src to bdry nderiv.
end  
Q = q.qfsco(eye(N));              % send in all poss dens: Q maps dens to co
A = B*Q;                          % the approx on-surf JR + D or S'=D^T operator
fprintf('N=%5d: cond(A_{QFS}) = %.3g\n',N,cond(A))
%svd(A) % note: for ext D (int S) case, is a single 0 eigval; else well cond.
lam = eig(A);
plot(lam,'+'); axis equal; hold on;
plot(-0.5*sign_from_side(interior)*JRsgn,0,'r*');   % supposed accum pt of spec
% since A = +-1/2 + D is what QFS approximates, in D case


