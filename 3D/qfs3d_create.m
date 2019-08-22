function q = qfs3d_create(b,interior,lpker,srcker,tol,o)
% QFS3D_CREATE.  kernel-indep QFS setup, single toruslike double global-PTR, 3D
%
% q = qfs3d_create(b,interior,lpker,srcker,tol,o)
%  sets up source and check surfaces, and dense matrices for computing
%  QFS density from a given surface layer potential density, global double-PTR.
%
% Inputs:
%  b = global torus-like boundary surface struct from BIE3D
%  interior = false (exterior eval) or true (interior eval)
%  lpker = density function handle, expecting interface [A An] = lpker(t,s)
%          where t = target pointset (or surf), s = source surf,
%                A = potential evaluation matrix, An = target-normal eval mat.
%  srcker = QFS source function handle of same interface as lpker.
%  tol = desired acc of eval
%  o = opts struct
%       o.verb = 0,1,... verbosity
%       o.srctype = monopoles, or CFIE (for Helm),... [not implemented]
%
% Outputs: QFS struct (object) q containing fields:
%  s - QFS source surf with s.x nodes, s.w weights, and s.nx normals.
%  qfsco - function handle returning QFS src density coeffs from bdry density
%          (also works for stacks of densities as cols)
%  Q1, Q2 - two matrix factors of sfb matrix, that are needed for qfsco
%           so that qfssrcdens = Q2*(Q1*givenbdrydens)
%
% With no arguments, self-test is done.

% Barnett 8/21/19, based on 2D qfs_create.
if nargin==0, test_qfs3d_create; return; end
if nargin<6, o=[]; end
if ~isfield(o,'verb'), o.verb = 0; end
if ~isfield(o,'factor'), o.factor = 's'; end
N = b.N; Nu=b.Nu; Nv=b.Nv;                     % nodes on input surf

aspects = (sqrt(sum(b.xu.^2,1))/b.Nu) ./ (sqrt(sum(b.xv.^2,1))/b.Nv);
if o.verb, fprintf('max h ratio = %.3g\n',max([aspects,1./aspects])); end

% QFS source surf
srcfac = 1.2;        % hack all for now
d = -sign_from_side(interior) * 0.2;   % hack
s = shiftedbdry(b,d,srcfac,o);
%if o.verb, fprintf('min(s.sp)/min(b.sp)=%.3g\n',min(s.sp)/min(b.sp)); end
% *** use for self-int test?

% QFS check (collocation) surf
dc = sign_from_side(interior) * 0.2;   % hack
c = shiftedbdry(b,dc,1.0,o);

% upsampled surf
bfac = 2.0;
Nuf = ceil(bfac*Nu/2)*2; Nvf = ceil(bfac*Nv/2)*2;  % fine bdry, insure even
bf = setupdoubleptr(b,[Nuf Nvf]);
q.b = b; q.bf = bf; q.s = s; q.c = c;
if o.verb, fprintf('QFS N=[%4d,%4d] tol=%6.3g\tsrc fac=%.2f,d=%6.3f\t  bdry fac=%.2f,d=%6.3f\n',Nu,Nv,tol,srcfac,dc,bfac,d); end

K = lpker(c,bf);       % eval at check from desired upsampled layer pot (c<-bf)
I = peri2dspecinterpmat([Nuf,Nvf],[Nu,Nv]);  % mat to upsample on surface
cfb = K*I;             % matrix giving check vals from original bdry dens
E = srcker(c,s);       % fill c<-s mat (becomes fat if src upsamp from invalid)
                       % (differs from David who keeps E N*N, then src upsamp)
% Now factor the cfb matrix...
reps = 1e-10;          % relative eps to set rank truncation
if o.factor=='s'       % trunc SVD - guaranteed  *** replace by faster sq LU?
  [U,S,V] = svd(E);
  r = sum(diag(S)>reps*S(1,1)); S = diag(S); S = S(1:r); iS = 1./S;  % r=rank
  q.Q2 = V(:,1:r)*diag(iS); q.Q1 = U(:,1:r)'*cfb;   % the 2 factors
else
end
q.qfsco = @(dens) q.Q2*(q.Q1*dens);                % func evals coeffs from dens


% ............................ helper functions .....................
function c = shiftedbdry(b,d,fac,o)
% Create quasi-parallel fac-bi-upsampled closed surf from given surf b.
% Uses normal displacement, using analytic functions in b.
% d = distance.
% d>0 is in interior, to match 2D.
if nargin<4, o=[]; end
Nuf = ceil(fac*b.Nu/2)*2; Nvf = ceil(fac*b.Nv/2)*2;  % pick new Ns, insure even
c = setupdoubleptr(b,[Nuf Nvf]);
c.x = c.x - bsxfun(@times,c.nx,d);   % use d to scale the n vectors
c = rmfield(c,{'Z','Zu','Zv','xu','xv','nx','sp','w'});

function qfs_show(q)                 % plot all 3D QFS geom on current fig
b=showsurf(q.b,'k'); s=showsurf(q.s,'r');
%c=showsurf(q.c,'b');
%o=[]; o.normals=0; f=showsurf(q.bf,'g',o);
%legend([b,s,c,f],'surf','QFS source','QFS colloc','fine surf');


% ............................... test function ............................
function test_qfs3d_create  % basic test at fixed N, vs plain adaptive integration
warning('off','MATLAB:integral:MaxIntervalCountReached');
verb = 1;
shape = 0;
a = 1.0; b = 0.5;       % baseline torus params
if shape==0
  disp('plain torus double PTR quadr test:')
else
  disp('cruller double PTR quadr test:')
  b = cruller(b,0.1,5,3);    % replaces b
end
N = [80 40];
b = setup_torus_doubleptr(a,b,N);
interior = false; lp='S';
if lp=='S'
  lpker = @Lap3dSLPmat;
  %lpfun = @(x,t) log(abs(x-b.Z(t)))/(-2*pi);  % Lap SLP formula
elseif lp=='D'
  lpker = @Lap3dDLPmat;
  %lpfun = @(x,t) real(conj(x-b.Z(t)).*bny(t))./(2*pi*abs(x-b.Z(t)).^2);  % Lap DLP formula (dot done in C)
end
srcker = @Lap3dSLPmat;  % fine for Laplace
tol = 1e-12;
o.verb = verb;
q = qfs3d_create(b,interior,lpker,srcker,tol,o);
if verb && lp=='S', figure; qfs_show(q); drawnow; end
