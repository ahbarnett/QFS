function q = qfs3d_create(b,interior,lpker,srcker,tol,o)
% QFS3D_CREATE.  kernel-indep QFS setup, single torus-like double global-PTR, 3D
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
%       o.factor = 's' (trunc SVD, will be v slow), 'l' (LU, 40x faster)
%
% Outputs: QFS struct (object) q containing fields:
%  s - QFS source surf with s.x nodes, s.w weights, and s.nx normals.
%  qfsco - function handle returning QFS src density coeffs from bdry density
%          (also works for stacks of densities as cols)
%
% With no arguments, self-test is done.

% Barnett 8/21/19, based on 2D qfs_create.
if nargin==0, test_qfs3d_create; return; end
if nargin<6, o=[]; end
if ~isfield(o,'verb'), o.verb = 0; end
if ~isfield(o,'factor'), o.factor = 'l'; end
N = b.N; Nu=b.Nu; Nv=b.Nv;                     % nodes on input surf

aspects = (sqrt(sum(b.xu.^2,1))/b.Nu) ./ (sqrt(sum(b.xv.^2,1))/b.Nv);
if o.verb, fprintf('max h ratio = %.3g\n',max([aspects,1./aspects])); end

% QFS source surf (only node x for now; would need nx for dipole or CFIE...)
srcfac = 1.2;        % hack all for now
d = -sign_from_side(interior) * 0.2;   % hack
s = shiftedbdry(b,d,srcfac,o);
s.w = ones(1,s.N);                     % dummy weights
%if o.verb, fprintf('min(s.sp)/min(b.sp)=%.3g\n',min(s.sp)/min(b.sp)); end
% *** use for self-int test?

% QFS check (collocation) surf (has only node x for now)
dc = sign_from_side(interior) * 0.2;   % hack
c = shiftedbdry(b,dc,1.0,o);

% upsampled surf
bfac = 2.0;
Nuf = ceil(bfac*Nu/2)*2; Nvf = ceil(bfac*Nv/2)*2;  % fine bdry, insure even
bf = setupdoubleptr(b,[Nuf Nvf]);
q.b = b; q.bf = bf; q.s = s; q.c = c;
if o.verb, fprintf('QFS N=[%4d,%4d] tol=%6.3g\tsrc fac=%.2f,d=%6.3f\t  bdry fac=%.2f,d=%6.3f\n',Nu,Nv,tol,srcfac,dc,bfac,d); end

tic; % fill some big matrices...
K = lpker(c,bf);       % eval at check from desired upsampled layer pot (c<-bf)
I = peri2dspecinterpmat([Nuf,Nvf],[Nu,Nv]);  % mat to upsample on surface
cfb = K*I;             % matrix giving check vals from original bdry dens
E = srcker(c,s);       % fill c<-s mat (becomes fat if src upsamp from invalid)
                       % (differs from David who keeps E N*N, then src upsamp)
if o.verb, fprintf('fill cfb and E in %.3g s\n',toc); end
tic; % Now factor the E matrix...
if o.factor=='s'       % trunc SVD - guaranteed  *** replace by faster sq LU?
  reps = 1e-15;        % relative eps to set rank truncation
  [U,S,V] = svd(E);
  r = sum(diag(S)>reps*S(1,1)); S = diag(S); S = S(1:r); iS = 1./S;  % r=rank
  q.Q2 = V(:,1:r)*diag(iS); q.Q1 = U(:,1:r)'*cfb;  % the 2 factors
  q.qfsco = @(dens) q.Q2*(q.Q1*dens);              % func evals coeffs from dens
elseif o.factor=='l'   % David's preferred. Q1,Q2 gets only 1e-9, sq or rect
  if diff(size(E))==0
    Is = eye(N);       % dummy for qfsco below
    [L,U,P] = lu(E); q.Q2 = inv(U); q.Q1 = L\(P*cfb);   % square (srcfac=1)
  else                 % rect case, David's projecting to NxN
    Is = peri2dspecinterpmat([s.Nu s.Nv],[Nu Nv]);   % mat upsamples N to N_src
    [L,U,P] = lu(E*Is);
    q.Q2 = Is*inv(U); q.Q1 = L\(P*cfb);           % Q1,Q2 *not* bkw stab: avoid
  end
  q.qfsco = @(dens) Is*(U\(L\(P*(cfb*dens))));    % func evals coeffs from dens
end
if o.verb, fprintf('factor in %.3g s\n',toc); end



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
b=showsurf(q.b,'k');
o=[]; o.alpha = 0.5; s=showsurf(q.s,'r',o);
o.alpha = 0.1; c=showsurf(q.c,'b',o);
o.normals=0; o.alpha = 0.3; f=showsurf(q.bf,'g',o);
lightangle(45,0);
legend([b,s,c,f],'surf','QFS source','QFS colloc','fine surf');


% ............................... test function ............................
function test_qfs3d_create  % basic test at fixed N, vs plain adaptive integr
warning('off','MATLAB:integral:MaxIntervalCountReached');
verb = 2;
shape = 1;          % 0: plain torus, 1: cruller.
a = 1.0; b = 0.5;       % baseline torus params
if shape==0, disp('plain torus double PTR quadr test:'), N = [40 20];
else, disp('cruller double PTR quadr test:'), N = [60 30];
  b = cruller(b,0.1,5,3);    % replaces b
end
b = setup_torus_doubleptr(a,b,N);
interior = false;
Zn = @(u,v) cross(b.Zu(u,v),b.Zv(u,v));                     % normal bdry func
sp = @(u,v) sqrt(sum(Zn(u,v).^2,1));                        % surf element
bny = @(u,v) Zn(u,v)./sp(u,v);                              % unit normal func
lp='S';
if lp=='S'
  lpker = @Lap3dSLPmat;
  lpfun = @(x,u,v) sqrt(sum((x-b.Z(u,v)).^2,1))/(4*pi);     % Lap SLP formula
elseif lp=='D'
  lpker = @Lap3dDLPmat;
  %lpfun = @(x,t) real(conj(x-b.Z(t)).*bny(t))./(2*pi*abs(x-b.Z(t)).^2);  % Lap DLP formula (dot done in C)
end
srcker = @Lap3dSLPmat;  % fine for Laplace
tol = 1e-6;
o.verb = verb; o.factor = 'l';
q = qfs3d_create(b,interior,lpker,srcker,tol,o);
densfun = @(u,v) 1.0 + sin(3*u - 2*v + cos(u+3*v));         % doubly-periodic
[buu bvv] = ndgrid(b.u,b.v); dens = densfun(buu(:),bvv(:)); % dens at surf nodes
dists = [1e-3 0.3];                % dists from bdry to test, must be row vec
u0=.1; v0=.2;    % params to base targs on; keep away from 0 for adaptive's sake
trg.x = b.Z(u0,v0) - bny(u0,v0)*sign_from_side(interior)*dists;    % targets
if verb>1 && lp=='S', figure; qfs_show(q); plot3(trg.x(1,:),trg.x(2,:),trg.x(3,:),'k*'); title(sprintf('int=%d',interior)); drawnow; end
ufar = lpker(trg,b) * dens;         % far field (smooth) rule, bad for near
uada = 0*ufar;             % adaptive integration of analytic func...
for i=1:numel(uada)
  uada(i) = integral2(@(u,v) lpfun(trg.x(:,i),u,v).*sp(u,v).*densfun(u,v),0,2*pi,0,2*pi,'abstol',1e-10,'reltol',1e-10);    % kernel * speed * dens
    end
tic; co = q.qfsco(dens); fprintf('apply QFS dens in %.3g s\n',toc)    % do QFS
uqfs = srcker(trg,q.s) * co;        % sum QFS srcs
if verb>1, fprintf('\t\tnear trg \t\tfar trg\n')
  fprintf('adaptive   \t'); fprintf('%.16g\t',uada)
  fprintf('\nnative (sm)\t'); fprintf('%.16g\t',ufar)
  fprintf('\nQFS        \t'); fprintf('%.16g\t',uqfs)
end
fprintf('\nQFS far rel err (vs native):\t%.3g\n',abs((uqfs(2)-ufar(2))/ufar(2)))
fprintf('QFS close rel err (vs adaptive):\t%.3g\n',abs((uqfs(1)-uada(1))/uada(1)))
