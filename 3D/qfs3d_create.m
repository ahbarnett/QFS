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
%       o.surfmeth = 'd' (given dists, needs o.param = [srcfac,dsrc,bfac,dchk])
%                    'a' (auto-chosen dists based on maxh & tol, needs
%                         o.param = [srcfac])
%                    'v' (auto-chosen variable dists based on maxh & tol, needs
%                         o.param = [srcfac])
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
if ~isfield(o,'surfmeth'), o.surfmeth='d'; o.param=[2.0,0.15,1.0,0.15]; end
N = b.N; Nu=b.Nu; Nv=b.Nv;                     % nodes on input surf

% get biggest inter-node spacing anywhere on original bdry
h1 = sqrt(sum(b.xu.^2,1))*(2*pi/b.Nu);  % h_1, h_2 at all given surf nodes
h2 = sqrt(sum(b.xv.^2,1))*(2*pi/b.Nv);
maxh = max([h1,h2]); if o.verb, fprintf('max h = %.3g\n',maxh); end
% crude max anisotropy of surf quadr: (should be cond of [xu;xv] or sth..)
as = h1./h2; if o.verb, fprintf('max h ratio = %.3g\n',max([as,1./as])); end

% QFS source surf
srcfac = o.param(1);
if o.surfmeth=='a'
  o.param(2) = (maxh/srcfac)/(2*pi) * log(1/tol);    % auto set src dist
end
ds = -sign_from_side(interior) * o.param(2);
s = constshiftbdry(b,ds,srcfac,o);
doverh = srcfac*abs(ds)/maxh;
if o.verb, fprintf('bfs min d/maxh ~ %.3g \t(corresp tol %.2g)\n',doverh,exp(-2*pi*doverh)); end
s.w = ones(size(s.w));      % dummy weights
%if o.verb, fprintf('min(s.sp)/min(b.sp)=%.3g\n',min(s.sp)/min(b.sp)); end
% *** use for self-int test?

% QFS check (collocation) surf: ensure no further than ratio bound
if o.surfmeth=='a'
  dc = ds*(1-log(eps)/log(tol));  % opp sgn from ds, in/out ratio not David 5.7-M
  bfac = maxh*log(1/eps) / (2*pi*abs(dc));
  if bfac < 1.0, bfac=1.0; dc = sign(dc)*maxh*log(1/eps)/(2*pi); end  % decr dc
elseif o.surfmeth=='d'
  bfac = o.param(3); dc = sign_from_side(interior) * o.param(4);
end
c = constshiftbdry(b,dc,1.0,o);

% upsampled physical surf
Nuf = ceil(bfac*Nu/2)*2; Nvf = ceil(bfac*Nv/2)*2;  % fine bdry, insure even
bf = setupdoubleptr(b,[Nuf Nvf]);
doverh = bfac*abs(dc)/maxh;
if o.verb, fprintf('cfb min d/maxh ~ %.3g \t(corresp tol %.2g)\n',doverh,exp(-2*pi*doverh)); end
if o.verb, fprintf('QFS N=[%3d,%3d] tol=%6.3g\tsrc fac=%.2f,d=%6.3f\t  bdry fac=%.2f,d=%6.3f\n',Nu,Nv,tol,srcfac,ds,bfac,dc); end
q.b = b; q.bf = bf; q.s = s; q.c = c;
  
t = tic; tic; % fill some big matrices...
K = lpker(c,bf);       % eval at check from desired upsampled layer pot (c<-bf)
if o.verb>1, fprintf('\tfill K\t\t%.3g s\n',toc); end, tic
if [Nuf Nvf]==[Nu Nv]  % no src upsampling
  cfb = K; clear K;
else
  I = peri2dspecinterpmat([Nuf,Nvf],[Nu,Nv]);  % mat to upsample on surface
  if o.verb>1, fprintf('\tfill I\t\t%.3g s\n',toc); end, tic
  cfb = K*I;           % matrix giving check vals from original bdry dens
  if o.verb>1, fprintf('\tcfb=K*I\t\t%.3g s\n',toc); end, tic
  clear K I
end
E = srcker(c,s);       % fill c<-s mat (becomes fat if src upsamp from invalid)
                       % (differs from David who keeps E N*N, then src upsamp)
if o.verb>1, fprintf('\tfill E\t\t%.3g s\n',toc); end, tic
% Now factor the E matrix...
if o.factor=='s'       % trunc SVD - guaranteed, slow
  reps = 1e-15;        % relative eps to set rank truncation
  [U,S,V] = svd(E);
  r = sum(diag(S)>reps*S(1,1)); S = diag(S); S = S(1:r); iS = 1./S;  % r=rank
  q.Q2 = V(:,1:r)*diag(iS); q.Q1 = U(:,1:r)'*cfb;  % the 2 factors
  q.qfsco = @(dens) q.Q2*(q.Q1*dens);              % func evals coeffs from dens
elseif o.factor=='l'   % David's preferred. (Q1,Q2 gets only 1e-9, sq or rect)
  if diff(size(E))==0  % square case
    [L,U,P] = lu(E);
    if o.verb>1, fprintf('\tLU(E)\t\t%.3g s\n',toc); end
    %q.Q2 = inv(U); q.Q1 = L\(P*cfb);   % square (srcfac=1), not bkw stab
    q.qfsco = @(dens) U\(L\(P*(cfb*dens)));        % func evals coeffs from dens
  else                 % rect case, David's projecting to NxN
    Is = peri2dspecinterpmat([s.Nu s.Nv],[Nu Nv]);   % mat upsamples N to N_src
    if o.verb>1, fprintf('\tfill Is\t\t%.3g s\n',toc); end, tic
    [L,U,P] = lu(E*Is);
    if o.verb>1, fprintf('\tLU(E*Is)\t%.3g s\n',toc); end
    %q.Q2 = Is*inv(U); q.Q1 = L\(P*cfb);           % Q1,Q2, not bkw stab: avoid
    q.qfsco = @(dens) Is*(U\(L\(P*(cfb*dens))));   % func evals coeffs from dens
  end
end
if o.verb, fprintf('QFS (N=%d) total setup %.3g s\n',N,toc(t)); end



% ............................ helper functions .....................
function c = constshiftbdry(b,d,fac,o)
% Create quasi-parallel fac-bi-upsampled closed surf from given surf b.
% Uses normal displacement, using analytic functions in b.
% d = distance.
% d>0 is in interior, to match 2D.
if nargin<4, o=[]; end
Nuf = ceil(fac*b.Nu/2)*2; Nvf = ceil(fac*b.Nv/2)*2;  % pick new Ns, insure even
c = setupdoubleptr(b,[Nuf Nvf]);
c.x = c.x - bsxfun(@times,c.nx,d);   % use d to scale the n vectors
c = rmfield(c,{'Z','Zu','Zv','xu','xv','sp'});   % leaves nx (nearly right), w

function qfs_show(q)                 % plot all 3D QFS geom on current fig
b=showsurf(q.b,'k');
o=[]; o.normals=0; o.alpha = 1.0; s=showsurf(q.s,'r',o);
o.alpha = 0.05; c=showsurf(q.c,'b',o); o.alpha = 0.3; f=showsurf(q.bf,'g',o);
legend([b,s,c,f],'surf','QFS source','QFS colloc','fine surf');
lightangle(45,0);



% test helpers: adaptive integral2 needs u,v to be any shape (eyeroll)...
function [sp ny] = speedfun(b,u,v)  % b is the surf object, (u,v) the params.
sz = size(u); u = u(:)'; v = v(:)';   % we can handle row-vectorization only
Zn = cross(b.Zu(u,v),b.Zv(u,v));                     % normal bdry func (3*n)
sp = sqrt(sum(Zn.^2,1));                             % surf element
ny = Zn./sp;                      % will vectorize 3xn ./ 1xn to give 3xn
sp = reshape(sp,sz);

function s = slpfun(b,x,u,v)      % bare kernel of SLP, handling u,v any shape
sz = size(u); u = u(:)'; v = v(:)';   % we can handle row-vectorization
s = (1/4/pi) ./ sqrt(sum((x-b.Z(u,v)).^2,1));        % note x is single 3x1
s = reshape(s,sz);

function s = dlpfun(b,x,u,v)      % bare kernel of DLP, handling u,v any shape
sz = size(u); u = u(:)'; v = v(:)';   % we can handle row-vectorization
[~,ny] = speedfun(b,u,v);
s = (1/4/pi)*sum((x-b.Z(u,v)).*ny,1) ./ sqrt(sum((x-b.Z(u,v)).^2,1)).^3;
s = reshape(s,sz);

% ............................... test function ............................
function test_qfs3d_create  % basic test at fixed N, vs plain adaptive integr
verb = 2;
shape = 0;                         % 0: plain torus, 1: cruller, 2: bent torus
o.surfmeth = 'a';
tol = 1e-5;                        % tol used in meth 'a' etc
a = 1.0; b = 0.5;                  % baseline torus shape params
if shape==0, disp('plain torus double PTR quadr test:')
  N = 50*[2 1]; o.param = [1.0,0.15,2.0,0.15];  % meth='d': srcfac,ds,bfac,dc
elseif shape==1, disp('cruller double PTR quadr test:')
  b = cruller(b,0.1,5,3);          % replaces b
  N = 72*[2 1]; o.param = [1.0,0.08,2.0,0.08]; % 20s; 1e-5; but nrms 1e5,1e7
  %N = 90*[2 1]; o.param = [1.0,0.07,2.0,0.07]; % 1min; 1e-6,1e-5; nrms 1e4,1e5
else, disp('bent torus double PTR quadr test:')
  b = benttorus(b,0.3,2);          % replaces b
  N = 60*[2 1]; o.param = [1.0,0.15,2.0,0.15];   % bfac needs to be 2 here to get 1e-6 in DLP
end
b = setup_torus_doubleptr(a,b,N);
interior = false;
for lp='SD' %'SD', lp             % .... loop over layer pot types
  if lp=='S',     lpker = @Lap3dSLPmat; lpfun = @slpfun;
  elseif lp=='D', lpker = @Lap3dDLPmat; lpfun = @dlpfun;
  end
  srcker = @Lap3dSLPmat;             % fine for Laplace
  o.verb = verb; o.factor = 'l';
  q = qfs3d_create(b,interior,lpker,srcker,tol,o);
  densfun = @(u,v) 1+cos(u+.4)+sin(3*u - 2*v + 2 + cos(u+3*v));  % doubly-peri
  [buu bvv] = ndgrid(b.u,b.v); dens = densfun(buu(:),bvv(:)); % dens at surf nodes
  dists = [1e-3 1];                  % dists from bdry to test, must be row vec
  u0=.1; v0=.2;    % params to base targ line on; avoid 0 for adaptive's sake
  [~,ny0] = speedfun(b,u0,v0);                                % unit normal
  trg.x = b.Z(u0,v0) - ny0*sign_from_side(interior)*dists;    % test targets
  if verb>2 && lp=='S', figure; qfs_show(q); plot3(trg.x(1,:),trg.x(2,:),trg.x(3,:),'k*'); title(sprintf('int=%d',interior)); drawnow; end
  ufar = lpker(trg,b) * dens;        % far field (smooth) rule
  ufar(1) = nan;                     % bad for near, ignore
  uada = 0*ufar;                     % adaptive integration of analytic func...
  warning('off','MATLAB:integral:MaxIntervalCountReached');
  for i=1:numel(uada)                % integrand = kernel * speed * dens...
    integrand = @(u,v) lpfun(b,trg.x(:,i),u,v).*speedfun(b,u,v).*densfun(u,v);
    tic; atol = 1e-10;
    uada(i) = integral2(integrand,0,2*pi,0,2*pi,'abstol',atol,'reltol',atol);
    fprintf('adaptive integral2 for targ %d in %.3g s\n',i,toc)
    n=100; [uu vv] = ndgrid((1:n)/n*2*pi,(1:n)/n*2*pi);
    MAD = mean(abs(integrand(uu(:),vv(:))));
    if verb, fprintf('\tintegrand cancellation metric (MAD/mean) %.2g\n',(2*pi)^2*MAD/uada(i)), end  % 1 if no cancellation, >1 otherwise
  end
  tic; co = q.qfsco(dens); fprintf('QFS get src (norm:%.3g), in %.3g s\n',norm(co),toc)    % do QFS
  uqfs = srcker(trg,q.s) * co;       % sum QFS srcs
  if verb, fprintf('\t\t near trg\t far trg\n')
    fprintf('adaptive   \t'); fprintf('%15.10f\t',uada)
    fprintf('\nnative (plain)\t'); fprintf('%15.10f\t',ufar)
    fprintf('\nQFS        \t'); fprintf('%15.10f\t',uqfs)
  end
  fprintf('\nnative far rel err (vs adapt):\t%.3g\n',abs((ufar(2)-uada(2))/uada(2)))
  fprintf('QFS far rel err (vs native):\t%.3g\n',abs((uqfs(2)-ufar(2))/ufar(2)))
  fprintf('QFS close rel err (vs adapt):\t%.3g\n',abs((uqfs(1)-uada(1))/uada(1)))
end                        % ....