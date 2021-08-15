function q = qfs3d_create(b,interior,lpker,srcker,tol,o)
% QFS3D_CREATE  kernel-indep QFS setup, torus- or sphere-like global quadr, 3D.
%
% q = qfs3d_create(b,interior,lpker,srcker,tol,o)
%  sets up source and check surfaces, and dense matrices for computing
%  QFS density from a given surface layer potential density, global quadr
%  (either double-PTR for torus or rings-over-interval for sphere).
%
% Inputs:
%  b = global torus-like or sphere-like boundary surface struct from BIE3D
%  interior = false (exterior eval) or true (interior eval)
%  lpker = density function handle, expecting interface [A An] = lpker(t,s)
%          where t = target pointset (or surf), s = source surf,
%                A = potential evaluation matrix, An = target-normal eval mat.
%          May also be cell array of function handles, in which case returns q
%          a cell array.
%  srcker = QFS source function handle of same interface as lpker.
%  tol = desired acc of eval (usused for o.surfmeth='d')
%  o = opts struct
%       o.verb = 0,1,... verbosity
%       o.srctype = monopoles, or CFIE (for Helm),... [not implemented]
%       o.factor = 's' (trunc SVD, most stable)
%                  'q' (attempt at trunc QR, fails; or as slow as SVD for p-QR),
%                  'B' (do CPQR backsolve w/ E matrix each time; slow, reliable)
%                  'b' (do a LU backsolve with E matrix each time, like 'l')
%                  'l' (piv-LU Gauss elim. Unstab: large soln norm if ill-cond)
%                  'r' (trunc rank-revealing randURV, a la Martinsson / Demmel)
%                  'u' (trunc randUTV, needs matlab interface to randutv)
%       o.surfmeth = 'd' given dists, needs o.param = [srcfac,dsrc,bfac,dchk]
%       o.surfmeth = 'c' confocal ellipsoid (assuming surf is ellipsoid a,b,c),
%                         needs o.param = [a,b,c,srcfac,dsrc,bfac,dchk]
%                    'a' auto-chosen dists based on maxh & tol, needs
%                         o.param = [srcfac]
%       o.dscale sets local distance scaling function: 'c' unity (const dists)
%                    'v' (auto-chosen variable dists based on maxh & tol, needs
%                         o.param = [srcfac])
%       o.I  passes in interpmat I, avoids its calc
%
% Outputs: q = QFS struct (object), or cell array of such, containing fields:
%  s - QFS source surf with s.x nodes, s.w weights, and s.nx normals.
%  qfsco - function handle returning QFS src density coeffs from bdry density
%          (also works for stacks of densities as cols)
%
% With no arguments, detailed self-test is done.

% Todo:
% * include HQRRP factorization once have matlab interface

% Expts on factor='r': normratio for 'l' 6e4,  'r' q=0 2e3, 'r' q=1 6e2,
% vs 's' 1e2, for the shape=2 case. Shows randURV not perfect, but close.
% randURV q=1 is 2x faster than SVD on i7, 3x on xeon.

% Barnett 8/21/19, based on 2D qfs_create. Sphere topo 9/5/19, multi-lp 9/6/19
% reporting 8/14/21
if nargin==0, test_qfs3d_create; return; end
if nargin<6, o=[]; end
if ~isfield(o,'verb'), o.verb = 0; end
if ~isfield(o,'factor'), o.factor = 'l'; end
if ~isfield(o,'surfmeth'), o.surfmeth='d'; o.param=[1.0,0.15,2.0,0.15]; end
N = b.N; Nu=b.Nu; Nv=b.Nv;                     % nodes on input surf

% get biggest inter-node spacing anywhere on original bdry
maxh = max(b.hmax); if o.verb, fprintf('max h = %.3g\n',maxh); end
% crude max anisotropy of surf quadr: (should be cond of [xu;xv] or sth..)
if o.verb, fprintf('max h ratio = %.3g\n',max(b.hmax./b.hmin)); end

% QFS source surf
srcfac = o.param(1);
if o.surfmeth=='a'
  o.param(2) = (maxh/srcfac)/(2*pi) * log(1/tol);    % auto set src dist
end
ds = -sign_from_side(interior) * o.param(2);         % src dist
s = constshiftbdry(b,ds,srcfac,o);
if o.verb
  doverh = srcfac*abs(ds)/maxh;
  fprintf('bfs min(d/h) ~ %.3g \t(corresp tol %.2g)\n',doverh,exp(-2*pi*doverh))
  fprintf('bfs max(d/h) ~ %.3g \t(>5 => bad LU soln norm)\n',ds./min(s.hmin))
end
s.w = ones(size(s.w));      % unit src weights
%if o.verb, fprintf('min(s.sp)/min(b.sp)=%.3g\n',min(s.sp)/min(b.sp)); end
% *** use for self-int test?

% QFS check (collocation) surf: ensure no further than ratio bound
if o.surfmeth=='a'                % dc will be check surf dist...
  dc = ds*(1-log(eps)/log(tol));  % opp sgn from ds. dc:ds ratio not David 5.7-M
  bfac = maxh*log(1/eps) / (2*pi*abs(dc));
  if bfac < 1.0, bfac=1.0; dc = sign(dc)*maxh*log(1/eps)/(2*pi); end  % decr dc
elseif o.surfmeth=='d'
  bfac = o.param(3); dc = sign_from_side(interior) * o.param(4);
end
c = constshiftbdry(b,dc,1.0,o);

bf = setupsurfquad(b,bfac*[max(Nu),Nv]);        % make upsampled physical surf
doverh = bfac*abs(dc)/maxh;
if o.verb, fprintf('cfb min d/maxh ~ %.3g \t(corresp tol %.2g)\n',doverh,exp(-2*pi*doverh)); end
if o.verb, fprintf('QFS N=[%3d,%3d] tol=%6.3g\tsrc fac=%.2f,d=%6.3f\t  bdry fac=%.2f,d=%6.3f\n',max(Nu),Nv,tol,srcfac,ds,bfac,dc); end
numlps = numel(lpker);
for i=1:numlps, q{i}.b = b; q{i}.bf = bf; q{i}.s = s; q{i}.c = c; end
if o.verb>4, qfs_show(q{1}); drawnow; end

ttot = tic; tic          % fill some big matrices...
if numlps==1, lpker = {lpker}; end      % pack so can use as cell arrays
if bfac~=1.0             % bdry upsampling (usual case)
  if ~isfield(o,'I')
    I = surfinterpmat(bf,b);              % mat to upsample bf <- b
  else, I = o.I; o=rmfield(o,'I'); end    % assume user passed in I
  if o.verb>1, fprintf('\tfill I\t\t%.3g s\n',toc); end, tic
  for i=1:numlps
    K = lpker{i}(c,bf);  % eval at check from desired upsampled layer pot (c<-bf)
    if o.verb>1, fprintf('\tfill K\t\t%.3g s\n',toc); end, tic
    cfb{i} = K*I;        % matrix giving check vals from original bdry dens
    if o.verb>1, fprintf('\tcfb=K*I\t\t%.3g s\n',toc); end, tic
  end
  q{1}.I = I;            % pass out
  clear K I
else                     % no bdry upsampling
  for i=1:numlps
    cfb{i} = lpker{i}(c,bf);  % eval at check from desired layer pot (c<-bf)
    if o.verb>1, fprintf('\tfill cfb\t\t%.3g s\n',toc); end, tic
  end  
end
E = srcker(c,s);       % fill c<-s mat (becomes fat if src upsamp from invalid)
                       % (differs from David who keeps E N*N, then src upsamp)
if o.verb>1, fprintf('\tfill E\t\t%.3g s\n',toc); end, tic
% Now factor the E matrix (apart from slow meth='b')...
reps = 1e-14;          % relative eps to set rank truncation (not for LU)
if o.factor=='s'       % trunc SVD - stablest, 20x slower than LU :(
  squarify = 0 && (diff(size(E))~=0);            % 0 manual override (since bad)
  if squarify          % make svd square, proj to lowest src modes (not faster!)
    Is = surfinterpmat(s,b);                       % mat upsamples N to N_src
    if o.verb>1, fprintf('\tfill Is\t\t%.3g s\n',toc); end, tic
    Es = E*Is;
    if o.verb>1, fprintf('\tEs=E*Is\t\t%.3g s\n',toc); end, tic
    [U,S,V] = svd(Es);
    r = sum(diag(S)>reps*S(1,1)); S = diag(S); S = S(1:r); iS = 1./S;  % r=rank
    if o.verb>1, fprintf('\tsvd(Es)\t\t%.3g s (rank=%d, minsig=%.3g)\n',toc,r,S(end)); end, tic
  else                 % full (poss rect) E
    [U,S,V] = svd(E);
    r = sum(diag(S)>reps*S(1,1)); S = diag(S); S = S(1:r); iS = 1./S;  % r=rank
    if o.verb>1, fprintf('\tsvd(E)\t\t%.3g s (rank=%d, minsig=%.3g)\n',toc,r,S(end)); end, tic
  end
  Q2 = V(:,1:r)*diag(iS); if squarify, Q2 = Is*Q2; end   % 2nd factor
  for i=1:numlps
    q{i}.Q1 = U(:,1:r)'*cfb{i}; q{i}.Q2 = Q2;  % 1st factor
    q{i}.qfsco = @(dens) q{i}.Q2*(q{i}.Q1*dens);  % func evals coeffs from dens
  end
  if o.verb>1, fprintf('\tQ1,Q2\t\t%.3g s\n',toc); end
elseif o.factor=='l'   % LU, David's pref.
  if diff(size(E))==0  % square case
    [L,U,P] = lu(E);
    if o.verb>1, fprintf('\tLU(E)\t\t%.3g s\n',toc); end
    for i=1:numlps
      q{i}.qfsco = @(dens) U\(L\(P*(cfb{i}*dens)));  % func taking dens to co
      %q{i}.Q2 = inv(U); q{i}.Q1 = L\(P*cfb{i});  % square (srcfac=1), not bkw stab ? only to 1e-6, etc
    end
  else                 % rect case, David's projecting to NxN then LU
    Is = surfinterpmat(s,b);                       % mat upsamples N to N_src
    if o.verb>1, fprintf('\tfill Is\t\t%.3g s\n',toc); end, tic
    [L,U,P] = lu(E*Is);
    if o.verb>1, fprintf('\tLU(E*Is)\t%.3g s\n',toc); end
    %q.Q2 = Is*inv(U); q.Q1 = L\(P*cfb);           % Q1,Q2, not bkw stab: avoid
    for i=1:numlps
      q{i}.qfsco = @(dens) Is*(U\(L\(P*(cfb{i}*dens)))); % func taking dens->co
    end
  end
elseif o.factor=='q'   % QR: plain no more stable than LU max(d/h)>5, p-QR is
  if diff(size(E))==0  % square case
    [Q,R,P] = qr(E);   % rank-revealing pivoted, slower than SVD!
    %[Q,R] = qr(E);     % plain QR, pair w/ hack below...
    r = sum(abs(diag(R))>reps*abs(R(1,1)));
    R = R(1:r,1:r); Q = Q(:,1:r); P = P(:,1:r);
    if o.verb>1, fprintf('\tpQR(E)\t\t%.3g s (rank=%d)\n',toc,r); end
    for i=1:numlps
      % now try a hack: truncate the plain QR, filling w/ 0s the rest
      %q{i}.qfsco = @(dens) [R\(Q'*(cfb{i}*dens));zeros(size(E,1)-r,size(dens,2))];
      q{i}.qfsco = @(dens) P*(R\(Q'*(cfb{i}*dens)));      % p-QR case
    end
  else                 % rect case, David's projecting to NxN
    % *** not needed, since QR a failure.
  end
elseif o.factor=='r'   % randURV (from Gunnar / Demmel), then trunc solve
  if diff(size(E))==0  % square case
    G = randn(size(E));
    G = E*(E'*G);   % q=1 : lowers norm in ill-cond case, rel to q=0
    %G = E*(E'*G);   % q=2
    [V,~] = qr(G);
    [U,R] = qr(E*V);
    r = sum(abs(diag(R))>reps*abs(R(1,1)));
    R=R(1:r,1:r); U=U(:,1:r); V=V(:,1:r);
    if o.verb>1, fprintf('\trandURV(E)\t%.3g s (rank=%d)\n',toc,r); end, tic
    for i=1:numlps
      q{i}.qfsco = @(dens) V*(R\(U'*(cfb{i}*dens))); % func taking dens to co
    end
  end
elseif o.factor=='u'   % randUTV, then trunc solve. Real E only (no Helmholtz)
  [U,T,V] = randutv(E);
  r = sum(abs(diag(T))>reps*abs(T(1,1)));
  T=T(1:r,1:r); U=U(:,1:r); V=V(:,1:r);
  if o.verb>1, fprintf('\trandUTV(E)\t%.3g s (rank=%d, minT=%.3g)\n',toc,r,abs(T(r,r))); end, tic
  for i=1:numlps
    q{i}.qfsco = @(dens) V*(T\(U'*(cfb{i}*dens))); % func taking dens to co
  end
elseif o.factor=='b'  % do a backsolve each call (uses LU, not bkw stab)
  for i=1:numlps
    q{i}.qfsco = @(dens) E\(cfb{i}*dens);        % func taking dens to co
  end
elseif o.factor=='B'  % do a backsolve each call (slow~SVD! guaranteed bkw stab)
  for i=1:numlps
    q{i}.qfsco = @(dens) linsolve(E,cfb{i}*dens,struct('RECT',true));
  end
end
if o.verb, fprintf('QFS (N=%d,Ns=%d,Nf=%d) total setup %.3g s\n',N,s.N,bf.N,toc(ttot)); end
if numlps==1, q = q{1}; end                        % don't ship a cell array


% ............................ helper functions .....................
function c = constshiftbdry(b,d,fac,o)
% Create quasi-parallel fac-bi-upsampled closed surf from given surf b.
% Uses normal displacement, using analytic functions in b.
% d = const distance. d>0 is in interior (convention matches 2D).
if nargin<4, o=[]; end
%o.minunodes=1;           % remove sources on little rings (default 8)
c = setupsurfquad(b,fac*[max(b.Nu) b.Nv],o);   % note max handles both topo's
c.x = c.x - bsxfun(@times,c.nx,d);   % use d to scale the n vectors
%c.x = circshift(c.x,1);             % *** hack to rot sphere src only!
c = rmfield(c,{'Z','Zu','Zv','xu','xv','sp'});   % leave old nx, w, hmin, hmax
% Notes: latter are nearly right if d << min rad of curvature.

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
verb = 2;                          % 0,1,2,3,4... (& passed to qfs3d_create)
shape = 3;                         % 0: plain torus, 1: cruller, etc (see below)
o.surfmeth = 'a';
tol = 1e-8;                        % tol used in meth 'a' etc
a = 1.0; b = 0.5;                  % baseline torus-like shape params
if shape==0, disp('plain torus double PTR quadr test:')
  N = 40*[2 1]; o.param = [1.0,0.25,1.5,0.2];  % meth='d': srcfac,ds,bfac,dc
  b = modulatedtorus(a,b);
elseif shape==1, disp('cruller double PTR quadr test:')
  b = cruller(b,0.1,5,3);          % replaces b
  N = 72*[2 1]; o.param = [1.0,0.08,2.0,0.08]; % 20s; 1e-5; but nrms 1e5,1e7
  %N = 90*[2 1]; o.param = [2.0,0.04,2.0,0.07];  % low nrm (<10) but srcfac^2 extra src pts!
  %N = 90*[2 1]; o.param = [1.0,0.07,2.0,0.07]; % 1min; 1e-6,1e-5; nrms 1e4,1e5
  b = modulatedtorus(a,b);
elseif shape==2, disp('bent torus double PTR quadr test:')
  b = benttorus(b,0.3,2);          % replaces b
  N = 60*[2 1]; o.param = [1.0,0.16,2.0,0.12];   % bfac needs to be 2 here to get 1e-6 in DLP
  %N = 60*[2 1]; o.param = [1.0,0.2,2.0,0.2];  % demo for 'd' of LU bad 1e4 nrm
  b = modulatedtorus(a,b);
elseif shape==3, disp('sphere interval x PTR test:');
  b = ellipsoid(1,1,1);
  N = 40*[2 1]; o.param = [1.0,0.2,1.5,0.3];   % tried srcfac=1.5
elseif shape==4, disp('ellipsoid interval x PTR test:');
  b = ellipsoid(0.8,1.3,2);
  N = 70*[2 1]; o.param = [1,0.17,2.0,0.17];  %surfmeth='d'
end
%o.minunodes = 16;      % helps sphere-like poles?
b = setupsurfquad(b,N,o);
interior = false;
for lp='D' %'SD', lp             % .... loop over layer pot types
  if lp=='S',     lpker = @Lap3dSLPmat; lpfun = @slpfun;   % lpfun for adaptive
  elseif lp=='D', lpker = @Lap3dDLPmat; lpfun = @dlpfun;
  end
  srcker = @Lap3dSLPmat;             % either S/D fine for Laplace
  o.verb = verb; o.factor = 'u';     % E factor method
  q = qfs3d_create(b,interior,lpker,srcker,tol,o);
  if b.topo=='t'
    densfun = @(u,v) 1+cos(u+.4)+sin(3*u - 2*v + 2 + cos(u+v+1));  % doubly-peri
    vlo=0.0; vhi=2*pi;                                             % v domain
  elseif b.topo=='s'
    %sth = @(v) sqrt(1-v.^2);  % P_1^1 assoc Legendre, so S^2 smooth at poles...
    %densfun = @(u,v) 1+v-v.^2 + 2.0*sin(u).*sth(v) + sin(2*u+1).*v.*sth(v).^2;
    fR3 = @(x,y,z) x.^2+y-0.5+exp(z)/2;  % smooth in R3 
    densfun = @(u,v) fR3(cos(u).*sqrt(1-v.^2), sin(u).*sqrt(1-v.^2), v); % restrict from R3 to S2
    vlo=-1.0; vhi=1.0;                                             % v domain
  end
  dens = fun2dquadeval(densfun,b)';  % node samples of density, col vec
  dists = [1e-6 0.7];                % dists from bdry to test, must be row vec
  u0=.1; v0=.2;    % params to base targ line on; avoid 0 for adaptive's sake
  [~,ny0] = speedfun(b,u0,v0);                                % unit normal
  trg.x = b.Z(u0,v0) - ny0*sign_from_side(interior)*dists;    % test targets
  if verb>3 && lp=='S'
    figure; qfs_show(q); plot3(trg.x(1,:),trg.x(2,:),trg.x(3,:),'k*');
    title(sprintf('int=%d',interior)); showsurffunc(b,dens); title('dens');
    drawnow; end
  ufar = lpker(trg,b) * dens;        % far field (smooth) rule
  ufar(1) = nan;                     % bad for near, ignore
  uada = 0*ufar;                     % adaptive integration of analytic func...
  warning('off','MATLAB:integral:MaxIntervalCountReached');
  for i=1:numel(uada)                % integrand = kernel * speed * dens...
    integrand = @(u,v) lpfun(b,trg.x(:,i),u,v).*speedfun(b,u,v).*densfun(u,v);
    tic; atol = 1e-10;               % for 1e-12 takes several secs :(
    uada(i) = integral2(integrand,0,2*pi,vlo,vhi,'abstol',atol,'reltol',atol);
    fprintf('adaptive integral2 for targ %d in %.3g s\n',i,toc)
    n=50; [uu vv] = ndgrid((0:n-1)/n*2*pi,vlo+(.5:n-.5)/n*(vhi-vlo));
    MAD = mean(abs(integrand(uu(:),vv(:))));                  % crude estimate
    if verb, fprintf('\tintegrand cancellation metric (MAD/mean) %.2g\n',...
      (2*pi)*(vhi-vlo)*MAD/uada(i)), end  % 1 if no cancellation, >1 otherwise
  end
  tic; co = q.qfsco(dens); t=toc;    % do QFS
  if verb>2, showsurffunc(q.s,co); title('QFS co'); end
  normratio = norm(co)/norm(dens.*b.w');    % for SLP, co growth fac due to QFS
  fprintf('QFS get src (normratio:%.3g), in %.3g s\n',normratio,t)
  uqfs = srcker(trg,q.s) * co;       % sum QFS srcs
  if verb, fprintf('\t\t near trg\t far trg\n')
    fprintf('adaptive   \t'); fprintf('%15.10f\t',uada)
    fprintf('\nnative (plain)\t'); fprintf('%15.10f\t',ufar)
    fprintf('\nQFS        \t'); fprintf('%15.10f\t',uqfs)
  end
  fprintf('\nnative far rel err (vs adapt):\t%.3g\n',abs((ufar(2)-uada(2))/uada(2)))
  fprintf('QFS far rel err (vs native):\t%.3g\n',abs((uqfs(2)-ufar(2))/ufar(2)))
  fprintf('QFS close rel err (vs adapt):\t%.3g\n',abs((uqfs(1)-uada(1))/uada(1)))
  if verb>4 && lp=='D'
    B = srcker(b,q.s);                % maps QFS src to bdry vals.
    disp('testing GMRES conv...');  % ext:poor since has small eigval~0 in 1/2+D
    [~,flag,relres,iter,resvec] = gmres(@(x) B*q.qfsco(x),dens,prod(N),tol,300)
    if verb>5, tic
      Q = q.qfsco(eye(prod(N)));    % send in all poss dens: Q maps dens to co
      fprintf('QFS get full Q mat in %.3g s\n',toc)
      A = B*Q;                          % the on-surf DLP eval operator
      %svd(A) % note: for exterior case, is a single 0 eigval; interior not.
      disp('computing eigvals...');
      lam = eig(A); figure; plot(lam,'+'); axis equal
      hold on; plot(-0.5*sign_from_side(interior),0,'r*');   % accum pt of spec?
    end
  end
end                        % ....
