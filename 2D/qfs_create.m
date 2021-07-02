function q = qfs_create(b,interior,lpker,srcker,tol,o)
% QFS_CREATE.  kernel-indep QFS setup, single closed global-PTR curve in 2D
%
% q = qfs_create(b,interior,lpker,srcker,tol,o)
%  sets up source and check boundaries, and dense matrices for computing
%  QFS source vector from a given boundary layer potential density.
%
% Inputs:
%  b = boundary curve struct in BIE2D format.
%  interior = false (exterior eval) or true (interior eval)
%  lpker = density function handle, expecting interface [A An] = lpker(t,s)
%          where t = target pointset (or bdry), s = source bdry,
%                A = potential evaluation matrix, An = target-normal eval mat.
%          If onsurf (QFS-B): used only for self-eval, needs JR term (+-lim).
%  srcker = QFS source kernel function handle of same interface as lpker.
%  tol = desired acc of eval
%  o = options struct:
%       o.verb = 0,1,... verbosity  (0=silent, >3 plots geom).
%       o.curvemeth = 'i' imaginary displ, 'n' normal displ, '2' 2nd-order displ
%       o.srcffac = extra multiplicative p/N srcfac (default 1.0).
%       o.srcfixd = if present enforce (override) imag displacement of src curve
%                   (needs correct sign), otherwise uses auto-choice (default).
%       o.chkfac = enforce check upsampling factor (QFS-D only), otherwise auto.
%       o.factor = dense factorization method:
%                  's' (trunc SVD, 2 mats), 'l' (LU), 'q' (QR; not working),
%                  'n' (naive pseudoinv, single mat).
%       o.onsurf = 0 use off-surf check pts (default), 1 use on-surf self-eval
%                  (1 assumes lpker(s,s) on-surf self-evaluates correctly).
%                  QFS-B is onsurf=1, QFS-D is onsurf=0.
%       o.extrarow = 1 add ncomp row(s) to linear system to fix total strength
%                  to total dens (handles bdry logcap=1), or 0 (default).
%
% Outputs: QFS struct (object) q containing fields:
%  s - QFS source curve with s.x nodes, s.w weights, and s.nx normals.
%  qfsco - function handle returning QFS src density coeffs from bdry density
%          (also works for stacks of densities as cols)
%  Q1, Q2 - two matrix factors of sfb matrix, maybe needed for qfsco
%           so that qfssrcdens = Q2*(Q1*givenbdrydens). LU doesn't use them.
%
% With no arguments, self-test is done.
%
% Notes: 1) in the self-test, adaptive for DLP sucks (can only get 1e-9 for
%           1e-6 dist), so Cauchy scheme still useful for close testing.
%           For Stokes 1e-4 is closest dist :(. Manas & I made better way using
%           stable ker evals and local interp of (x-y) and (x-y).ny, to do.

% Barnett 8/15/19. ncomp for Sto 3/19/21, src/chk logic 6/18/21.
if nargin==0, test_qfs_create; return; end
if nargin<6, o=[]; end          % defaults...
if ~isfield(o,'verb'), o.verb = 0; end
if ~isfield(o,'onsurf'), o.onsurf = 0; end
if ~isfield(o,'srcffac'), o.srcffac=1.0; end
if ~isfield(o,'factor'), o.factor = 's'; end
if ~isfield(o,'curvemeth'), o.curvemeth='i'; end
if ~isfield(o,'extrarow'), o.extrarow=0; end
N = b.N; if mod(N,2), error('b.N must be even for now!'); end    % # user nodes
if tol>1, error('tol should be <1!'); end

% QFS source: curv choice, then p (# src)...
tola = tol;                               % method-adjusted tol
if o.curvemeth=='n', tola = tol/10; end   % equiv David's FF fudge fac ~ 0.3
if ~isfield(o,'srcfixd')    % auto-choose source curve displ imds
  srcsgn = -sign_from_side(interior);
  imds = srcsgn * log(1/tola)/N;          % imag displ of src, for tol@p=N
  imdbig = 2.0*imds;                    % purely for checking self-int
  if intersectfun(b,imdbig,o)==-1       % then |d| > |d_int|/2, maybe not safe
    if ~interior, dmin = 0; dmax = imdbig; else, dmin = imdbig; dmax = 0; end
    [imdint, nsteps] = bisect(@(d) intersectfun(b,d,o), dmin, dmax, 1e-3); % crude tol for src loc only
    badimdfrac = 1.0;               % must be >0.5 due to defn of imdbig above
    imdbad = badimdfrac*imdint;             % borderline bad imag displ
    if abs(imds)>abs(imdbad)
      if o.verb, fprintf('\t qfs_create, src self-int: changing imag displ %.3g to %.3g...\n',imds,imdbad); end
      imds = imdbad;
    end
  end
else
  imds = o.srcfixd;           % user has to override w/ correct sign too :)
  if o.verb && intersectfun(b,imds,o)==-1, fprintf('qfs_create: o.srcfixd=%.3g self-intersects!\n',imds); end
end
ptol = abs(1/imds)*log(1/tola);        % criterion for # src to get adj tol
srcfac = o.srcffac * max(1.0, ptol/N); % will control p = # srcs (never p<N)
if o.verb && N*srcfac<0.95*ptol      % here 0.95 fudge factor avoids many msgs
  fprintf('qfs_create: srcfac=%.3g not pred enough for adj tol %.3g at imds=%.3g! (pred err=%.3g)\n',srcfac, tola, imds, exp(-srcfac*N*abs(imds)));
end
s = shiftedbdry(b,imds,srcfac,o);      % make src curve (2 params: imds, srcfac)
q.b = b; q.s = s;                      % copy out

if o.onsurf     % basic QFS-B, no chk pts needed...
  if o.verb, fprintf('QFS-B N=%3d tol=%5.3g\tsfac=%.2f (p=%d), d=%6.3f\n',N,tol,srcfac,s.N,imds); end
else   % off-surf: QFS check (colloc) curve, specify by its imag shift (imd)..
  imd = imds*(1-log(eps)/log(tola));    % <0, use in/out ratio, not David 5.7-M
  if ~isfield(o,'chkfac'), o.chkfac=1.1*srcfac; end   % > # src, for spec(A) Sto
  valid = false; imdo=imd;
  while ~valid                    % move in check curve until valid...
    c = shiftedbdry(b,imd,o.chkfac,o);
    valid = isempty(selfintersect(real(c.x),imag(c.x)));
    if ~valid
      if o.verb>2, fprintf('\tqfs_create, chk: imd=%.3g self-intersects, reducing...\n',imd); end
      imd = imd/1.05;                % *** crude, should update to bisect
    end
  end
  if o.verb>1 && imdo~=imd, fprintf('chk: adjusted imag dist from %.3g to %.3g...\n',imdo,imd); end
  bfac = log(1/eps)/abs(imd) / N;     % bdry c<-bf upsampling: NB emach not tol!
  if bfac<1.0, bfac=1.0; end          % since why bother
  Nf = ceil(bfac*N/2)*2;               % fine bdry, insure even
  bf = setupquad(b,Nf);                % fine (upsampled) bdry object
  q.bf=bf; q.c=c;
  if o.verb, fprintf('QFS-D N=%3d tol=%5.3g\tsfac=%.2f (p=%d) d=%6.3f\t  bfac=%.2f,d=%6.3f\n',N,tol,srcfac,s.N,imds,bfac,imd); end
end

if o.verb>3, figure(17); clf; qfs_show(q); axis equal tight; drawnow; end   % geom

if o.onsurf            % simpler QFS-B
  cfb = lpker(b,b);    % assumes singular on-surf eval (incl JR +- limit), "A"
  ncomp = size(cfb,1)/N;        % # vector components in the kernel
  E = srcker(b,s);
else
  K = lpker(c,bf);     % eval at check from desired upsampled layer pot (c<-bf)
  ncomp = size(K,2)/Nf;         % # vector components in the kernel
  I = perispecinterpmat(Nf,N);  % mat to upsample one cmpt by bfac
  I = kron(eye(ncomp), I);      % upsample all cmpts
  cfb = K*I;           % matrix giving check vals from original bdry dens
  E = srcker(c,s);     % fill c<-s mat (becomes fat if src upsamp from invalid)
end                    % (differs from David who keeps E N*N, then src upsamp)
if o.extrarow
  srctotrows = s.w'; bdrytotrows = b.w';                   % Laplace case
  % (note in the paper the srctotrow is ones, but we have weights in the QFS
  %  rep not notated in the paper)
  if ncomp==2, srctotrows = [s.w',0*s.w'; 0*s.w',s.w'];    % Stokes (loses nr
    bdrytotrows = [b.w',0*b.w'; 0*b.w',b.w']; end          %  digits but why?)
  E = [E; srctotrows]; cfb = [cfb; bdrytotrows]; 
end

% Now factor the E matrix...
if o.factor=='s'       % trunc SVD - guaranteed, but slow
  reps = 1e-15;        % relative eps to set rank truncation
  [U,S,V] = svd(E);
  minsingval = min(diag(S));
  r = sum(diag(S)>reps*S(1,1)); S = diag(S); S = S(1:r); iS = 1./S;  % r=rank
  if o.verb>2, fprintf('o.factor=s: E is %dx%d, minsigval=%.3g, epstrunc=%.3g, r=%d\n',size(E,1),size(E,2),minsingval,reps,r); end
  q.Q2 = V(:,1:r)*diag(iS); q.Q1 = U(:,1:r)'*cfb; % the 2 factors
  q.qfsco = @(dens) q.Q2*(q.Q1*dens);             % func evals coeffs from dens
  if 0                % experimental only (spit out Qnul)...
    if o.onsurf
      q.Qnul = U(:,r+1:end);                          % keep onb for Nul E^T
    else
      Unul = U(:,r+1:end);                          % onb for Nul E^T
      Qnul = cfb\Unul;    % dens vecs that cannot be produced
      [q.Qnul,R] = qr(Qnul);  % onb for that subspace
      %diag(R)
      %svd(q.Qnul)
    end
  end
elseif o.factor=='n'     % naive version of 's', unstable since 1 matvec.
  X = E\cfb;                    % *** bad way (math correct)
  q.Q1 = eye(N*ncomp); q.Q2 = X;
  q.qfsco = @(dens) X*dens;     % func evals coeffs from dens, bad way
elseif o.factor=='l'   % David's preferred. Q1,Q2 gets only 1e-9, sq or rect
  if diff(size(E))==0           % square case
    Is = eye(N*ncomp);          % dummy for qfsco below
    [L,U,P] = lu(E); q.Q2 = inv(U); q.Q1 = L\(P*cfb);   % square (srcfac=1)
  else                 % rect case, David's projecting to NxN
    Is = perispecinterpmat(s.N,N);     % mat to upsample from N to N_src
    [L,U,P] = lu(E*Is);
    q.Q2 = Is*inv(U); q.Q1 = L\(P*cfb);           % Q1,Q2 *not* bkw stab: avoid
  end
  q.qfsco = @(dens) Is*(U\(L\(P*(cfb*dens))));    % func evals coeffs from dens
elseif o.factor=='q'   % QR, was not good for fat nor square...
  [Q,R] = qr(E); q.Q2 = pinv(R); q.Q1 = Q'*cfb;   % *** test
  %[Q,R] = qr(E',0); q.Q2 = Q; q.Q1 = (R\cfb')';  % turn fat into tall ** debug
end

if 0  % if the factors Q1,Q2 present, check they solve needed matrix lin sys:
  % when P=eye(N), is as bad as 1e-3 (SLP), or O(1) (DLP) - why?
  % not fixed for ncomp yet
  P = perispecinterpmat(N,round(N/4)*2);   % projector tests low modes only
  fprintf('rel ||cfb P - E Q2 Q1 P||=%.3g\n',norm(cfb*P - E*(q.Q2*(q.Q1*P)))/norm(cfb*P))
end


% ............................ helper functions .....................
function c = shiftedbdry(b,imagd,fac,o)
% Create quasi-parallel fac-upsampled closed curve from given curve b.
% Uses imaginary displacement of complexification of boundary parameterization,
% or simple speed-scaled normal displacement.
% imagd - controls distance (either actual imag shift, or equiv normal shift).
%         imagd>0 is in interior. imagd = (2pi/N)*delta from paper.
% Opts:
% o.curvemeth = 'i' (imag shift in complex param), 'n' (speed * normal displ),
%               '2' (as n but w/ 2nd-order term too, using xpp)
% o.forcenum field forces numerical shift via x,nx, even if analytic avail.
if nargin<4, o=[]; end
if isfield(o,'forcenum'), b = rmfield(b,'Zp'); end  % local

% *** todo: add case making Z,Zp from b.x only (no analytic)
% via fft then rebuilding Z, Zp - see larrycup.m

Nf = round(fac*b.N/2)*2;         % pick new N, insure even
if o.curvemeth=='n' || o.curvemeth=='2' % speed-normal displacement variants
  cf.x = perispecinterp(b.x,Nf); % resample bdry
  cf = setupquad(cf);            % rebuild node info (normals) spectrally
  c.x = cf.x - imagd*cf.sp.*cf.nx;  % now just use moved nodes, plain normal
  if o.curvemeth=='2', c.x = c.x - 0.5*imagd^2*cf.xpp; end  % 2nd-ord imag shift
elseif o.curvemeth=='i'              % imag shift, needs analytic curve defn
  c.Z = @(t) b.Z(t + 1i*imagd);
  c.Zp = @(t) b.Zp(t + 1i*imagd);
end
c = setupquad(c,Nf);

function selfint = intersectfun(b,d,o)   % -1 if bdry b self-int, +1 if not
s = shiftedbdry(b,d,1.0,o);
selfint = -1 + 2*isempty(selfintersect(real(s.x),imag(s.x)));

function qfs_show(q)                 % plot all QFS geom on current fig
b=showcurve(q.b,'k'); s=showcurve(q.s,'r');
if isfield(q,'c'), c=showcurve(q.c,'b'); f=plot(q.bf.x,'g.');  % off-surf
  legend([b,s,c,f],'bdry','QFS source','QFS colloc','fine bdry');
else, legend([b,s],'bdry','QFS source'); end

function h=showcurve(s,c)        % simplified showsegment with color control c
if nargin<2, c='k'; end
if ~isnumeric(c), c = rgb(c); end   % convert char to 3-vector of color
h = plot([s.x; s.x(1)],'.-','color',c);               % curve (h = line handle)
h =h(1);       % just keep first obj for legending
hold on; plot([s.x, s.x+0.05*s.nx].','-','color',c);  % normals
  

% ............................... test function ............................
function test_qfs_create  % basic test at fixed N, vs plain adaptive integration
% Laplace LP eval only
a = .3; w = 5;         % smooth wobbly radial shape params
tol = 1e-12;   % 1e-12 N=380;   1e-16 N=180;
N = 380; b = wobblycurve(1,a,w,N);  % really N here needs to grow as log(1/tol)
o.srcffac = 1.05;                     % enforce some upsampling
figure(10); clf; h = showcurve(b); hold on; interior=true;
meths='n2i'; cols = 'bgr';
for i=1:3, o.curvemeth=meths(i); o.verb = 1;
  q = qfs_create(b,interior,@LapSLP,@LapSLP,tol,o);
  h = [h; showcurve(q.s,cols(i))];
end
axis equal tight; legend(h,'bdry','n','2','i');

verb = 1; o.verb = verb; o.curvemeth = '2';      % 'n' int=1 poor, why?
o.onsurf = 0;      % default 0 (QFS-D), otherwise 1 (QFS-B)
srcker = @LapSLP;  % fine for Laplace; will need CFIE for Helmholtz

for interior = [false true], interior  % ------- loop over topologies
  %o.srcfixd = -0.1*sign_from_side(interior);   % test the override
  bny = @(t) b.Zp(t)./(1i*abs(b.Zp(t)));  % unit bdry normal func
  for lp = 'SD', lp          % .... loop over layer pot types
    if lp=='S'
      lpker = @LapSLP;
      selfker = lpker;    % no JR since Dir BC only
      lpclose = @LapSLP_closeglobal; b.a = 0;     % merely for comparison
      lpfun = @(x,t) log(abs(x-b.Z(t)))/(-2*pi);  % Lap SLP formula
    elseif lp=='D'
      lpker = @LapDLP;
      selfker = @(b,varargin) LapDLP(b,varargin{:}) - 0.5*sign_from_side(interior)*eye(b.N);  % JR for DLP on-surf
      lpclose = @LapDLP_closeglobal;
      lpfun = @(x,t) real(conj(x-b.Z(t)).*bny(t))./(2*pi*abs(x-b.Z(t)).^2);  % Lap DLP formula (dot done in C)
    end
    if o.onsurf
      q = qfs_create(b,interior,selfker,srcker,tol,o);   % needs JR
    else
      q = qfs_create(b,interior,lpker,srcker,tol,o);
    end
    densfun = @(t) 1+sin(2+7*t+cos(t)); %exp(sin(t + 1));  % analytic density wrt param
    dens = densfun(b.t);        % density samples
    dists = [1e-3 0.3]';          % distances from bdry to test, must be col vec
    t0 = -0.1;   % param to base targs on; keep away from 0 for adaptive's sake
    trg.x = b.Z(t0) - bny(t0)*sign_from_side(interior)*dists;    % targets
    if verb>1 && lp=='S', figure; qfs_show(q); plot(trg.x,'k*'); axis equal tight; title(sprintf('int=%d',interior)); axis([.5 1.5 -.5 .5]); drawnow; end
  
    ufar = lpker(trg,b,dens);  % far field (smooth) rule
    ufar(1) = nan;             % bad for near, ignore
    uada = 0*ufar;             % adaptive integration of analytic func...
    warning('off','MATLAB:integral:MaxIntervalCountReached');
    for i=1:numel(trg.x)
      uada(i) = integral(@(t) lpfun(trg.x(i),t).*abs(b.Zp(t)).*densfun(t),0,2*pi,'abstol',1e-14,'reltol',1e-14);    % kernel * speed * dens
    end
    % f = @(t) lpfun(trg.x(1),t).*abs(b.Zp(t)).*densfun(t); figure(2); tt=2*pi*(1:1e5)/1e5; plot(tt,f(tt),'-'); title('integrand');  % nasty integrand
    co = q.qfsco(dens);                        % do QFS
    uqfs = srcker(trg,q.s,co);                 % sum QFS srcs
    side = 'e'; if interior, side='i'; end     % compare to barycentric...
    ucau = lpclose(trg,b,dens,side);
    if verb>1, fprintf('\t\t near trg\t\t far trg\n')
      fprintf('adaptive   \t'); fprintf('%18.16g\t',uada)
      fprintf('\nnative (plain)\t'); fprintf('%18.16g\t',ufar)
      fprintf('\nCauchy     \t'); fprintf('%18.16g\t',ucau)
      fprintf('\nQFS        \t'); fprintf('%18.16g\t',uqfs)
    end
    fprintf('\nQFS far rel err (vs native):\t%.3g\n',abs((uqfs(2)-ufar(2))/ufar(2)))
    fprintf('QFS close rel err (vs Cauchy):\t%.3g\n',abs((uqfs(1)-ucau(1))/ucau(1)))
    fprintf('Cau close rel err vs adaptive:\t%.3g\n',abs((ucau(1)-uada(1))/uada(1)))
  end                   % ....
end                                  % ---------
