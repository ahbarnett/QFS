function q = qfs_create(b,interior,lpker,srcker,tol,o)
% QFS_CREATE.  kernel-indep QFS setup, single closed global-PTR curve in 2D
%
% q = qfs_create(b,interior,lpker,srcker,tol,o)
%  sets up source and check boundaries, and dense matrices for computing
%  QFS density from a given boundary layer potential density.
%
% Inputs:
%  b = boundary curve struct from BIE2D.
%  interior = false (exterior eval) or true (interior eval)
%  lpker = density function handle, expecting interface [A An] = lpker(t,s)
%          where t = target pointset (or bdry), s = source bdry,
%                A = potential evaluation matrix, An = target-normal eval mat.
%  srcker = QFS source function handle of same interface as lpker.
%  tol = desired acc of eval
%  o = opts struct
%       o.verb = 0,1,... verbosity
%       o.srctype = monopoles, or CFIE (for Helm),... [not implemented]
%       o.curvemeth = 'i' imaginary displacement, 'n' normal displacement.
%       o.factor = 's' (trunc SVD), 'l' (LU), 'q' (QR; not working).
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
% Notes: 1) adaptive for DLP sucks (can only get 1e-9 for 1e-6 dist), so Cauchy
% still useful for close testing.

% Barnett 8/15/19
if nargin==0, test_qfs_create; return; end
if nargin<6, o=[]; end
if ~isfield(o,'verb'), o.verb = 0; end
if ~isfield(o,'factor'), o.factor = 's'; end
if ~isfield(o,'curvemeth'), o.curvemeth='i'; end
N = b.N;                      % nodes on input bdry

% QFS source curve
FF = 0.0;                       % David's fudge "factor" (0.37 gains 1 digit)
if o.curvemeth=='n', FF = 0.5;  end   % since worse curves (0.5 good for ext)
imds = -sign_from_side(interior) * log(1/tol)/N;   % initial guess, tol @ fac=1
valid = false;
while ~valid                    % move in & upsample src curve until valid...
  srcfac = (log(1/tol)+2*pi*FF)/abs(imds) / N;
  s = shiftedbdry(b,imds,srcfac,o);
  valid = isempty(selfintersect(real(s.x),imag(s.x)));
  %valid=1;   % *** force no src upsampling
  if o.verb, fprintf('min(s.sp)/min(b.sp)=%.3g\n',min(s.sp)/min(b.sp)); end
  if o.curvemeth=='n', valid = valid & min(s.sp)>0.5*min(b.sp); end  % David's
  if ~valid, imds = imds/1.1; end   % bring closer, from which fac will be set
end

% QFS check (collocation) curve, specify by its imag shift (imd)
imd = imds*(1-log(eps)/log(tol)); % <0, use in/out ratio instead of David 5.7-M
valid = false;
while ~valid                    % move in check curve until valid...
  c = shiftedbdry(b,imd,1.0,o);
  valid = isempty(selfintersect(real(c.x),imag(c.x)));
  if ~valid, imd = imd/1.1; end
end
bfac = log(1/eps)/abs(imd) / N;       % bdry c<-bf upsampling: NB emach not tol!
if bfac<1.0, bfac=1.0; end            % since why bother
%bfac = 3.0*srcfac;  % old plain fixed bdry upsampling
Nf = ceil(bfac*N/2)*2;               % fine bdry, insure even
bf = setupquad(b,Nf);  % fine (upsampled) bdry
q.b = b; q.bf = bf; q.s = s; q.c = c;  % copy out (q.s is only crucial)
if o.verb, fprintf('QFS N=%4d tol=%6.3g\tsrc fac=%.2f,d=%6.3f\t  bdry fac=%.2f,d=%6.3f\n',N,tol,srcfac,imds,bfac,imd); end
if o.verb>2, figure; qfs_show(q); axis equal tight; drawnow; end

K = lpker(c,bf);       % eval at check from desired upsampled layer pot (c<-bf)
I = perispecinterpmat(Nf,N);  % mat to upsample by bfac
cfb = K*I;             % matrix giving check vals from original bdry dens
E = srcker(c,s);       % fill c<-s mat (becomes fat if src upsamp from invalid)
                       % (differs from David who keeps E N*N, then src upsamp)
% Now factor the E matrix...
if o.factor=='s'       % trunc SVD - guaranteed, but slow
  reps = 1e-15;        % relative eps to set rank truncation
  [U,S,V] = svd(E);
  r = sum(diag(S)>reps*S(1,1)); S = diag(S); S = S(1:r); iS = 1./S;  % r=rank
  q.Q2 = V(:,1:r)*diag(iS); q.Q1 = U(:,1:r)'*cfb; % the 2 factors
  q.qfsco = @(dens) q.Q2*(q.Q1*dens);             % func evals coeffs from dens
elseif o.factor=='l'   % David's preferred. Q1,Q2 gets only 1e-9, sq or rect
  if diff(size(E))==0
    Is = eye(N);       % dummy for qfsco below
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
  P = perispecinterpmat(N,round(N/4)*2);   % projector tests low modes only
  fprintf('rel ||cfb P - E Q2 Q1 P||=%.3g\n',norm(cfb*P - E*(q.Q2*(q.Q1*P)))/norm(cfb*P))
end


% ............................ helper functions .....................
function c = shiftedbdry(b,imagd,fac,o)
% Create quasi-parallel fac-upsampled closed curve from given curve b.
% Uses imaginary displacement of complexification of boundary parameterization,
% or simple normal displacement.
% imagd - controls distance (either actual imag shift, or equiv normal shift).
%         imagd>0 is in interior.
% Opts:
% o.curvemeth = 'i' (imag shift in complex param), 'n' (speed * normal displ).
% o.forcenum field forces numerical shift via x,nx, even if analytic avail.
if nargin<4, o=[]; end
if isfield(o,'forcenum'), b = rmfield(b,'Zp'); end  % local

% *** todo: add case making Z,Zp from b.x only (no analytic)
% via fft then rebuilding Z, Zp - see larrycup.m

Nf = round(fac*b.N/2)*2;         % pick new N, insure even
if o.curvemeth=='i'              % imag shift analytic curve defn
  c.Z = @(t) b.Z(t + 1i*imagd);
  c.Zp = @(t) b.Zp(t + 1i*imagd);
elseif o.curvemeth=='n'          % plain speed-normal displacement
  cf.x = perispecinterp(b.x,Nf); % resample bdry
  cf = setupquad(cf);            % rebuild node info spectrally
  c.x = cf.x - imagd*cf.sp.*cf.nx;  % now just use moved nodes
end
c = setupquad(c,Nf);

function qfs_show(q)                 % plot all QFS geom on current fig
b=showcurve(q.b,'k'); s=showcurve(q.s,'r'); c=showcurve(q.c,'b');
f=plot(q.bf.x,'g.');
legend([b,s,c,f],'bdry','QFS source','QFS colloc','fine bdry');

function h=showcurve(s,c)        % simplified showsegment with color control c
h = plot([s.x; s.x(1)],[c '.-']);                     % curve (h = line handle)
hold on; plot([s.x, s.x+0.05*s.nx].',[c '-']);        % normals


% ............................... test function ............................
function test_qfs_create  % basic test at fixed N, vs plain adaptive integration
warning('off','MATLAB:integral:MaxIntervalCountReached');
verb = 1;
a = .3; w = 5;         % smooth wobbly radial shape params
N = 400; b = wobblycurve(1,a,w,N);  % really N here needs to grow as log(1/tol)
for interior = [false true], interior  % ------- loop over topologies
  bny = @(t) b.Zp(t)./(1i*abs(b.Zp(t)));  % unit bdry normal func
  for lp = 'SD', lp          % .... loop over layer pot types
    if lp=='S'
      lpker = @LapSLP;
      lpclose = @LapSLP_closeglobal; b.a = 0;     % merely for comparison
      lpfun = @(x,t) log(abs(x-b.Z(t)))/(-2*pi);  % Lap SLP formula
    elseif lp=='D'
      lpker = @LapDLP;
      lpclose = @LapDLP_closeglobal;
      lpfun = @(x,t) real(conj(x-b.Z(t)).*bny(t))./(2*pi*abs(x-b.Z(t)).^2);  % Lap DLP formula (dot done in C)
    end
    srcker = @LapSLP;  % fine for Laplace; will need CFIE for Helmholtz
    tol = 1e-12;
    o.verb = verb; o.curvemeth = 'i';   % 'n' int=1 poor, needs N larger - why?
    q = qfs_create(b,interior,lpker,srcker,tol,o);
    densfun = @(t) 1+sin(7*t+cos(t)); %exp(sin(t + 1));  % analytic density wrt param
    dens = densfun(b.t);        % density samples
    dists = [1e-3 0.3]';          % distances from bdry to test, must be col vec
    t0 = -0.1;   % param to base targs on; keep away from 0 for adaptive's sake
    trg.x = b.Z(t0) - bny(t0)*sign_from_side(interior)*dists;    % targets
    if verb>1 && lp=='S', figure; qfs_show(q); plot(trg.x,'k*'); axis equal tight; title(sprintf('int=%d',interior)); axis([.5 1.5 -.5 .5]); drawnow; end
  
    ufar = lpker(trg,b,dens);  % far field (smooth) rule
    ufar(1) = nan;             % bad for near, ignore
    uada = 0*ufar;             % adaptive integration of analytic func...
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
