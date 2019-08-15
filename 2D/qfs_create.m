function q = qfs_create(b,interior,lpker,srcker,tol,o)
% QFS_CREATE.  kernel-indep QFS setup, single closed global-PTR curve in 2D
%
% q = qfs_mats(b,interior,lpker,srcker,tol,o)
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
%  o = opts struct ...
%       o.srctype = monopoles, or CFIE (for Helm),...
%
% Outputs: qfs struct (object) q containing fields:
%  s - QFS source curve with s.x nodes, s.w weights, and s.nx normals.
%  qfsco - function handle returning QFS src density coeffs from bdry density.
%  Q1, Q2 - two matrices that are needed for getdens
%           so that srcdens = Q2*(Q1*dens)
%
% With no arguments, self-test is done

% Barnett 8/15/19
if nargin==0, test_qfs_create; return; end
if nargin<6, o=[]; end
N = b.N;                      % nodes on input bdry

% QFS source curve
srcfac = 1.0;
FF = 0.3;
alpha = log(1/tol)/(2*pi) + FF; % # h-units away needed for tol b<-s (David's M)
imd = -sign_from_side(interior) * alpha * (2*pi)/(srcfac*N);   % 2pi facs cancel
s = shiftedbdry(b,imd,srcfac);
valid = isempty(selfintersect(real(s.x),imag(s.x))); if ~valid, 'src invalid', end
% *** bring closer if not valid, should also scale the

% QFS check (collocation) curve
bfac = 3.0;  % bdry upsampling  *** fix
alpha = log(1/tol)/(2*pi);   % # hfine-units away needed for tol c<-bf
Nf = ceil(bfac*N);
imd = sign_from_side(interior) * alpha * 2*pi/Nf;  % 2pi facs cancel
c = shiftedbdry(b,imd,1.0);
valid = isempty(selfintersect(real(c.x),imag(c.x))); if ~valid, 'chk invalid', end
fprintf('QFS: srcfac=%.3g, bdryfac=%.3g\n',srcfac,bfac)

bf = setupquad(b,Nf);  % fine (upsampled) bdry
K = lpker(c,bf);       % eval at check from desired upsampled layer pot (c<-bf)
I = perispecinterpmat(Nf,N);  % upsample by bfac
E = srcker(c,s);       % fill c<-s mat
[U,S,V] = svd(E);
reps = 1e-15;
r = sum(diag(S)>reps); S = diag(S); S = S(1:r); iS = 1./S;  % r=rank, trunc
q.Q2 = V(:,1:r)*diag(iS); q.Q1 = U(:,1:r)'*(K*I);  % the 2 factors
q.qfsco = @(dens) q.Q2*(q.Q1*dens);                % func evals coeffs from dens

q.b = b; q.bf = bf; q.s = s; q.c = c;  % copy out


function c = shiftedbdry(b,imagd,fac)
% use imaginary displacement of complexification of boundary parameterization,
% by given imag dist, with given upsampling factor, from boundary object b.
% (imagd>0 is in interior)

Nf = ceil(fac*b.N);
if isfield(b,'Zp')               % wrap analytic curve defn
  c.Z = @(t) b.Z(t + 1i*imagd);
  c.Zp = @(t) b.Zp(t + 1i*imagd);
  c = setupquad(c,Nf);
else
  % *** add case of purely from nodes b.x
end



function a = sign_from_side(interior)  % David's convention
a = -1; if interior, a = 1; end

function qfs_show(q)                 % plot all QFS geom on current fig
b=showcurve(q.b,'k'); s=showcurve(q.s,'r'); c=showcurve(q.c,'b');
f=plot(q.bf.x,'g.');
legend([b,s,c,f],'bdry','QFS source','QFS colloc','fine bdry');

function h=showcurve(s,c)        % simplified showsegment with color control c
h = plot([s.x; s.x(1)],[c '.-']);                     % curve (h = line handle)
hold on; plot([s.x, s.x+0.05*s.nx].',[c '-']);   % normals

%%%%%%%%%
function test_qfs_create  % basic test at fixed N, vs plain adaptive integration
verb = 0;
a = .3; w = 5;         % smooth wobbly radial shape params
N = 240; b = wobblycurve(1,a,w,N);
interior = false;
for lp = 'SD'          % .... loop over layer pot types
  lp
  if lp=='S'
    lpker = @LapSLP;
    lpclose = @LapSLP_closeglobal; b.a = 0;     % merely for comparison
    lpfun = @(x,t) log(abs(x-b.Z(t)))/(-2*pi);  % Lap SLP formula
  elseif lp=='D'
    lpker = @LapDLP;
    lpclose = @LapDLP_closeglobal;
    bny = @(t) b.Zp(t)./(1i*abs(b.Zp(t)));  % unit bdry normal func
    lpfun = @(x,t) real(conj(x-b.Z(t)).*bny(t))./(2*pi*abs(x-b.Z(t)).^2);  % Lap DLP
  end
  srcker = @LapSLP;  % fine for Laplace; will need CFIE for Helmholtz
  tol = 1e-12;
  q = qfs_create(b,interior,lpker,srcker,tol);
  if verb, figure(1); qfs_show(q); axis equal tight; end
  densfun = @(t) exp(sin(t + 1));       % analytic density wrt param
  dens = densfun(b.t);        % density samples
  trg.x = b.x(1) + b.nx(1)*10.^[-6 0]';   % targets in C cplane, must be col vec
  ufar = lpker(trg,b,dens);  % far field (smooth) rule
  uada = 0*ufar;                     % adaptive integration of analytic funcs
  for i=1:numel(trg.x)
    uada(i) = integral(@(t) lpfun(trg.x(i),t).*abs(b.Zp(t)).*densfun(t),0,2*pi,'abstol',1e-15);    % kernel * speed * dens
  end
  co = q.qfsco(dens);
  uqfs = srcker(trg,q.s,co);
  side = 'e'; if interior, side='i'; end          % compare to barycentric for fun
  ucau = lpclose(trg,b,dens,side);
  fprintf('\t\tnear trg \t\tfar trg\n')
  fprintf('adaptive   \t'); fprintf('%.16g\t',uada)
  fprintf('\nnative (sm)\t'); fprintf('%.16g\t',ufar)
  fprintf('\nCauchy     \t'); fprintf('%.16g\t',ucau)
  fprintf('\nQFS        \t'); fprintf('%.16g\t',uqfs)
  fprintf('\nQFS far err (vs native):\t%.3g\n',abs((uqfs(2)-ufar(2))/ufar(2)))
  fprintf('QFS close err (vs Cauchy):\t%.3g\n',abs((uqfs(1)-ucau(1))/ucau(1)))
  fprintf('Cau close err (vs adaptive):\t%.3g\n',abs((ucau(1)-uada(1))/uada(1)))
end