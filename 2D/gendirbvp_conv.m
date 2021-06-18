function [r g] = gendirbvp_conv(pde,interior,known,qfs,o)
% GENDIRBVP_CONV  test convergence of QFS on single-body Dir int/ext BVP, 3 PDEs
%
% gendirbvp_conv(pde,interior,known,onsurf) tests convergence, makes a figure,
% plots text output.
% Solve one of 3 PDEs, for int or ext Dirichlet BVP, via QFS-D/QFS-B and Kress.
% Check soln conv + spectrum of QFS A and of Kress A (called A0).
%
% Inputs:
%  pde      - 'L' Laplace, 'H' Helmholtz, or 'S' Stokes
%  interior - 0=ext, 1=int
%  known    - data: 0=scatt (use self-conv), 1=src on non-phys side (known)
%  qfs      - optional struct for QFS params w/ optional fields:
%             qfs.tol - tolerance
%             qfs.onsurf - BIE system matrix filling: 1=QFS-B, 0=QFS-D (default)
%                 etc
%  o        - optional struct controlling tests
%             o.verb = 0,1,2...
%             o.Ns = list of N values to test convergence over
%             o.grid = struct with 1d lists grid.x, grid.y to eval prod grid
% Outputs:
%  r        - results struct with fields... (when columns, are for the Ns vals)
%             Ns - row of N values used
%             d1 - density soln at 1st fixed node (2 cols): Kress, QFS
%             ed - max over bdry of density diffs btw Kress & QFS
%             eu - pt pot errs (4 cols): far Kress, nr K+adap, far QFS, nr QFS
%             fd - Fourier decay ratio at N/2 for the RHS
%             kA - (2 cols) cond(A_Kress), cond(A_QFS)
%             eA - max abs elementwise A_Kress - A_QFS
%  g - (only if o.grid input), plot grid object with fields...
%             us - u scatt or rep (using the last Ns value)
%             ui - incident or known
%             press, presi - p scatt and p incident (or known), for Sto only
%             x  - target points as C-#s used
%             ii - indices in std meshgrid (array) from grid.x, grid.y
%             q - QFS struct used (for last Ns value)
%
% Notes: Mashup of fig_specBIO.m and GRF_conv.m.  Barnett 3/26/21, 6/15/21.
if nargin<1, test_gendirbvp_conv; return; end
if nargin<5, o=[]; end
if ~isfield(o,'verb'), o.verb=0; end

% BVP test case setup
a = .3; w = 5;   % smooth wobbly radial shape params
khelm = 20.0;    % wavenumber (Helm only)
mu=0.7;          % viscosity (Sto only)

% QFS params...
if nargin<4, qfs=[]; end
qfs.verb = o.verb;
if ~isfield(qfs,'tol'), qfs.tol = 1e-10; end
if ~isfield(qfs,'onsurf'), qfs.onsurf = 0; end
if ~isfield(qfs,'curvemeth'), qfs.curvemeth = '2'; end    % 2nd-ord displ pts

% BVP data we use
b = wobblycurve(1,a,w,100);    % only to access b.Z
sgn = -1+2*interior;              % David convention (+1 if interior)
if known
  t0 = 0.5;   % bdry param to base the rhs pt src on
  N0=250; imt0 = -log(eps)/N0;   % emach by N0 @ the N-rate or 2N0 @ N/2-rate
  z0 = b.Z(t0 - 1i*sgn*imt0);    % data src pt, given imag dist, sets conv rate
  if pde=='L'                    % (see GRF_conv.m for fx,fy Neu data...)
    f0 = exp(1i*4.0);            % strength to give data size O(1), phase orient
    %fholom = @(z) f0./(z-z0);    % holom ext soln, with no log nor const
    fholom = @(z) 1.0*log(z-z0);     % ext soln w/ log (net charge), but zero offset still
    f = @(z) real(fholom(z));    % use real part as known Laplace soln
  elseif pde=='H'
    f = @(z) 2.0 * besselh(0,khelm*abs(z-z0));    % point source at z0
  elseif pde=='S'
    f0 = [0.6;0.8];              % source strength force vector
    f = @(z) StoSLPvelker(mu,z,z0,NaN) * f0;     % Stokeslet at z0, strength f0
    fpres = @(z) StoSLPpresker(mu,z,z0,NaN) * f0;
  end
else              % scattering. f gives bdry data, so is -u_inc on bdry
  if interior, error('interior scatt not implemented!'); end
  if pde=='L'                    % linear potential
    ang = pi/7;
    f = @(z) real(exp(-1i*ang)*z);
  elseif pde=='H'                % plane wave
    ang = pi/7;
    f = @(z) exp(1i*khelm*real(exp(-1i*ang)*z));
  elseif pde=='S'                % backgnd flow
    %f = @(z) 0.6*[-real(z);imag(z)]; fpres = @(z) 0*z;  % stagnation point flow
    f = @(z) [imag(z);0*z]; fpres = @(z) 0*z;  % shear flow
    %f = @(z) 0.6*[-imag(z);real(z)]; fpres = @(z) 0*z;  % pure rot flow
    %f = @(z) [imag(z).^2;0*z]; fpres = @(z) 2*mu*real(z); % leftwards Poisseuil
    if norm(applyStokesPDEs(@(z) [f(z);fpres(z)], 1.2-1.5i, mu, 1e-4))>1e-6, warning('f,fpres not a Stokes soln!'); end
  end
end

ncomp = 1 + (pde=='S');           % # vector cmpts dep on PDE type
if interior, trg.x = -0.1+0.2i;   % BVP soln tests: far int target point
else trg.x = 1.5-0.5i; end        % far ext target point (for int known)
nrdist = 1e-4;                    % adaptive DLP dies any closer than 1e-6, sad
s=4.0; trg.x(2) = b.Z(s) -sgn*nrdist * (b.Zp(s)/1i)/abs(b.Zp(s));  % near test pt
trg.x=trg.x(:); uex = cmpak(f(trg.x),ncomp);   % u_exact vals (or pack as C-#)
if o.verb>1, figure(1); clf; plot([b.x; b.x(1)],'-'); hold on;
  if exist('z0','var'), plot(z0,'r*'); end
  plot(trg.x, 'k+'); axis equal; drawnow; end
  

% PDE-specific setup, the BVP data...
eta = 0.0 + ~interior;      % eta = amount of S mix:  ext D+S, int D alone
if pde=='L'
  lpker = @(varargin) LapDLP(varargin{:}) + eta*LapSLP(varargin{:});
  refker = @(varargin) LapDLPpotker(varargin{:}) + eta*LapSLPpotker(varargin{:});  % ref needed for adaptive eval
elseif pde=='H'
  eta = 0.0 - 1i*khelm*(~interior);   % CFIE for ext only
  lpker = @(varargin) HelmDLP(khelm,varargin{:}) + eta*HelmSLP(khelm,varargin{:});
  refker = @(varargin) HelmDLPpotker(khelm,varargin{:}) + eta*HelmSLPpotker(khelm,varargin{:});
elseif pde=='S'
  lpker = @(t,s,varargin) StoDLP(t,s,mu,varargin{:}) + eta*StoSLP(t,s,mu,varargin{:});  % arg list mu wrong place :(
  refker = @(varargin) StoDLPvelker(mu,varargin{:}) + eta*StoSLPvelker(mu,varargin{:});
end
srcker = lpker;    % what we claim about QFS: uses the same rep as BIE

% convergence in # Nystrom pts ..................
if isfield(o,'Ns'), Ns=o.Ns; else, Ns = 70:30:500; end
% save stuff...
eA = nan(numel(Ns),1); ed=eA; fd=eA; d1=[eA,eA]; u=d1; u0=d1; eu=[d1,d1]; kA=d1; % >=1cols
for i=1:numel(Ns); N=Ns(i);
  b = wobblycurve(1,a,w,N);
  selfker = @(varargin) lpker(varargin{:}) - 0.5*sign_from_side(interior)*eye(ncomp*N);  % JR for Dirichlet BVP; N inside qfs_create fixed
  if qfs.onsurf, q = qfs_create(b,interior,selfker,srcker,qfs.tol,qfs);  % QFS-B
  else, q = qfs_create(b,interior,lpker,srcker,qfs.tol,qfs); end         % QFS-D
  B = srcker(b,q.s);   % s2b
  A = (B*q.Q2)*q.Q1;   % QFS nyst (kinda redundant for QFS-B, but useful test)
  A0 = selfker(b,b);   % gold-standard Kress nyst
  eA(i) = max(abs(A(:)-A0(:)));   % matrix el err - doesn't have to be small
  if pde=='S' && interior    % int Sto needs nullspace 1s-mat correction...
    nvec = [real(b.nx);imag(b.nx)];  % nor = nullvec?
    wnvec = nvec.*[b.w;b.w];      % weighted
    %norm(A*nvec), norm(A0*nvec)   % not in R nullspace!!
    %norm(wnvec'*A0)  % is 1e-16, left nullvec known
    %norm(wnvec'*A)  % is 1e-6, left nullvec not good
    %A(1,:) = A(1,:) + nvec'; A0(1,:) = A0(1,:) + nvec'; end  % 1s-mat fix row
    A(:,1) = A(:,1) + wnvec; A0(:,1) = A0(:,1) + wnvec;  % 1s-mat fix col
  end
  rhs = f(b.x);                        % Dirichlet data
  fhat=fft(rhs(1:N)); fd(i)=abs(fhat(N/2+1)/max(fhat));  % RHS_1 Fou n/2 decay
  dens = A\rhs; dens0=A0\rhs;          % solves (us and Kress)
  %[svd(A), svd(A0)]
  d1(i,:)=[dens0(1), dens(1)];         % dens soln at 1st (fixed) node, 1st cmpt
  kA(i,:) = [cond(A0), cond(A)];
  ed(i) = norm(dens-dens0,inf);        % could be high-freq, not relevant?
  cod = q.qfsco(dens);                 % get QFS src coeffs
  u(i,:) = cmpak(srcker(trg,q.s,cod),ncomp); % QFS eval (all trg) from QFS dens
  u0(i,:) = cmpak(lpker(trg,b,dens0),ncomp);  % native eval from Kress dens0
  dens0fun = interpfunfromgrid(dens0,ncomp); % scalar or vector
  u0(i,2) = cmpak(lpevaladapt(trg.x(2), refker, dens0fun, b, 1e-12),ncomp); % adaptive (slow)
  %real([u0(i,:), u(i,:)])
  eu(i,:) = abs([u0(i,:)-uex.',u(i,:)-uex.']);   % 4 pt errs
  fprintf('\td10=%.12g\td1=%.12g\tfd=%.3g ed=%.3g\n',d1(i,1),d1(i,2),fd(i),ed(i))
  if known   % (notes removed eA(i,1) showing matrix el diff btw Kress, QFS)
    fprintf('\tK(A0)=%.8g\tK(A)=%.8g\teu0=%.2g,%.2g eu=%.2g,%.2g\n',kA(i,1),kA(i,2),eu(i,1),eu(i,2),eu(i,3),eu(i,4))  % Kress and QFS u errs (far,nr)
  else
    fprintf('\tK(A0)=%.8g\tK(A)=%.8g\t\teu=%.2g,%.2g\n',kA(i,1),kA(i,2),abs(u(i,1)-u0(i,1)),abs(u(i,2)-u0(i,2)))    % QFS errs rel to Kress (far,nr)
  end
end
% if no exact soln known, estim err from final u0 vals...
if ~known, uex = u0(end,:); eu = abs([u0,u] - [uex,uex]); end

% pass out
r.Ns=Ns; r.eu=eu; r.eA=eA; r.ed=ed; r.d1=d1; r.kA=kA; r.A=A; r.A0=A0; r.fd=fd;
r.u=u; r.u0=u0;

g.q = q; if known, g.z0 = z0; end
if isfield(o,'grid')    % eval on 2D grid using last N value.......
  if ~isfield(o.grid, 'x'), o.grid.x=linspace(-2,2,201); end   % quantize 0.02
  if ~isfield(o.grid, 'y'), o.grid.y=linspace(-1.8,1.8,181); end    % "
  g.grid = o.grid;
  [xx yy]=meshgrid(o.grid.x,o.grid.y); zz = xx+1i*yy;
  g.ii = xor(~interior, b.inside(zz));   % indices for plotting
  tp.x = zz(g.ii); g.x=tp.x;             % ext targets for plotting
  g.ui = -f(g.x);                        % u_inc
  tic;                                   % slow eval u_scatt on grid (via qfs-b)
  if pde~='S'
    g.us = srcker(tp,q.s,cod);
  else
    g.presi = -fpres(g.x);               % p_inc
    [uD pD] = StoDLP(tp,q.s,mu,cod);     % by hand split out S+D parts of scatt
    [uS pS] = StoSLP(tp,q.s,mu,cod);
    g.us = uD+eta*uS; g.press = pD+eta*pS;  % (easier than func handle wrappers)
  end
  toc
end
%%%%%%%%%%%%%%%%%

function z = cmpak(v,ncomp)          % helpers:  if ncomp=2 pack pairs as C-#s
if ncomp==1, z=v;
elseif ncomp==2, z=v(1:end/2)+1i*v(end/2+1:end);
else, error('cmpak only for ncomp=1 or 2!');
end

function fun=interpfunfromgrid(y,ncomp)  % gives a func handle (scalar or vec)
% this is for 1d grids, periodic only.
if ncomp==1
  [~,fun] = perispecinterparb(y,nan);   % scalar spectral interpolant
elseif ncomp==2
  N = numel(y)/ncomp;
  [~,fun1] = perispecinterparb(y(1:N),nan);   % vector spectral interpolant
  [~,fun2] = perispecinterparb(y(N+1:end),nan);   % vector spectral interpolant
  fun = @(t) [fun1(t); fun2(t)];                 % stack cmpts
else, error('interpfunfromgrid only for ncomp=1 or 2!');
end

function rhs = applyStokesPDEs(f,x,mu,eps)   % from BIE2D/test/testStokernels
% check if func f (returning 3 components: u_1, u_2, and p, given C-# loc x)
% obeys mu-Stokes PDE. Outputs the RHS (3 components: f_1, f_2, and rho)
if nargin<4, eps=1e-5; end
up = nan(3,5);  % three rows are u1,u2,p at each of 5 stencil pts
up(:,1) = f(x); up(:,2) = f(x-eps); up(:,3) = f(x+eps);
up(:,4) = f(x-1i*eps); up(:,5) = f(x+1i*eps);  % do 5-pt stencil evals
gradp = [up(3,3)-up(3,2); up(3,5)-up(3,4)]/(2*eps);
stencil = [-4 1 1 1 1]'/eps^2;
lapu = up(1:2,:)*stencil;
divu = (up(1,3)-up(1,2) + up(2,5)-up(2,4))/(2*eps);
rhs(1:2) = -mu*lapu + gradp;  % Stokes PDE pair
rhs(3)   =  divu;



%%%%%%%%%%%%%%%%%%%%%
function test_gendirbvp_conv
pde='L';
interior=0;
known=1;
qfs.tol = 1e-12;
qfs.onsurf = 1;  % 1 makes QFS-B: only seems to affect ed (dens err)
qfs.srcffac = 1.05;   % 1.2 for Sto, 1.05 for others
%qfs.chkfac = 1.0*qfs.srcfac;  % 1.1 for Sto
o.verb = 1;
o.grid = [];      % tells to do grid eval
%o.Ns = 400;
[r g] = gendirbvp_conv(pde,interior,known,qfs,o);  % do conv -> results struct
N=r.Ns;

% conv plot ........
figure(2);clf; semilogy(N,r.eu(:,1),'b+-', N,r.eu(:,2),'b.-', N,r.eu(:,3),'k+-', N,r.eu(:,4),'k.-', N,r.fd,'g.-', N,r.ed,'ro-', N,r.kA,'--');
legend('u0 far plain','u0 nr adap','u far QFS','u nr QFS', 'RHS Fou decay','dens err','cond A0','cond A')
set(gca,'ylim',[1e-15 1e3])
xlabel('n'); hline(qfs.tol)
title(sprintf('%s: int=%d known=%d tol=%.3g',pde,interior,known,qfs.tol))
drawnow
%lam=eig(r.A); lam0=eig(r.A0);
%figure; imagesc(r.A-r.A0); axis equal; colorbar

if isfield(o,'grid')  % ...... 2D image plot (using g struct out) ..............
figure(3);clf; colormap(jet(256)); up = nan*g.ii;  % plot vals
if pde=='L'
  u0 = 2.0; up(g.ii) = g.us + (1-known)*g.ui;   % u or utot (or NaN)
  contourf(g.grid.x,g.grid.y,up,linspace(-u0,u0,20));
elseif pde=='H'
  u0 = 2.0; up(g.ii) = g.us + (1-known)*g.ui;   % u or utot (or NaN)
  surf(g.grid.x,g.grid.y,real(up),'alphadata',~isnan(up));   % NaN->white
elseif pde=='S'
  u0 = 2.0;
  pp = nan*g.ii; pp(g.ii) = g.press + (1-known)*g.presi;   % p or ptot (or NaN)
  contourf(g.grid.x,g.grid.y,pp,linspace(-u0,u0,20)); hold on;
  % pack vel into C#, subsample the grid for vel arrow plot...
  m = numel(g.x);  % how many targs, their indices in the g.x array
  dxq = 0.1; iiq = find(abs(mod(real(g.x)+dxq/2,dxq)-dxq/2)+abs(mod(imag(g.x)+dxq/2,dxq)-dxq/2)<1e-6);
  up = g.us(iiq)+1i*g.us(m+iiq) + (1-known)*(g.ui(iiq)+1i*g.ui(m+iiq));
  quiver(real(g.x(iiq)),imag(g.x(iiq)), real(up),imag(up), 2.0,'k-')
end
view(2); shading interp; caxis(u0*[-1 1]); grid off;
title(sprintf('%s: int=%d known=%d',pde,interior,known),'interpreter','latex');
hold on; plot([g.q.b.x; g.q.b.x(1)],'k-');
text(0,0,1,'$\Omega$','interpreter','latex','fontsize',20);  % NB z=1 3D lift
text(0.5,0.5,1,'$\partial\Omega$','interpreter','latex','fontsize',20);
axis xy equal tight; axis([min(g.grid.x),max(g.grid.x),min(g.grid.y),max(g.grid.y)]);
end

%keyboard
