function r = gendirbvp_conv(pde,interior,known,qfs,o)
% GENDIRBVP_CONV  test convergence of QFS for single-body BVP int/ext, 3 PDEs
%
% gendirbvp_conv(pde,interior,known,onsurf) tests convergence, makes a figure,
% plots text output.
% Solve one of 3 PDEs, for int or ext Dir BVP, via QFS-D/QFS-B and Kress.
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
%             o.Ns  ***
% Outputs:
%  r        - results struct
%             *** to doc
%
% Notes: Mashup of fig_specBIO.m and GRF_conv.m.  Barnett 3/26/21
if nargin<1, test_gendirbvp_conv; return; end
if nargin<5, o=[]; end
if ~isfield(o,'verb'), o.verb=0; end

% BVP setup
a = .3; w = 5;   % smooth wobbly radial shape params
khelm = 10.0;    % wavenumber (Helm only)
mu=0.7;          % viscosity (Sto only)

% QFS params...
if nargin<4, qfs=[]; end
qfs.verb = o.verb;
if ~isfield(qfs,'tol'), qfs.tol = 1e-10; end
if ~isfield(qfs,'onsurf'), qfs.onsurf = 0; end
if ~isfield(qfs,'curvemeth'), qfs.curvemeth = '2'; end

% BVP data we use
b = wobblycurve(1,a,w,100);    % only to access b.Z
sgn = -1+2*interior;              % David convention (+1 if interior)
if known
  t0 = 0.5;   % bdry param to base the rhs pt src on
  N0=200; imt0 = -log(eps)/N0;   % emach by N0 @ the N-rate or 2N0 @ N/2-rate
  z0 = b.Z(t0 - 1i*sgn*imt0);    % data src pt, given imag dist, sets conv rate
  if pde=='L'                    % (see GRF_conv.m for fx,fy Neu data...)
    f0 = exp(1i*4.0);          % data strength to give data size O(1)
    fholom = @(z) f0./(z-z0);       % holom ext soln, with no log nor const
    f = @(z) real(fholom(z));       % use real part as known Laplace soln
  elseif pde=='H'
    f0 = 2.0;
    f = @(z) f0*besselh(0,k*abs(z-z0));    % point source at z0
  elseif pde=='S'
    % ***

  end
else              % scattering. L:linear pot, H:plane wave, S:shear flow
  
  
end
  
if interior, trg.x = -0.1+0.2i;   % BVP soln tests: far int target point
else trg.x = 1.5-0.5i; end        % far ext target point
nrdist = 1e-4;                    % adaptive DLP dies any closer than 1e-6, sad
s=4.0; trg.x(2) = b.Z(s) + nrdist * (b.Zp(s)/1i)/abs(b.Zp(s));  % near test pt
trg.x=trg.x(:); uex = f(trg.x);   % u_exact vals
if o.verb>1, figure(1); clf; plot([b.x; b.x(1)],'-'); hold on; plot(z0,'r*');
  plot(trg.x, 'k+'); axis equal; end

% PDE-specific setup, the BVP data...
ncomp = 1 + (pde=='S');              % # vector cmpts
eta = 0.0 + ~interior;      % eta = amount of S mix:  ext D+S, int D alone
if pde=='L'
  lpker = @(varargin) LapDLP(varargin{:}) + eta*LapSLP(varargin{:});
  refker = @(varargin) LapDLPpotker(varargin{:}) + eta*LapSLPpotker(varargin{:});  % ref needed for adaptive eval
elseif pde=='H'
  eta = 0.0 - 1i*khelm*(~interior);   % CFIE for ext only
  lpker = @(varargin) HelmDLP(khelm,varargin{:}) + eta*HelmSLP(khelm,varargin{:});
  refker = @(varargin) HelmDLPpotker(khelm,varargin{:}) + eta*HelmSLPpotker(khelm,varargin{:});
elseif pde=='S'
  lpker = @(varargin) StoDLP(varargin{:},mu) + eta*StoSLP(varargin{:},mu);
end
srcker = lpker;    % what we claim about QFS

Ns = 100:30:400;           % convergence in # Nystrom pts ..................
% save stuff...
eA = nan(numel(Ns),1); ed=eA; fd=eA; d1=[eA,eA]; eu=[d1,d1]; kA=d1; % >=1cols
for i=1:numel(Ns); N=Ns(i);
  b = wobblycurve(1,a,w,N);
  selfker = @(varargin) lpker(varargin{:}) - 0.5*sign_from_side(interior)*eye(ncomp*N);  % JR; N inside qfs_create fixed
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
  fhat=fft(rhs); fd(i)=abs(fhat(N/2+1)/max(fhat));   % RHS(f) Fou n/2 decay
  dens = A\rhs; dens0=A0\rhs;          % solves (us and Kress)
  d1(i,:)=[dens0(1), dens(1)];         % dens soln vec at 1st (fixed) node
  kA(i,:) = [cond(A0), cond(A)];
  ed(i) = norm(dens-dens0,inf);
  % (but note this is usually high-freq, so how affects distant soln?)
  cod = q.qfsco(dens);
  u(i,:) = srcker(trg,q.s,cod);        % QFS eval from QFS dens
  u0(i,:) = lpker(trg,b,dens0);        % native eval from Kress dens0
  [~,dens0fun] = perispecinterparb(dens0,nan);   % spectral interpolant
  u0(i,2) = lpevaladapt(trg.x(2), refker, dens0fun, b, 1e-12); % adaptive (slow)
  eu(i,:) = abs([u0(i,:)-uex.',u(i,:)-uex.']);   % 4 pt errs
  fprintf('\td10=%.12g\td1=%.12g\tfd=%.3g ed=%.3g\n',d1(i,1),d1(i,2),fd(i),ed(i))
  fprintf('\teA=%.3e\tK(A0)=%.8g\tK(A)=%.8g\teu=%.3g,%.3g\n',eA(i),kA(i,1),kA(i,2),eu(i,3),eu(i,4))
end

% pass out
r.Ns=Ns; r.eu=eu; r.eA=eA; r.ed=ed; r.d=d1; r.kA=kA; r.A=A; r.A0=A0; r.fd=fd;



%%%%%%%%%%%%%%%%%%%%%
function test_gendirbvp_conv
pde='L';
interior=1;
known=1;
qfs.tol = 1e-12;
%qfs.onsurf = 1;
%qfs.srcfac=1.0;
o.verb = 1;
r = gendirbvp_conv(pde,interior,known,qfs,o);
N=r.Ns;
figure(2);clf; semilogy(N,r.eu(:,1),'b+-', N,r.eu(:,2),'b.-', N,r.eu(:,3),'k+-', N,r.eu(:,4),'k.-', N,r.fd,'g.-', N,r.ed,'ro-', N,r.kA,'--');
legend('u0 far plain','u0 nr adap','u far QFS','u nr QFS', 'RHS Fou decay','dens err','cond A0','cond A')
set(gca,'ylim',[1e-15 1e3])
xlabel('n'); hline(qfs.tol)
title(sprintf('%s: int=%d known=%d tol=%.3g',pde,interior,known,qfs.tol))
%lam=eig(r.A); lam0=eig(r.A0);
%figure; imagesc(r.A-r.A0); axis equal; colorbar
