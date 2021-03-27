% simpler test of spec of BIO from QFS, various PDEs, ext Dir BVP.
% Barnett 3/19/21
clear
a = .3; w = 5;         % smooth wobbly radial shape params
tol = 1e-10;        % qfs
interior = 1; % 1     % 0=ext,1=int  (note: int Sto has nullspace)
o.verb=1;
o.onsurf=0;           % 1=QFS-B, 0=QFS-D
o.srcfac = 1.4;  % 'auto'; % fix, or auto checks src curve self-int
o.chkfac = 1.6;
o.curvemeth='2';
pde = 'S';        % PDE: Lap or Helm or Sto

ncomp = 1;              % default # vector cmpts
eta = 0.0 + ~interior;  % default mix of S (ext needs it)
if pde=='L'
  lpker = @(varargin) LapDLP(varargin{:}) + eta*LapSLP(varargin{:});  % D+S
elseif pde=='H'
  k = 10;              % wavenumber
  eta = 0.0 - 1i*k*(~interior);   % CFIE for ext only
  lpker = @(varargin) HelmDLP(k,varargin{:}) + eta*HelmSLP(k,varargin{:});
elseif pde=='S'
  mu=0.7;       % viscosity
  lpker = @(varargin) StoDLP(varargin{:},mu) + eta*StoSLP(varargin{:},mu); % D+S
  ncomp = 2;
end
srcker = lpker;    % what we claim about QFS

fig=figure;
Ns = 100:50:500;
for i=1:numel(Ns); N=Ns(i); %if N==Ns(end), o.verb=3; end
  b = wobblycurve(1,a,w,N);
  selfker = @(varargin) lpker(varargin{:}) - 0.5*sign_from_side(interior)*eye(ncomp*N);  % JR; N is always same inside qfs_create
  if o.onsurf, q = qfs_create(b,interior,selfker,srcker,tol,o);      % QFS-B
  else, q = qfs_create(b,interior,lpker,srcker,tol,o); end           % QFS-D
  B = srcker(b,q.s);   % s2b
  A = (B*q.Q2)*q.Q1;   % QFS nyst (kinda redundant for QFS-B, but useful test)
  A0 = selfker(b,b);     % gold-standard Kress nyst
  Aerr = max(abs(A(:)-A0(:)));
  if pde=='S' && interior, nvec = [real(b.nx);imag(b.nx)];  % nor = nullvec?
    wnvec = nvec.*[b.w;b.w];      % weighted
    %norm(A*nvec), norm(A0*nvec)   % not in R nullspace!!
    %norm(wnvec'*A0)  % is 1e-16, left nullvec known
    %norm(wnvec'*A)  % is 1e-6, left nullvec not good
    %A(1,:) = A(1,:) + nvec'; A0(1,:) = A0(1,:) + nvec'; end  % 1s-mat fix row
    A(:,1) = A(:,1) + wnvec; A0(:,1) = A0(:,1) + wnvec; end  % 1s-mat fix col
  rhs = sin(17*real(b.x)); if ncomp==2, rhs = [rhs; cos(6*imag(b.x))]; end % fake data
  x = A\rhs; x0=A0\rhs;
  %[x(1), x0(1)]   % check dens conv: not stable for int Stokes
  denserr = norm(x-x0,inf);    % fake solve, difference in dens
  % but note this is usually high-freq, so how affects distant soln?
  fprintf('\tAerr=%.3g\tdenserr=%.3g\tK(A)=%g\tK(A0)=%g\n',Aerr,denserr,cond(A),cond(A0))
  lam=eig(A); lam0=eig(A0);
  %figure(fig);
  plot(N,real(lam),'rx', N,real(lam0),'kx'); hold on; %set(gca,'xlim',[min(Ns),max(Ns)]); drawnow
end
%figure; imagesc(A-A0); axis equal; colorbar
