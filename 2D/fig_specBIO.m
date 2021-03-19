% simpler test of spec of BIO from QFS, various PDEs, ext Dir BVP.
% Barnett 3/19/21
clear
a = .3; w = 5;         % smooth wobbly radial shape params
tol = 1e-10;        % qfs
interior = 0; % 1     % 0=ext,1=int  (note: Sto int has nullspace)
o.verb=1;
o.onsurf=1;           % 1=QFS-B, 0=QFS-D
o.srcfac = 1.0; %'auto';
o.curvemeth='2';
pde = 'S';        % Lap or Helm or Sto, 1 letter
ncomp = 1;     % # vector cmpts, by default
eta = 0.0 + ~interior;    % mix of S certainly needed for ext, by default
if pde=='L'
  lpker = @(varargin) LapDLP(varargin{:}) + eta*LapSLP(varargin{:});  % D+S
elseif pde=='H'
  k = 10;
  eta = 0.0 - 1i*k*(~interior);   % CFIE for ext only
  lpker = @(varargin) HelmDLP(k,varargin{:}) + eta*HelmSLP(k,varargin{:});
elseif pde=='S'
  mu=0.7;       % viscosity
  lpker = @(varargin) StoDLP(varargin{:},mu) + eta*StoSLP(varargin{:},mu); % D+S
  ncomp = 2;
end
srcker = lpker;    % what we claim about QFS

figure;
Ns = 100:50:400;
for i=1:numel(Ns); N=Ns(i);
  b = wobblycurve(1,a,w,N);
  selfker = @(varargin) lpker(varargin{:}) - 0.5*sign_from_side(interior)*eye(ncomp*N);  % JR; N is always same inside qfs_create
  if o.onsurf, q = qfs_create(b,interior,selfker,srcker,tol,o);      % QFS-B
  else, q = qfs_create(b,interior,lpker,srcker,tol,o); end           % QFS-D
  B = srcker(b,q.s);   % s2b
  A = (B*q.Q2)*q.Q1;   % QFS nyst (kinda redundant for QFS-B, but useful test)
  A0 = selfker(b,b);     % gold-standard Kress nyst
  Aerr = max(abs(A(:)-A0(:)));
  fprintf('\tAerr=%.3g\tK(A)=%.12g\tK(A0)=%.12g\n',Aerr,cond(A),cond(A0))
  lam=eig(A);
  plot(N,real(lam),'kx'); hold on;
  %drawnow
end
%figure; imagesc(A-A0); axis equal; colorbar
