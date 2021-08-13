% 2D ext kappa(A) plot over upsampling parameters, 3 PDEs, on single shape.
% QFS-D Nystrom A. Exterior only.
% Barnett 8/11/21.

clear; %close all
interior = 0;
a = .3; w = 5;   % smooth wobbly radial shape params
curve = @(N) wobblycurve(1.0,a,w,N);    % wrap the discretized bdry struct-maker
%dexpt=0.2; curve = @(N) wobblycurve(exp(dexpt),0,w,N); qfs.srcfixd = dexpt; % everride expt to test src a unit circle thus w/ C_gamma=1: S-rep no probs
%qfs.srcfixd=0.2;
khelm = 20.0;    % wavenumber (Helm only)
mu=0.7;          % viscosity (Sto only)

sgn = -1+2*interior;              % David convention (+1 if interior)
qfs.onsurf = 0;                   % QFS-D

qfs.factor = 's'; qfs.meth='2'; qfs.verb=1;          % QFS meth pars
pdes = 'LHS';                                        % PDE types, 1 char each
pdenam = {'Laplace', sprintf('Helm.%s        k=%g',char(10),khelm), 'Stokes'};
tol = 1e-12;
N = 200;
b = curve(N);           % setup discretization

for ipde=3                                       % which PDEs to test
  pde = pdes(ipde); fprintf('PDE=%s: -----------------------\n',pdenam{ipde})
  ncomp = 1 + (pde=='S');           % # vector cmpts dep on PDE type
      
  if pde=='L'                       % wrap the LPs in PDE-indep way
    SLP = @LapSLP; DLP = @LapDLP;
  elseif pde=='H'
    SLP = @(varargin) HelmSLP(khelm,varargin{:});
    DLP = @(varargin) HelmDLP(khelm,varargin{:});
  elseif pde=='S'
    SLP = @(t,s,varargin) StoSLP(t,s,mu,varargin{:});
    DLP = @(t,s,varargin) StoDLP(t,s,mu,varargin{:});
  end
  
  if pde=='L'          % PDE-dep choice of QFS rep, and of rep to test A for...
    qfsrep = @(varargin) SLP(varargin{:});   % plain S rep, robust for Lap
    LP = @(varargin) DLP(varargin{:}) + SLP(varargin{:});   % test completed rep
    lpnam = 'D+S';
  elseif pde=='H'      % CFIE to avoid gamma interior resonances
    qfsrep = @(varargin) DLP(varargin{:}) + 1i*khelm*SLP(varargin{:});
    LP = @(varargin) DLP(varargin{:}) + 1i*khelm*SLP(varargin{:}); % test CFIE
    lpnam = 'D-ikS';
  elseif pde=='S'      % QFS rep: use S+D (pure S enough when source-free)
    qfsrep = @(varargin) SLP(varargin{:}) + DLP(varargin{:});  % D or no?
    LP = @(varargin) DLP(varargin{:}) + SLP(varargin{:});   % test D+S completed
    lpnam = 'D+S';
  end
  
  srcfacs=1:0.02:1.5;
  chkfacs=1:0.02:1.6;
  A0 = LP(b,b) - 0.5*sgn*eye(ncomp*N);     % Kress A matrix
  k0 = cond(A0)
  kk=nan(numel(srcfacs),numel(chkfacs));        % alloc conv metrics
  for isf = 1:numel(srcfacs)   % sweep over upsampling plane .................
    for icf = 1:numel(chkfacs)
      qfs.srcffac = srcfacs(isf); qfs.chkfac = chkfacs(icf);
      q = qfs_create(b,interior,LP,qfsrep,tol,qfs);
      B = qfsrep(b,q.s);   % bdry from src
      A = (B*q.Q2)*q.Q1;   % QFS nyst approx (BY)Z
      kk(isf,icf) = cond(A);
      fprintf('srcfac=%.3g\tchkfac=%.3g\tkappa(A)=%.6g\n',qfs.srcffac,qfs.chkfac,kk(isf,icf))
    end
  end
  figure; subplot(1,2,1);
  h=imagesc(srcfacs,chkfacs,log10(kk' / k0));
  colormap(goodbw);
  colorbar;
  axis xy tight;
  xlabel('$\upsilon$','interpreter','latex');
  ylabel('$\upsilon_c$','interpreter','latex');
  chkfac0 = 1.5;   % slice val
  hline(chkfac0);
  text(1.25,1.55, '(a) $\log_{10} (\kappa(\tilde A)/\kappa(A))$','interpreter','latex');
  drawnow

%  srcfacs=1:0.05:1.5;
  qfs.chkfac = chkfac0;   % ................. slice along srcfac only
  eA = nan(ncomp*N,numel(srcfacs));
  for isf = 1:numel(srcfacs)   % sweep over upsampling plane .................
    qfs.srcffac = srcfacs(isf);
    q = qfs_create(b,interior,LP,qfsrep,tol,qfs);
    B = qfsrep(b,q.s);   % bdry from src
    A = (B*q.Q2)*q.Q1;   % QFS nyst approx (BY)Z
    eA(:,isf) = eig(A);
    fprintf('srcfac=%.3g\tchkfac=%.3g\tmaxeig(A)=%.6g\n',qfs.srcffac,qfs.chkfac,max(abs(eA(:,isf))))
  end
  subplot(1,2,2);
  plot(srcfacs,abs(eA),'b+');
  xlabel('$\upsilon$','interpreter','latex');
  ylabel('$|\lambda_j(\tilde A)|$','interpreter','latex');
  axis([1 max(srcfacs) 0 3.1]);
  text(1.25,2.8, '(b) spec($\tilde A$) magnitudes, $\upsilon_c=1.5$','interpreter','latex');
end

set(gcf,'paperposition',[0 0 12 3]);
print -dpng tmp.png
file = sprintf('kappaplane.png');
system(['convert tmp.png -trim ' file ' && rm -f tmp.png']);
