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
tol = 1e-6;
N = 100;
b = curve(N);           % setup discretization
srcfacs=1:0.05:1.5;
chkfacs=1:0.05:1.6;

fig=figure;
for ipde=3%1:3                                         % which PDEs to test
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
  elseif pde=='S'      % pure SLP, enough when source-free
    qfsrep = @(varargin) SLP(varargin{:}) + DLP(varargin{:});  % D or no? uh-oh
    LP = @(varargin) DLP(varargin{:}) + SLP(varargin{:});   % test D+S completed
    lpnam = 'D+S';
  end
  
  A0 = LP(b,b) - 0.5*sgn*eye(ncomp*N);
  k0 = cond(A0)
  subplot(1,numel(pdes),ipde);   % one subplot per PDE...
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
  h=imagesc(srcfacs,chkfacs,log10(abs(kk'-k0)));
  colormap(goodbw); caxis([-10 5]);
  colorbar;
  axis xy tight;
  xlabel('$\upsilon$','interpreter','latex');
  ylabel('$\upsilon_c$','interpreter','latex');

  CR = ''; if pde=='H', CR=char(10); end  % insert carriage return
  text(1.0,max(chkfacs), sprintf('%s(%s) %s %s',CR,char(96+ipde),lpnam,pdenam{ipde}));
  drawnow
end

stop

set(gcf,'paperposition',[0 0 12 3]);
print -dpng tmp.png
file = sprintf('kappaplane.png');
system(['convert tmp.png -trim ' file ' && rm -f tmp.png']);
