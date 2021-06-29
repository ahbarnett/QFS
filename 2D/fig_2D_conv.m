% 2D ext QFS-B conv plots, 3 PDEs, on single shape, each of S, D, separately.
% This is more primitive that GRF test - could add GRF easily.
% Start with: QFS-B for Lap {S or D}LP, vs gold std eval.

% Barnett 6/29/21
clear; %close all
interior = 0;
a = .3; w = 5;   % smooth wobbly radial shape params
curve = @(N) wobblycurve(1,a,w,N);    % only to access b.Z
khelm = 20.0;    % wavenumber (Helm only)
mu=0.7;          % viscosity (Sto only)

t0 = 0.5;   % bdry param to base the rhs data src on
imt0 = 0.15;   % choose source dist
b = curve(100);    % only to access b.Z
sgn = -1+2*interior;              % David convention (+1 if interior)
z0 = b.Z(t0 - 1i*sgn*imt0);       % data src pt, given imag dist, sets conv rate
if interior, trg.x = -0.1+0.2i;   % far int target point
else trg.x = 1.5-0.5i; end        % far ext target point
nrdist = 1e-4;                % adaptive DLP dies any closer than 1e-6, sad
s=4.0; trg.x(2) = b.Z(s) - sgn*nrdist * (b.Zp(s)/1i)/abs(b.Zp(s));  % nr test pt
trg.x=trg.x(:);

lp = 'D';    % or 'D'.  LP flavor to test
qfs.onsurf = 1;       % 1 makes QFS-B
verb = 1;
pdes = 'L'; %'LHS';     % which PDEs to test, 1 char each
pdenam = {'Laplace', 'Helmholtz k=20', 'Stokes'};
tols = [1e-4 1e-8 1e-12];
Ns = 30:30:420;
fig=figure;

for ipde=1:numel(pdes)
  pde = pdes(ipde); fprintf('PDE=%s: -----------------------\n',pdenam{ipde})
%  qfs.srcffac=1.05; if pde=='S', qfs.srcffac = 1.2; end     % Sto bump up  
  ncomp = 1 + (pde=='S');           % # vector cmpts dep on PDE type
  if pde=='L'                       % wrap the LPs in PDE-indep way
    SLP = @LapSLP; DLP = @LapDLP;
    SLPker = @LapSLPpotker; DLPker = @LapDLPpotker;
  elseif pde=='H'
    
  elseif pde=='S'
    
  end
  subplot(1,numel(pdes),ipde);
  
  % density func wrt t param, analytic w/ known singularity
  densfun = @(t) (0.5+sin(3*t+1)).*real(exp(4i)./(b.Z(t)-z0));
  if pde=='S'
    densfun = @(t) [(0.5+sin(3*t+1)).*real(exp(4i)./(b.Z(t)-z0)); cos(2*t-1).*real(exp(5i)./(b.Z(t)-z0))];
  end
  
  if lp=='S'   % works for any PDE. For Dirichet (value or vel) data match
    LP = SLP;
    ker = SLPker;
    JRterm = 0;
  elseif lp=='D'
    LP = DLP;
    ker = DLPker;
    JRterm = -0.5*sign_from_side(interior);
  else   % *** also test combos, CFIE, this way? (but not GRF?)
  end
  
  dd=nan(numel(Ns),1); u=[dd,dd]; u0=u;          % alloc conv metrics
  for itol=1:numel(tols)   % *** to add in, legend only for itol=1
    tol = tols(itol); fprintf('qfs tol=%.3g:.......\n',tol)
    for i=1:numel(Ns), N=Ns(i); fprintf('\tN=%d:\n',N);
      b = curve(N);           % set up discretization
      dens = densfun(b.t);    % tau vec
      if itol==1              % do the reference (non-QFS) stuff
        dhat=fft(dens(1:N)); dd(i)=abs(dhat(N/2+1)/max(dhat));  % dens_1 decay
        u0(i,:) = cmpak(LP(trg,b,dens),ncomp);          % smooth rule @ all trgs
        interpdensfun = interpfunfromgrid(dens,ncomp);  % scalar or vector
        u0(i,2) = cmpak(lpevaladapt(trg.x(2), ker, interpdensfun, b, 1e-14),ncomp);  % adap for nr targ only
      end
      if qfs.onsurf                 % QFS-B (needs on-surf A matrix)
        LPqfs = @(varargin) LP(varargin{:}) + JRterm*eye(ncomp*N);
      else, LPqfs = LP; end         % QFS-D
      q = qfs_create(b,interior,LPqfs,LP,tol,qfs);
      cod = q.qfsco(dens);                 % get QFS src coeffs (str) for dens
      u(i,:) = cmpak(LP(trg,q.s,cod),ncomp); % QFS eval (all trg) from QFS co
    end
    eu0 = abs(u0-u0(end,:)); eu = abs(u-u0(end,:));    % errors (vs ref u0)
    fdpred = exp(-imt0*Ns/2);
    
    if itol<numel(tols)
      semilogy(Ns,eu(:,1),'k+-', Ns,eu(:,2),'k.-'); hold on;
    else
      hp = semilogy(Ns,eu(:,1),'k+-', Ns,eu(:,2),'k.-', Ns,eu0(:,1),'b+-', Ns,eu0(:,2),'b.-', Ns,dd,'r-', Ns,fdpred,'m--');  % grab handle for legend
      h = legend(hp,'far, QFS-B','nr, QFS-B', 'far, plain','nr, adap','$\hat \tau$ decay','$e^{-\delta_\ast n/2}$');
      set(h,'interpreter','latex');
      xlabel('N');
      axis([min(Ns), max(Ns), 1e-15 1e0])
    end  
    hline(tol,'b:');
  end
end
