% 2D ext QFS-B/D conv plots, 3 PDEs, on single shape, each of S, D, separately.
% This is more primitive than GRF test (could add GRF easily, no need yet).
% Compares to gold std eval and LSC2D meth.
% Barnett 6/29/21

clear; %close all
interior = 0;
a = .3; w = 5;   % smooth wobbly radial shape params
curve = @(N) wobblycurve(1,a,w,N);    % wrap the discretized bdry struct-maker
khelm = 20.0;    % wavenumber (Helm only)
mu=0.7;          % viscosity (Sto only)

t0 = 0.5;   % bdry param to base the rhs data src on
imt0 = 0.15;   % choose source dist
sgn = -1+2*interior;              % David convention (+1 if interior)
tsing = t0 - 1i*sgn*imt0;         % singularity in complex t
Nwhatever = 100; b = curve(Nwhatever); clear Nwhatever   % only to access b.Z
z0 = b.Z(tsing);                  % data src pt, given imag dist, sets conv rate
if interior, trg.x = -0.1+0.2i;   % far int target point
else trg.x = 1.5-0.5i; end        % far ext target point
nrdist = 1e-4;                % adaptive DLP dies any closer than 1e-5, sad
s=4.0; trg.x(2) = b.Z(s) - sgn*nrdist * (b.Zp(s)/1i)/abs(b.Zp(s));  % nr test pt
trg.x=trg.x(:);

lp = 'D';             % or 'D'.  LP flavor to test (changes output filename)
qfs.onsurf = 1;       % 1 makes QFS-B, 0 for QFS-D
qfs.factor = 's'; qfs.meth='2'; qfs.verb=1; % QFS meth pars
pdes = 'LHS';         % which PDEs to test, 1 char each
pdenam = {'Laplace', sprintf('Helm.%s        k=%g',char(10),khelm), 'Stokes'};
tols = [1e-4 1e-8 1e-12];
Ns = 30:30:420;
fig=figure;

for ipde=1:numel(pdes)
  pde = pdes(ipde); fprintf('PDE=%s: -----------------------\n',pdenam{ipde})
  ncomp = 1 + (pde=='S');           % # vector cmpts dep on PDE type
  side = 'e'; if interior, side = 'i'; end  % LSC2D bary setup
  qfseta = 1.0;                      % D+eta.S mixing for QFS rep
  if pde=='S', qfs.srcffac = 1.2; end     % Sto bump up for QFS-B even
  if pde=='L'                       % wrap the LPs in PDE-indep way
    SLP = @LapSLP; DLP = @LapDLP;
    SLPker = @LapSLPpotker; DLPker = @LapDLPpotker;
    SLPclo = @(t,s,d) LapSLP_closeglobal(t,s,d,side);
    DLPclo = @(t,s,d) LapDLP_closeglobal(t,s,d,side);
  elseif pde=='H'
    SLP = @(varargin) HelmSLP(khelm,varargin{:});
    DLP = @(varargin) HelmDLP(khelm,varargin{:});
    SLPker = @(varargin) HelmSLPpotker(khelm,varargin{:});
    DLPker = @(varargin) HelmDLPpotker(khelm,varargin{:});
    qfseta = 1i*khelm;              % use Kress CFIE for QFS's src rep
  elseif pde=='S'
    SLP = @(t,s,varargin) StoSLP(t,s,mu,varargin{:});
    DLP = @(t,s,varargin) StoDLP(t,s,mu,varargin{:});
    SLPker = @(varargin) StoSLPvelker(mu,varargin{:});
    DLPker = @(varargin) StoDLPvelker(mu,varargin{:});
    SLPclo = @(t,s,d) StoSLP_closeglobal(t,s,mu,d,side);
    DLPclo = @(t,s,d) StoDLP_closeglobal(t,s,mu,d,side);
  end
  
  % density func wrt t param, real (for L,S) analytic w/ known singularity...
  densfun = @(t) (0.5+sin(3*t+1)).*cot((t-tsing)/2);  % cmplx, w/ peri t-sing
  if pde=='L', densfun = @(t) real(densfun(t)); end
  if pde=='S'
    densfun = @(t) [(0.5+sin(3*t+1)).*real(exp(4i)*cot((t-tsing)/2)); cos(2*t-1).*real(exp(5i)*cot((t-tsing)/2))];
  end
  
  if lp=='S'   % works for any PDE. For Dirichet (pot or vel) data match
    LP = SLP;
    ker = SLPker;
    LPclo = SLPclo;
    JRterm = 0;     % Id term for on-surf A fill (pot or vel, QFS needs)
  elseif lp=='D'
    LP = DLP;
    ker = DLPker;
    LPclo = DLPclo;
    JRterm = -0.5*sgn;
  else   % *** also test combos, CFIE, this way? (but not GRF?)
  end
  
  subplot(1,numel(pdes),ipde);   % one subplot per PDE...
  dd=nan(numel(Ns),1); u=[dd,dd]; u0=u; uc=u;          % alloc conv metrics
  for itol=1:numel(tols)
    tol = tols(itol); fprintf('qfs tol=%.3g:.......\n',tol)
    for i=1:numel(Ns), N=Ns(i); %fprintf('\tN=%d:\n',N);
      b = curve(N);           % setup discretization
      dens = densfun(b.t);    % tau vec is our input
      if itol==1              % do the reference (non-QFS) stuff
        dhat=fft(dens(1:N)); dd(i)=abs(dhat(N/2+1)/max(dhat));  % dens_1 decay
        u0(i,:) = cmpak(LP(trg,b,dens),ncomp);          % smooth rule @ all trgs
        interpdensfun = interpfunfromgrid(dens,ncomp);  % scalar or vector
        u0(i,2) = cmpak(lpevaladapt(trg.x(2), ker, interpdensfun, b, 1e-12),ncomp);  % adap for nr targ only
        uc(i,:) = cmpak(LPclo(trg,b,dens),ncomp);       % LSC2D bary close eval
      end
      % QFS setup & do...
      if qfs.onsurf, qfsnam = 'QFS-B';          % needs on-surf A matrix
        LPA = @(varargin) LP(varargin{:}) + JRterm*eye(ncomp*N);
      else, LPA = LP; qfsnam = 'QFS-D'; end     % off-surf A, aka E, matrix
      qfsrep = @(varargin) DLP(varargin{:}) + qfseta*SLP(varargin{:});
      q = qfs_create(b,interior,LPA,qfsrep,tol,qfs);
      cod = q.qfsco(dens);                 % get QFS src coeffs (str) for dens
      u(i,:) = cmpak(qfsrep(trg,q.s,cod),ncomp);   % QFS: eval (all trg), done
    end
    eu0 = abs(u0-u0(end,:)); eu = abs(u-u0(end,:));    % errors (vs ref u0)
    euc = abs(uc-u0(end,:));
    
    if itol<numel(tols)
      semilogy(Ns,eu(:,1),'k+-', Ns,eu(:,2),'k.-'); hold on;
    else
      fdpred = exp(-imt0*Ns/2);             % Nyq Fou decay prediction
      if pde=='H'             % no bary for Helm, or would be new research :(
        hp = semilogy(Ns,eu(:,1),'k+-', Ns,eu(:,2),'k.-', Ns,eu0(:,1),'g+-', Ns,eu0(:,2),'g.-', Ns,dd,'r-', Ns,fdpred,'m--');  % grab handle for legend
        h = legend(hp,['far, ',qfsnam],['nr, ',qfsnam], 'far, plain','nr, adap.','$\hat \tau$ decay','$e^{-\delta_\ast N/2}$');
      else
        hp = semilogy(Ns,eu(:,1),'k+-', Ns,eu(:,2),'k.-', Ns,eu0(:,1),'g+-', Ns,eu0(:,2),'g.-', Ns,euc(:,2),'co-', Ns,dd,'r-', Ns,fdpred,'m--');  % grab handle for legend
        h = legend(hp,['far, ',qfsnam],['nr, ',qfsnam], 'far, plain','nr, adap.','nr, bary.','$\hat \tau$ decay','$e^{-\delta_\ast N/2}$');
      end
      set(h,'interpreter','latex');
      xlabel('N'); ylabel('abs error');
      axis([Ns(1), Ns(end-1), 1e-15 1e0])
    end  
    hline(tol,'b:');
    text(Ns(1)+10,tol*0.3,sprintf('$\\epsilon_a=$%.0e',tol),'color',[0 0 1],'interpreter','latex');
  end
  text(Ns(1)+20, 0.3, sprintf('(%s) %sLP %s',char(96+ipde),lp,pdenam{ipde}));
  drawnow
end

set(gcf,'paperposition',[0 0 12 4]);
print -dpng tmp.png
file = sprintf('2DB_conv_%sLP.png',lp);
system(['convert tmp.png -trim ' file ' && rm -f tmp.png']);
