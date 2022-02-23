% simple Helmholtz exterior scatt BVP by QFS demo, one boundary
% Barnett 2/22/22 for resub.

clear % BVP params...
a = .3; w = 5;         % smooth wobbly radial shape params (a=.3,w=5)
k = 20; 
incang = pi/5;
kcomp = k*exp(1i*incang);   % inc wavevector as complex number
ui = @(z) exp(1i*real(conj(kcomp)*z));    % inc plane wave
% *** add uix,uiy, *** try Neu BVP
%  fprintf('inc wave checkgrad: %.3g\n', checkgrad(@(x) ui(x(1,:)+1i*x(2,:)), @(x) [uix(x(1,:)+1i*x(2,:)); uiy(x(1,:)+1i*x(2,:))], [-0.7;0.2]))
% BIE rep...
eta_CFIE = k;             % amount of iS_k; 0 for no CFIE
lpker = @(varargin) HelmDLP(k,varargin{:}) -1i*eta_CFIE*HelmSLP(k,varargin{:});
refker = @(varargin) HelmDLPpotker(k,varargin{:}) -1i*eta_CFIE*HelmSLPpotker(k,varargin{:});  % ref needed for adaptive eval
% QFS params...
tol = 1e-12;
o.verb = 1; o.curvemeth='2'; %o.srcffac = 1.05;  % makes little diff
o.factor = 's';        % 'l' is possible too but need to change form A below.
qfsbker = @(b,varargin) lpker(b,varargin{:}) + 0.5*eye(b.N);  % JR, QFS-B only
srcker = @(varargin) HelmDLP(k,varargin{:}) -1i*eta_CFIE*HelmSLP(k,varargin{:});
%srcker = @(varargin) HelmSLP(k,varargin{:});  % plain charge sources

% solution testing...
t.x = [-0.2-2i; nan];       % far test point (& reserve as 2x1)
nrdist = 1e-4;                % adaptive DLP dies any closer than 1e-6, sad
b = wobblycurve(1,a,w,100);   % just to setup b.Z etc anal formulae
s=4.0; t.x(2) = b.Z(s) + nrdist * (b.Zp(s)/1i)/abs(b.Zp(s));  % near test pt

% domain's theoretical exp conv rates...
ss = newtzero(b.Zp,0);        % Schwartz singularities of domain (sqrt type)
ss = ss(find(real(ss)<=1 & real(ss)>=0 & imag(ss)<0.3)); disp(ss)
sw = min(abs(imag(ss)));      % param-plane strip half-width (alpha)

Ns=60:20:440;           % N (bdry nodes) convergence study
us=nan(numel(Ns),8);  % dens via Kress, cols:
% 1 nat far, 2 nat nr, 3 adap-nr, 4 qfs-b far, 5 qfs-b nr, 6 rhsinterr, 7 dens(0), 8 denshatdecay
vs=nan(numel(Ns),4);  % dens via QFS-D, cols:
% 1 nat far, 2 qfs-d far, 3 qfs-d nr, 4 dens(0)
for i=1:numel(Ns); N=Ns(i);   % ..................
  b = wobblycurve(1,a,w,N);
  f = -ui(b.x);           % RHS
  fs = perispecinterparb(f,s); us(i,5)=fs+ui(b.Z(s));   % plain RHS interp err
  fhat=fft(f); us(i,5) = abs(fhat(N/2+1)/max(fhat));   % RHS Fou n/2 decay
  % Kress for dens...
  A0 = lpker(b,b) + 0.5*eye(b.N);    % "0" = Kress on-surf: ext JR for Dir data
  dens0 = A0\f;
  us(i,7) = dens0(1);               % 1st node @ s=0 always (no dens interp)
  d0hat=fft(dens0); us(i,8) = abs(d0hat(N/2+1)/max(d0hat));  % n/2 Fou decay
  us(i,1:2) = lpker(t,b,dens0);     % u_scatt far native Nystrom eval
  [~,densfun] = perispecinterparb(dens0,nan);   % spectral interpolant
  us(i,3) = lpevaladapt(t.x(2), refker, densfun, b, 1e-12);  % adaptive (slow)
  % QFS-B...
  o.onsurf=1; q = qfs_create(b,false,qfsbker,srcker,tol,o);
  if 1, figure(1); clf; plot([b.x; b.x(1)],'k.-'); hold on; plot(q.s.x,'r.');
    plot(b.Z(ss),'*'); axis xy equal tight; title(sprintf('N=%d',N)); drawnow;
  end
  co = q.qfsco(dens0);
  us(i,4:5) = srcker(t,q.s,co);       % QFS eval
  fprintf('N=%d:\tfar: Re u nat %.12g \tqfs-b %.12g\tdiff %.3g\n',N,real(us(i,1)), real(us(i,4)), abs(us(i,1)-us(i,4)))
  fprintf('\tnr:  Re u ada %.12g \tqfs-b %.12g\tdiff %.3g\n',real(us(i,3)), real(us(i,5)), abs(us(i,3)-us(i,5)))
  % QFS-D...
  o.onsurf=0; qd = qfs_create(b,false,lpker,srcker,tol,o);
  A = (srcker(b,qd.s)*qd.Q2)*qd.Q1;          % BIE matrix - stable for SVD fact
  %A = srcker(b,qd.s)*(qd.Q2*qd.Q1);  %  Atilde=BX stops at 1e-10 for tol=1e-14
  %A = (srcker(b,qd.s)/qd.U)*qd.iLPcfb;          % stable BIE matrix for LU fact
  dens = A\f;                         % solve w/ QFS-D mat
  vs(i,4) = dens(1);               % 1st node @ s=0 always (no dens interp)
  cod = qd.qfsco(dens);
  vs(i,1:2) = lpker(t,b,dens);     % u_scatt far native Nystrom eval
  vs(i,2:3) = srcker(t,qd.s,cod);       % QFS eval, overwrites col 2
  fprintf('\tqfs-d:     nr %.12g (%.3g) \t far %.12g (%.3g)\n',real(vs(i,2)),abs(vs(i,2)-us(i,1)),real(vs(i,3)),abs(vs(i,3)-us(i,3)))
end            % .................................




figure;      % ================ paper fig =========================== (resub)
subplot(1,3,1);   % soln image
g0 = 2.0; ng = 200; u0 = 2.0;  % potential images:  box size, max u value
g = linspace(-g0,g0,ng); [xx yy]=meshgrid(g,g); zz = xx(:)+1i*yy(:);
ii = ~b.inside(zz); tp.x = zz(ii);       % ext targets for plotting
%tic; u = lpker(tp,b,dens0); toc       % eval u_scatt on grid (native)
tic; u = lpker(tp,q.s,co); toc       % eval u_scatt on grid (via qfs-b)
up = nan*xx; up(ii) = real(u + ui(zz(ii)));   % eval u_inc and tot
surf(g,g,up,'alphadata',~isnan(up)); view(2); shading interp;   % NaN->white
caxis(u0*[-1 1]); grid off;
title('(a) Re $(u_{inc} + u)$, acoustic pressure','interpreter','latex');
hold on; plot([b.x; b.x(1)],'k-');
text(0,0,1,'$\Omega$','interpreter','latex','fontsize',20);  % NB z=1 3D lift
text(1,-.1,1,'$\partial\Omega$','interpreter','latex','fontsize',20);
axis xy equal tight; v=axis;
plot3(real(t.x),imag(t.x),1+0*t.x, 'k*'); text(real(t.x(1))+0.05,imag(t.x(1)),1,'far');
text(real(t.x(2))+0.05,imag(t.x(2)),1,'near');
%
vz = [0 .9 0.5 1.4];  % zoomed view box
plot3(vz([1 1 2 2 1]),vz([3 4 4 3 3]),ones(1,5), '-','color',0.6*[1 1 1]);  % zoom box
text(vz(2),vz(4),1,'zoom in (b), (c)','color',0.6*[1 1 1]);

subplot(1,3,2);  % QFS-B geom (same tol)
Np=160; bp = wobblycurve(1,a,w,Np);
o.onsurf=1; qp = qfs_create(bp,false,qfsbker,srcker,tol,o);
plot([bp.x; bp.x(1)],'k+'); hold on; plot(qp.s.x,'r.');
axis equal; axis(vz);
title(sprintf('(b) zoom of QFS-B sources $(N=%d)$',Np),'interpreter','latex');

subplot(1,3,3);  % QFS-D geom (same tol)
Np=160; bp = wobblycurve(1,a,w,Np);
o.onsurf=0; qp = qfs_create(bp,false,lpker,srcker,tol,o);
plot([bp.x; bp.x(1)],'k+'); hold on; plot(qp.s.x,'r.');
plot(qp.bf.x,'b.'); plot(qp.c.x,'g.');
axis equal; axis(vz);
title(sprintf('(c) zoom of QFS-D geom $(N=%d)$',Np),'interpreter','latex');

%stop
% Output
set(gcf,'paperposition',[0 0 12 7]); print -dpng setup2.png
system('convert setup2.png -trim setup2_trim.png')
system('cp setup2_trim.png ../../qmf_paper/figs/');
% then import to xfig setup2_lab, export as PDF.
%
% Or, used to do: export to EPS, finally:
% convert -density 95.4 setup_lab.eps setup_lab.png
% (here the resolution matches the size of the original png)
% ================================================================


