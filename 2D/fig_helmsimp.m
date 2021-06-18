% simple Helmholtz exterior scatt BVP by QFS demo, one boundary
% Barnett 2/25/21

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
  A = (srcker(b,qd.s)*qd.Q2)*qd.Q1;          % BIE matrix - best way? ***
  dens = A\f;                         % solve w/ QFS-D mat
  vs(i,4) = dens(1);               % 1st node @ s=0 always (no dens interp)
  cod = qd.qfsco(dens);
  vs(i,1:2) = lpker(t,b,dens);     % u_scatt far native Nystrom eval
  vs(i,2:3) = srcker(t,qd.s,cod);       % QFS eval, overwrites col 2
  fprintf('\tqfs-d:     nr %.12g (%.3g) \t far %.12g (%.3g)\n',real(vs(i,2)),abs(vs(i,2)-us(i,1)),real(vs(i,3)),abs(vs(i,3)-us(i,3)))
end            % .................................


figure;      % ================ paper fig ===========================
subplot(2,3,3);   % QFS-B conv
ufex = us(end,1); unex = us(end,3);
semilogy(Ns,abs(us(:,1)-ufex),'gs-', Ns,abs(us(:,4)-ufex),'k+-', Ns,abs(us(:,2)-unex),'go-', Ns,abs(us(:,3)-unex),'g*-', Ns,abs(us(:,5)-unex),'k.-');
legend('far, plain', 'far, QFS-B', 'near, plain', 'near, adaptive', 'near, QFS-B');
xlabel('n'); ylabel('abs error in u'); axis([Ns(1) Ns(end-1) 1e-15 1]);
hline(tol,'r:');
text(Ns(2),tol*5,'QFS $\epsilon$','interpreter','latex','color',[1 0 0]);
title('(e)   $n$-node Kress $\tau$ solve; various $u$ eval','interpreter','latex');

subplot(2,3,6);   % QFS-D conv
semilogy(Ns,abs(vs(:,2)-ufex),'k+-', Ns,abs(vs(:,3)-unex),'k.-');
legend('far, QFS-D', 'near, QFS-D');
xlabel('n'); ylabel('abs error in u'); axis([Ns(1) Ns(end-1) 1e-15 1]);
hline(tol,'r:');
text(Ns(2),tol*5,'QFS $\epsilon$','interpreter','latex','color',[1 0 0]);
title('(f)   $n$-node QFS-D $\tau$ solve and $u$ eval','interpreter','latex');

subplot(2,3,1);   % soln image
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
text(0.5,0.5,1,'$\partial\Omega$','interpreter','latex','fontsize',20);
axis xy equal tight; v=axis;
%

vz = [0 .9 0.5 1.4];  % zoomed view box

subplot(2,3,4);  % nystrom discr
Np=70; bp = wobblycurve(1,a,w,Np);  % for plotting only
plot([bp.x; bp.x(1)],'k.'); hold on;
text(real(bp.x(1))+0.1,imag(bp.x(1)),'$\bf{x}_1$','interpreter','latex');
text(real(bp.x(2))+0.1,imag(bp.x(2)),'$\bf{x}_2$','interpreter','latex');
plot(t.x, 'm*'); text(real(t.x(1)),imag(t.x(1)),'far');
text(real(t.x(2)),imag(t.x(2)),'near'); axis equal; axis(v);
title(sprintf('(b) Nystrom discretization $(n=%d)$',Np),'interpreter','latex');
plot(vz([1 1 2 2 1]),vz([3 4 4 3 3]), '-','color',0.8*[1 1 1]);  % zoom box
text(vz(2),vz(4),'zoom in (c), (d)','color',0.8*[1 1 1]);
axis off

subplot(2,3,2);  % QFS-B geom (same tol)
Np=160; bp = wobblycurve(1,a,w,Np);
o.onsurf=1; qp = qfs_create(bp,false,qfsbker,srcker,tol,o);
plot([bp.x; bp.x(1)],'k+'); hold on; plot(qp.s.x,'r.');
axis equal; axis(vz);
title(sprintf('(c) zoom of QFS-B sources $(n=%d)$',Np),'interpreter','latex');

subplot(2,3,5);  % QFS-D geom (same tol)
Np=160; bp = wobblycurve(1,a,w,Np);
o.onsurf=0; qp = qfs_create(bp,false,lpker,srcker,tol,o);
plot([bp.x; bp.x(1)],'k+'); hold on; plot(qp.s.x,'r.');
plot(qp.bf.x,'b.'); plot(qp.c.x,'g.');
axis equal; axis(vz);
title(sprintf('(d) zoom of QFS-D geom $(n=%d)$',Np),'interpreter','latex');

set(gcf,'paperposition',[0 0 12 7]); print -dpng setup.png
system('convert setup.png -trim setup_trim.png')
% now copy to ../../qmf-paper/figs/
% then import to xfig setup_lab, export as EPS, finally:
% convert -density 95.4 setup_lab.eps setup_lab.png
% (here the resolution matches the size of the original png)
% ================================================================







if 0 % old research conv plots......
figure; subplot(1,2,1);
ufex = us(end,1); unex = us(end,5); dens0ex = us(end,7);
semilogy(Ns,abs([us(:,1)-ufex, us(:,4)-ufex, us(:,3)-unex, us(:,5)-unex, us(:,6), us(:,7)-dens0ex, us(:,8)]),'+-');
hold on; plot(Ns,[1e-4*exp(-sw*Ns/2); exp(-sw*Ns); exp(-imag(ss(1))*Ns/2)],'--');
xlabel('n'); ylabel('u abs error'); axis([min(Ns) max(Ns) 1e-15 1]);
legend('Kress far-nat','Kress far-qfs-b','Kress nr-adap','Kress nr-qfs-b',...
       'rhs interp err','Kress den(0)','Kress denhatdecay','ext SS n/2 rate',...
       'ext SS n rate','int SS n/2 rate');
title('Kress-solved, conv'); %,'interpreter','latex');
subplot(1,2,2);
semilogy(Ns,abs([vs(:,1)-ufex, vs(:,2)-ufex, vs(:,3)-unex]),'+-');
xlabel('n'); ylabel('u abs error'); axis([min(Ns) max(Ns) 1e-15 1]);
legend('far-nat','far-qfs','nr-qfs');
title(sprintf('k=%g: QFS-D-solved, conv',k));

figure; semilogy(abs([fhat, d0hat])); hold on;
nn=1:max(Ns)/2; plot(nn,exp([-sw*nn/2; -imag(ss(1))*nn/2; -sw*nn]),'--');
legend('fhat','denshat','ext SS n/2','int SS n/2', 'ext SS n'); axis([0 max(nn) -1e-15 10]); title('Fourier decay');
% int SS n/2 controls denshat decay tail (once fhat died first)
end

if 0   % three potentials...
g0 = 2.0; ng = 200; u0 = 2.0;  % potential images:  box size, max u value
g = linspace(-g0,g0,ng); [xx yy]=meshgrid(g,g); zz = xx(:)+1i*yy(:);
figure(2); clf; subplot(1,3,1);
ii = ~b.inside(zz); uip = nan*xx; uip(ii) = real(ui(zz(ii)));   % eval ext pts
imagesc(g,g,uip); caxis(u0*[-1 1]); title('Re u_i'); axis equal tight xy;
hold on; plot([b.x; b.x(1)],'w.-'); plot(t.x, 'k.'); drawnow
tp.x = zz(ii);       % ext targets for plotting
%tic; u = lpker(tp,b,dens0); toc       % eval u_scatt on grid
tic; u = lpker(tp,q.s,co); toc       % eval u_scatt on grid
up = nan*xx; up(ii) = real(u);
subplot(1,3,2); imagesc(g,g,up); caxis(u0*[-1 1]); title('Re u');
axis xy equal tight; hold on; plot([b.x; b.x(1)],'w.-'); plot(q.s.x,'w.');
subplot(1,3,3); imagesc(g,g,up+uip); caxis(u0*[-1 1]); title('Re u_{tot}');
axis xy equal tight; hold on; plot([b.x; b.x(1)],'w.-');
end

