% simple Helmholtz exterior scatt BVP by QFS demo, one boundary
% Barnett 2/25/21

clear % BVP params...
a = .3; w = 5;         % smooth wobbly radial shape params
k = 10;
tol = 1e-14;   % QFS
incang = pi/5;
kcomp = k*exp(1i*incang);   % inc wavevector as complex number
ui = @(z) exp(1i*real(conj(kcomp)*z));    % inc plane wave
% *** add uix,uiy, *** try Neu BVP
%  fprintf('inc wave checkgrad: %.3g\n', checkgrad(@(x) ui(x(1,:)+1i*x(2,:)), @(x) [uix(x(1,:)+1i*x(2,:)); uiy(x(1,:)+1i*x(2,:))], [-0.7;0.2]))
% BIE rep...
eta_CFIE = k;             % amount of iS_k; 0 for no CFIE
lpker = @(varargin) HelmDLP(k,varargin{:}) -1i*eta_CFIE*HelmSLP(k,varargin{:});
refker = @(varargin) HelmDLPpotker(k,varargin{:}) -1i*eta_CFIE*HelmSLPpotker(k,varargin{:});  % ref needed for adaptive eval
% QFS internal...
o.onsurf = 1; o.srcfac = 1.0;
qfsker = @(b,varargin) lpker(b,varargin{:}) + o.onsurf*0.5*eye(b.N);  % JR
srcker = @(varargin) HelmDLP(k,varargin{:}) -1i*eta_CFIE*HelmSLP(k,varargin{:});
%srcker = @(varargin) HelmSLP(k,varargin{:});  % plain charge sources

% solution testing...
t.x = [-0.2-1.5i; nan];       % far test point (& reserve as 2x1)
nrdist = 1e-4;                % adaptive DLP dies any closer than 1e-6, sad
b = wobblycurve(1,a,w,100);   % just to setup b.Z etc anal formulae
s=0.7; t.x(2) = b.Z(s) + nrdist * (b.Zp(s)/1i)/abs(b.Zp(s));  % near test pt

Ns=100:20:400;           % N (bdry nodes) convergence study
us=nan(numel(Ns),4);     % cols: far, nr, qfsfar, qfsnr
%profile clear; profile on
for i=1:numel(Ns); N=Ns(i);
  b = wobblycurve(1,a,w,N);
  f = -ui(b.x);           % RHS
  A0 = lpker(b,b) + 0.5*eye(b.N);    % "0" = Kress on-surf: ext JR for Dir data
  dens0 = A0\f;
  us(i,1:2) = lpker(t,b,dens0);     % u_scatt far native Nystrom eval
  [~,densfun] = perispecinterparb(dens0,nan);   % spectral interpolant
  us(i,2) = lpevaladapt(t.x(2), refker, densfun, b, 1e-13);  % adaptive (slow)
  q = qfs_create(b,false,qfsker,srcker,tol,o);         % setup
  co = q.qfsco(dens0);
  us(i,3:4) = lpker(t,q.s,co);       % QFS eval
  fprintf('N=%d:\tfar: Re u nat %.12g \tqfs %.12g\tdiff %.3g\n',N,real(us(i,1)), real(us(i,3)), abs(diff(us(i,[1 3]))))
  fprintf('\tnr:  Re u nat %.12g \tqfs %.12g\tdiff %.3g\n',real(us(i,2)), real(us(i,4)), abs(diff(us(i,[2 4]))))
end
%profile off; profile viewer
figure(1); clf; title('convergence in $n$','interpreter','latex');
ufex = us(end,1); unex = us(end,4);
semilogy(Ns,abs([us(:,1)-ufex, us(:,3)-ufex, us(:,2)-unex, us(:,4)-unex]),'+-');
xlabel('n'); ylabel('u abs error'); legend('far nat','far qfs','nr adap','nr qfs');


g0 = 2.0; ng = 100; u0 = 2.0;  % potential images:  box size, max u value
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
axis xy equal tight; hold on; plot([b.x; b.x(1)],'w.-');
subplot(1,3,3); imagesc(g,g,up+uip); caxis(u0*[-1 1]); title('Re u_{tot}');
axis xy equal tight; hold on; plot([b.x; b.x(1)],'w.-');


