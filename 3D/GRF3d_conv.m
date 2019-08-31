function [Nvs es cns] = GRF3d_conv(shape,tol,verb,surfmeth,param,interior)
% QFS convergence test for GRF, 3D, one smooth surface.
%
% [Nvs es cns] = GRF3d_conv(shape,tol,verb,surfmeth,param,interior)
%  produces a plot in current figure, convergence of on-surface GRF.
%
%  shape - 0 (torus), 1 (cruller), 2 (bent torus)
%  tol - requested QFS tol
%  verb - verbosity of text output (fig always produced)
%  surfmeth, param - passed to qfs3d_create
%  interior - (boolean) whether test interior/exterior GRF
%
% Example:
%  GRF3d_conv(0,1e-6,2,'a',[1.0],false)

% Barnett 8/31/19
if nargin<1, shape = 0; end
if nargin<2, tol = 1e-5; end          % QFS (don't make too close to emach!)
if nargin<3, verb = 0; end
o = []; o.factor='s';             % QFS opts
o.verb = (verb>1);
if nargin>=4, o.surfmeth=surfmeth; o.param = param; end
if nargin<6, interior = false; end

a = 1.0; bpar = 0.5;                  % baseline torus shape params
if shape==0, nam='plain torus';
  Nvs = 20:5:35;
elseif shape==1, nam='cruller';
  bpar = cruller(bpar,0.1,5,3);          % replaces bpar
  Nvs = 40:5:60;
else, nam='bent torus';
  bpar = benttorus(bpar,0.3,2);          % replaces bpar
  Nvs = 20:5:40;
end

if interior
  trg.x = [0.9;0.5;0.2];          % far int target point
  z0 = [1.8;-0.7;1.0]; f0=1.0;    % far ext src pt
else
  trg.x = [1.8;-0.7;1.0];         % far ext targ, or [1.6;0;0] nr.
  z0 = [0.9;0.3;0.1]; f0=0.3;     % far int src point
end
f = @(x) f0./sqrt(sum((x-z0).^2,1));               % pt src Laplace soln
gradf = @(x) -f0*(x-z0)./sqrt(sum((x-z0).^2,1)).^3;
fprintf('did we get analytic grad f right? yes: %.3g\n',checkgrad(f,gradf))
srcker = @Lap3dSLPmat;                     % choose QFS src rep
sgn = sign_from_side(interior);
o.factor = 's';                   % 'l' fails, has O(1) errors, since max(d/h)>5

es = nan(3,numel(Nvs)); cns=nan*Nvs;                  % save errors, etc
for i=1:numel(Nvs); N=Nvs(i)*[2 1];       % --------- N convergence
  b = setup_torus_doubleptr(a,bpar,N);                % bdry quadr
  qs = qfs3d_create(b,interior,@Lap3dSLPmat,srcker,tol,o); % make QFS for SLP
  qd = qfs3d_create(b,interior,@Lap3dDLPmat,srcker,tol,o); % make QFS for DLP
  tau = -sgn*f(b.x)';                                 % GRF densities (tau=DLP)
  sig = sgn*sum(gradf(b.x).*b.nx,1)';                 % SLP = n-deriv
  if i==1 & verb>1, figure(2); clf; oo.nofig=1; showsurffunc(b,sig,oo); hold on;
    plot3(trg.x(1),trg.x(2),trg.x(3),'.'); plot3(z0(1),z0(2),z0(3),'r*');
    title('tau, z_0 and trg'); drawnow; end
  co = qs.qfsco(sig) + qd.qfsco(tau);                 % add the QFS coeffs
  cns(i) = sqrt(mean(co.^2));                         % coeff rms
  ub = srcker(b,qs.s) * co;        % eval the QFS rep (srces same locs)
  uerr = ub - f(b.x)';             % compare to bdry potential
  es(1,i) = sqrt(sum(b.w .* sum(uerr.^2,1)));         % QFS L_2 bdry val err
  es(2,i) = srcker(trg,qs.s)*co - f(trg.x);           % QFS far val err
  es(3,i) = Lap3dSLPmat(trg,b)*sig + Lap3dDLPmat(trg,b)*tau - f(trg.x);  % native far val err
  fprintf('  GRF [%3d,%3d]: QFSsurf %.3g, QFSfar %.3g, nativefar %.3g, nrm:%.3g\n',N(1),N(2),es(1,i),es(2,i),es(3,i),cns(i))
end                               % ----------
figure(1); clf;  % for a script
semilogy(Nvs,abs(es),'+-'); hold on; plot(Nvs,cns,'o-');
xlabel('N_v'); ylabel('error');
legend('QFS on-surf','far QFS', 'far native','QFS co rms','location','northwest');
title(sprintf('GRF3dconv, %s, int=%d, tol=%.3g',nam,interior,tol));
axis tight; %v=axis; v(3:4)=[1e-10 1e3]; axis(v);            % cement y-range
%hold on; plot(Nvs,tol+0*Nvs,'b:'); text(Nvs(1),1.5*tol,'QFS tol'); % horiz line

