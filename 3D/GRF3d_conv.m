function [Nvs es cns] = GRF3d_conv(shape,tol,verb,interior,o)
% QFS convergence test for GRF, 3D, one smooth surface.
%
% [Nvs es cns] = GRF3d_conv(shape,tol,verb,interior,o)
%  produces a plot in current figure, convergence of on-surface GRF.
%
%  shape - 0 (torus), 1 (cruller), 2 (bent torus)
%  tol - requested QFS tol
%  verb - verbosity of text output (fig always produced)
%  interior - (boolean) whether test interior/exterior GRF
%  o - opts struct, passed to qfs3d_create (see its opts)
%
% Examples:
%  o.surfmeth='a'; o.param=1.0; o.factor='l'; GRF3d_conv(3,1e-9,2,false,o)
%  o.surfmeth='d'; o.param=[1,0.12,1.5,0.27]; o.factor='b'; GRF3d_conv(4,1e-4,2,false,o)
% o.surfmeth='d'; o.param=[1,0.1,1.5,0.2]; o.factor='b'; GRF3d_conv(0,1e-5,2,false,o)

% Barnett 8/31/19, spheres 9/6/19, debugged spheres 12/11/19
if nargin<1, shape = 0; end
if nargin<2, tol = 1e-5; end          % QFS (don't make too close to emach!)
if nargin<3, verb = 0; end
if nargin<4, interior = false; end
if nargin<5, o = []; end
if ~isfield(o,'surfmeth'), o.surfmeth='a'; end   % QFS opts...
if ~isfield(o,'factor'), o.factor='s';  end      % default SVD, slow
% NB 'l' fails (has O(1) errors) if max(d/h)>5
if ~isfield(o,'param'), o.param=[1,0.2,2.0,0.2]; end   % srcfac,ds,bfac,dc

trg.x = [0.9;0.5;0.2];          % far int point (for tori)
z0 = [1.8;-0.7;0.5]; f0=1.0;    % far ext pt

qo = [];
if shape==0, nam='plain torus';
  b = modulatedtorus(1.0,0.5);
  Nvs = 24:6:48;
elseif shape==1, nam='cruller';
  b = modulatedtorus(1.0,cruller(0.5,0.1,5,3));
  Nvs = 40:5:60;
elseif shape==2, nam='bent torus';
  b = modulatedtorus(1.0,benttorus(0.5,0.3,2));
  Nvs = 20:5:40;
elseif shape==3, nam='sphere';
  b = ellipsoid(1,1,1);
  trg.x = [0.3;-0.4;0.2];        % inside sphere
  %trg.x = [0;0;0.5];        % inside sphere, axial
  Nvs = 16:8:40;
elseif shape==4, nam='ellipsoid (aspect 2.5)';
  b = ellipsoid(.8,1.3,2);
  %qo.minunodes = 20;   % ??
  trg.x = [0.3;-0.4;0.6];        % inside ellipsoid
  Nvs = 40:8:64;
end
if ~interior, f0 = 0.2; [z0 trg.x] = deal(trg.x,z0); end    % swap src<->trg

if 1         % 1: default: all m-modes excited.
  f = @(x) f0./sqrt(sum((x-z0).^2,1));               % 1 pt src Laplace soln
  gradf = @(x) -f0*(x-z0)./sqrt(sum((x-z0).^2,1)).^3;
else   % (debug) 2 equal pt sources, to kill the m=1 (in fact, all m odd) terms
  z1 = [-z0(1);-z0(2);z0(3)];          % rot pi about z axis
  f = @(x) f0*(1./sqrt(sum((x-z0).^2,1)) + 1./sqrt(sum((x-z1).^2,1)));
  gradf = @(x) -f0*((x-z0)./sqrt(sum((x-z0).^2,1)).^3 + (x-z1)./sqrt(sum((x-z1).^2,1)).^3);
end
if verb>1, fprintf('grad f formula good? yes: %.3g\n',checkgrad(f,gradf)); end
srcker = @Lap3dSLPmat;                     % choose QFS src rep
sgn = sign_from_side(interior);
es = nan(3,numel(Nvs)); cns=nan*Nvs;                  % save errors, etc
for i=1:numel(Nvs); N=Nvs(i)*[2 1];       % --------- N convergence
  b = setupsurfquad(b,N,qo);                          % bdry nodes
  %b.w = b.w.*(1 + 1e-9*randn(size(b.w)));   % test perturbations!
  o.verb = 1;
  q = qfs3d_create(b,interior,{@Lap3dSLPmat,@Lap3dDLPmat},srcker,tol,o);
  qs = q{1}; qd = q{2};                               % get QFS for {S,D}LP
  tau = -sgn*f(b.x)';                                 % GRF densities (tau=DLP)
  sig = sgn*sum(gradf(b.x).*b.nx,1)';                 % SLP = n-deriv
  if verb>1, figure(2); clf; oo.nofig=1;
    showsurffunc(b,sig,oo); hold on; plot3(trg.x(1),trg.x(2),trg.x(3),'.');
    plot3(z0(1),z0(2),z0(3),'r*'); title('GRF tau, z_0 and trg'); drawnow; end
  co = qs.qfsco(sig) + qd.qfsco(tau);                 % add the QFS coeffs
  cns(i) = sqrt(mean(abs(co).^2));                    % coeff rms
  ub = srcker(b,qs.s) * co;        % eval the QFS rep (srces same locs)
  uerr = ub - f(b.x)';             % compare to bdry potential
  if verb>1, figure(3); clf; oo.nofig=1; showsurffunc(b,uerr,oo); title('GRF u err');
    figure(4); clf; oo.nofig=1; showsurffunc(qs.s,co,oo); title('QFS coeffs');
    drawnow; end
  %es(1,i) = sqrt(sum(b.w .* sum(uerr.^2)));           % QFS L_2 bdry val err
  es(1,i) = max(abs(uerr));                    % QFS max bdry val err
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

