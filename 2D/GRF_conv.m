function GRF_conv(tol,verb,curvemeth,interior)
% QFS convergence test for GRF, 2D, one smooth curve. Compares to sing quad.
%
% GRF_conv(tol,verb) produces a plot in current figure.
%  tol - requested QFS tol
%  verb - verbosity of text output (fig always produced)
%
% See: GRF_conv_multitol.m for example usage
%
% todo: ext case. Helm (needs new BIE2D kernels)
%
% Barnett 8/15/19

if nargin<1, tol = 1e-10; end          % QFS (don't make too close to emach!)
if nargin<2, verb = 0; end
o = []; o.factor='l';             % QFS opts
o.verb = (verb>1);
if nargin>=3, o.curvemeth=curvemeth; end
if nargin<4, interior = false; end

k = 0;                 % wavenumber (0 for now)
a = .3; w = 5;         % smooth wobbly radial shape params

t0 = 1.042;            % bdry pt to base the rhs pt src (match z0=.2+.7i, Feb'19)
% Here's Newton to solve for b.Z(t)=z0: t = t - (b.Z(t)-z0)./b.Zp(t), repeat!
N0 = 250;              % desired N by which Nystrom (exp rate N) fully converged
imt0 = -log(eps)/N0;              % imag dist of src so Nystrom emach at N=N0
                                  % (or Kress SLP conv, at rate N/2, by N=2*N0)
sgn = -1+2*interior;              % David convention (+1 if interior)
b = wobblycurve(1,a,w,100);       % to access to b.Z
z0 = b.Z(t0 - 1i*sgn*imt0);       % data src pt, given imag dist, sets conv rate
if interior, trg.x = 0.3+0.2i;    % far int target point
else trg.x = 1.5-0.5i; end        % far ext target point
f0 = exp(1i*4.0);      % data strength to give data size O(1)
if k==0                % Laplace
  fholom = @(z) f0./(z-z0);       % holomorphic ext soln, with no log nor const
  fpholom = @(z) -f0./(z-z0).^2;  % its complex deriv
  f = @(z) real(fholom(z));       % use real part as known Laplace soln
  fx = @(z) real(fpholom(z)); fy = @(z) -imag(fpholom(z)); % partials, NB sign!
  SLP = @LapSLP; DLP = @LapDLP;
else    % ... Helm
  SLP = @HelmSLP; DLP = @HelmDLP; % don't yet exist, *** need f too
end
  
srcker = SLP;                     % choose QFS src rep (*** fix when Helm)

Ns = 100:30:600;
es = nan(1,numel(Ns)); cns=nan*Ns;                    % save errors, etc
for i=1:numel(Ns); N=Ns(i);       % ---------- N convergence
  b = wobblycurve(1,a,w,N);
  qs = qfs_create(b,interior,SLP,srcker,tol,o);       % make QFS objects for SLP
  if verb==1, o.verb = (i==numel(Ns)); end            % verb=1: print last only
  qd = qfs_create(b,interior,DLP,srcker,tol,o);       % DLP
  tau = -sgn*f(b.x);                                  % GRF densities (tau=DLP)
  sig = sgn*(fx(b.x).*real(b.nx) + fy(b.x).*imag(b.nx));  % SLP = n-deriv
  co = qs.qfsco(sig) + qd.qfsco(tau);                 % add the QFS coeffs
  cns(i) = sqrt(mean(co.^2));                         % coeff rms
  ub = srcker(b,qs.s,co);          % eval the QFS rep (srces same locs)
  uerr = ub - f(b.x);              % compare to bdry potential
  es(1,i) = sqrt(sum(b.w .* abs(uerr).^2));           % QFS L_2 bdry err
  ubk = SLP(b,b,sig) +  DLP(b,b,tau) - 0.5*sgn*tau;   % Kress sing quad PV + JR
  ukerr = ubk - f(b.x);
  es(2,i) = sqrt(sum(b.w .* abs(ukerr).^2));          % its L_2 bdry err
  es(3,i) = srcker(trg,qs.s,co) - f(trg.x);           % QFS far err
  es(4,i) = LapSLP(trg,b,sig) + LapDLP(trg,b,tau) - f(trg.x);  % native far err
end                               % ----------
%disp([Ns', es'])
%figure(1); clf;  % for a script
semilogy(Ns,abs(es),'+-'); hold on; plot(Ns,cns,'o-');
xlabel('N'); ylabel('error, L^2(\Gamma)');
legend('QFS on-surf','Kress sing quad','far QFS', 'far native','QFS co rms','location','northwest');
title(sprintf('GRF conv, int=%d, tol=%.3g',interior,tol));
axis tight; v=axis; v(3:4)=[1e-16 1e2]; axis(v);                % cement y-range
hold on; plot(Ns,tol+0*Ns,'b:'); text(Ns(1),1.5*tol,'QFS tol'); % horiz line

