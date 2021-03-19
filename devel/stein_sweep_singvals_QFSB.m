% test SLP Stokes interior QFS sweep lowest sing val hitting zero issue,
% when source pts are always 5h away but n changes. Also eigval case.
% Barnett 3/15/21, for QFS w/ David Stein. Only needs: BIE2D

clear
mu = 0.7;  % Stokes viscosity param (>0, irrelevant)
ns = 10:1:300;                   % kill the thing
R=1;     % bdry radius
hscale = 5.4; %log(1/eps)/(2*pi)        % how many h to place away (5.7 for epsmach)

fprintf('last three sigma_j(A) for Laplace or Stokes bdry-from-src matrix:\n')
n0 = max(ns);                     % offset within the ss array
ss = nan(3*n0,numel(ns));         % save all sing vals
for i=1:numel(ns), n=ns(i);       % loop over n, also controls source radius
  b = wobblycurve(R,0,1,n);                       % bdry
  s = wobblycurve(R*(1+(2*pi/n)*hscale),0,1,n);   % src pts
  A = LapSLP(b,s);                % bdry-from-src mat
  ss(1:n,i) = eig(A); %svd(A);
  fprintf('n=%d\t   Lap:',n); fprintf('\t%.2g',ss(n-2:n,i));

  A = StoSLP(b,s,mu);             % bdry-from-src mat
  nv = [real(s.nx);imag(s.nx)];   % col vec of the source normals
  %A(1,:) = A(1,:) + nv';          % ones-matrix to kick away nullity 1
  ss(n0+1:n0+2*n,i) = eig(A); %svd(A);
  fprintf('\t   Sto:'); fprintf('\t%.2g',ss(n0+2*n-2:n0+2*n,i));
  fprintf('\n')
end
%figure; plot(b.x,'.'); hold on; plot(s.x,'.'); axis equal tight; title('geom')
figure;
subplot(2,2,1); semilogy(ns,abs(ss(1:n0,:)),'k.','markersize',10);
xlabel('n'); ylabel('\sigma_j'); title('Lap');
subplot(2,2,2); semilogy(ns,abs(ss(n0+1:end,:)),'k.','markersize',10);
xlabel('n'); ylabel('\sigma_j'); title('Sto');
subplot(2,2,3); plot(ns,ss(1:n0,:),'k.','markersize',10);
xlabel('n'); ylabel('\sigma_j'); title('Lap'); set(gca,'ylim',1e-8*[-1 1]);
subplot(2,2,4); plot(ns,ss(n0+1:end,:),'k.','markersize',10);
xlabel('n'); ylabel('\sigma_j'); title('Sto'); set(gca,'ylim',1e-8*[-1 1]);
set(gcf,'GraphicsSmoothing','on')
%print -dpng stein_sweep_singvals_QFSB.png
%print -dpng stein_sweep_eigvals_QFSB.png

% check all eigvecs w/ small eigvals around +-1e-8 are oscillatory: yes...
%[V,D]=eig(A); figure; imagesc(real(V)); colorbar
