% 2D QFS-B conv plots, 3 PDEs.
% Barnett 6/22/21
clear
interior = 0;
known = 1;       
qfs.onsurf = 1;  % 1 makes QFS-B
o.verb = 1;
pdes = 'LHS';     % 1 char each
pdenam = {'Laplace', 'Helmholtz k=20', 'Stokes'};
tols = [1e-4 1e-8 1e-12];
N = 30:30:420;
o.Ns = N;
o.imt0 = 0.15;
fig=figure;
for ipde=1:3
  pde = pdes(ipde); fprintf('PDE=%s: ----------\n',pde)
%  qfs.srcffac=1.05; if pde=='S', qfs.srcffac = 1.2; end     % Sto bump up  
  subplot(1,3,ipde);
  for itol=1:numel(tols)
    qfs.tol = tols(itol); fprintf('tol=%.3g:.......\n',qfs.tol)
    [r g] = gendirbvp_conv(pde,interior,known,qfs,o);  % do QFS
    figure(fig);
    if itol==1
      fdpred = exp(-o.imt0*N/2);
      semilogy(N,r.eu(:,1),'g+-', N,r.eu(:,2),'g.-', N,r.d0d,'r-', N,fdpred,'m--');
      hold on;
    end
    semilogy(N,r.eu(:,3),'k+-', N,r.eu(:,4),'k.-');   % QFS-B   
    axis([min(N), max(N), 1e-15 1e0])
    xlabel('n');
    hline(qfs.tol,'b:');
    text(min(N)+10,qfs.tol*0.3,sprintf('$\\epsilon=$%.0e',qfs.tol),'color',[0 0 1],'interpreter','latex');
%    plot([300 max(N)],qfs.tol*[1 1],'k:');   % partial hline
  end
  h=legend('far, plain','nr, adap','$\hat \tau$ decay',...
         '$e^{-\delta_\ast n/2}$', 'far, QFS-B','nr, QFS-B');
  set(h,'interpreter','latex');
  % *** to do add eps labels?  
  text(min(N)+20, 0.3, sprintf('(%s) %s',char(96+ipde),pdenam{ipde}));
  drawnow
end
set(gcf,'paperposition',[0 0 12 4]);
print -dpng tmp.png
system('convert tmp.png -trim 2DB_conv.png && rm -f tmp.png');
