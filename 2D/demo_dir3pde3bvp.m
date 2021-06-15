% demo all 3 PDEs x 3 Dirichlet BVPs: known int / known ext / scatt ext: 9 pics
% No convergence, just an acceptable N. Driver for: gen_dirbvpconv.m
% Barnett 6/15/21.
clear
qfs.tol = 1e-12;
%qfs.onsurf = 1;  % 1 makes QFS-B: only seems to affect ed (dens err)
o.Ns = 300;
o.verb = 1;
o.grid = [];      % tells to do grid eval
pdes = 'LHS';     % 1 char each
figure;
for ipde=1:3
  pde = pdes(ipde); fprintf('PDE=%s: ----------\n',pde)
  qfs.srcfac=1.0; if pde=='S', qfs.srcfac = 1.2; end     % Sto bump up  
  for ibvp=1:3
    interior = (ibvp==1);
    known = (ibvp<=2);
    fprintf('interior=%d known=%d: ....... \n',interior,known)
    
    [r g] = gendirbvp_conv(pde,interior,known,qfs,o);  % do QFS
    
    tsubplot(3,3,ipde+3*(ibvp-1));
    colormap(jet(256)); up = nan*g.ii;  % plot vals
    u0 = 2.0; if known, u0=1.0; end     % color scale
    if pde=='L'
      up(g.ii) = g.us + (1-known)*g.ui;   % u or utot (or NaN)
      contourf(g.grid.x,g.grid.y,up,linspace(-u0,u0,20));
    elseif pde=='H'
      up(g.ii) = g.us + (1-known)*g.ui;   % u or utot (or NaN)
      surf(g.grid.x,g.grid.y,real(up),'alphadata',~isnan(up));   % NaN->white
    elseif pde=='S'
      pp = nan*g.ii; pp(g.ii) = g.press + (1-known)*g.presi;  % p, ptot (or NaN)
      contourf(g.grid.x,g.grid.y,pp,linspace(-u0,u0,20)); hold on;
      dxq = 0.1; iiq = find(abs(mod(real(g.x)+dxq/2,dxq)-dxq/2)+abs(mod(imag(g.x)+dxq/2,dxq)-dxq/2)<1e-6);  % subsample the grid for vel arrow plot
      m = numel(g.x);   % # targs, indices in the g.x array. Now pack as C#...
      up = g.us(iiq)+1i*g.us(m+iiq) + (1-known)*(g.ui(iiq)+1i*g.ui(m+iiq));
      quiver(real(g.x(iiq)),imag(g.x(iiq)), real(up),imag(up), 2.0,'k-')
    end
    view(2); shading interp; caxis(u0*[-1 1]); grid off;
    title(sprintf('%s: int=%d known=%d',pde,interior,known),'interpreter','latex');
    hold on; plot([g.q.b.x; g.q.b.x(1)],'k-');
    if known, plot(g.z0,'r*'); end
    text(0,0,1,'$\Omega$','interpreter','latex','fontsize',20);  % NB z=1 3D lift
    text(0.5,0.5,1,'$\partial\Omega$','interpreter','latex','fontsize',20);
    axis xy equal tight;
    axis([min(g.grid.x),max(g.grid.x),min(g.grid.y),max(g.grid.y)]);
  end
end
