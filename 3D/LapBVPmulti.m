% try Laplace ext multibody BVP in 3D, ellipsoids only.
% Barnett 8/15/21
%startup(3);
clear
lpker = @(t,s) Lap3dDLPmat(t,s) + Lap3dSLPmat(t,s);       % D+S rep
verb = 0;     % plots in this driver
interior = 0;
sgn = sign_from_side(interior);
qfsker = @Lap3dSLPmat;                     % choose QFS src rep (SLP is robust)
twosided = true;   % true is better for cond#
tol = nan;       % QFS tol, dummy since not used for surfmeth='d'
gmrestol = 1e-8;     % make nan for direct solve
precond = false;

E0 = [.5,1,1.5];     % baseline semiaxes
K=10;   % build cluster of K same shape ellipsoids
dmin = 0.1; %0.03;    % target dist of each to prev ellipsoids
rng(0);    % seed
[E R t xclo] = ellipsoid_cluster(E0,K,dmin);        % grow a DLA cluster w/ dmin
if verb>1, figure(1); for j=1:K
    colorvec = mod(j*[.42,.1,.29],1);   % some irrationals so colors unique
    show_ellipsoid(E{j},R{j},t{j},colorvec); hold on;
    text(t{j}(1),t{j}(2),t{j}(3),sprintf('%d',j), 'fontsize',20);
  end
  axis equal tight vis3d; set(gca,'Clipping','off');
  title('ellipsoids cluster (pts for viz only)');
end

Edir = [1;2;3]; Edir = Edir/norm(Edir);  % uniform field: inc dir, col vec
uinc = @(x) 0.3 * Edir'*x;               % make uinc only size 1 at cluster
%uinfty = 0.0; uinc = @(x) uinfty+0*x(1,:);   % const u_infty (not "field")
%dinc=0.1; z0 = Edir*(1+dinc);
%uinc = @(x) 2*dinc./sqrt(sum((x-z0).^2,1));  % ... or pt src @ z0
voltage = rand(K,1)-1/2; voltage(1)=0;        % body voltages for RHS
%voltage = 0*voltage;  % zero voltages for an easy prob, no singular densities

b = ellipsoid(E0(1),E0(2),E0(3)); Nufac = 4/3;  % ellipsoid (0<u<2pi, -1<v<1)
%qfso.surfmeth = 'd'; qfso.param = [1.0,0.1,3,0.1];  % tuned for axes .5:1:1.5
qfso.surfmeth = 'd'; qfso.param = [1.0,0.08,3,0.1];  % tuned for axes .5:1:1.5

xfar = [1;-1;2]; %[1.5;0.4;-0.3]; % good far point out to K=10
nrdist = 0*1e-6; unr = 2.0; vnr = 0.7;     % pick near targ above (u,v) surf pt
nnr = cross(b.Zu(unr,vnr),b.Zv(unr,vnr)); nnr=nnr/norm(nnr);  % unit nor
xnear = b.Z(unr,vnr) + nrdist*nnr;           % note on body 1
trg.x = [xfar, xnear]; trg.N = 2;
if verb>1, plot3(trg.x(1,:),trg.x(2,:),trg.x(3,:),'k*'); drawnow; end

zsli = xclo{2}(3);    % slice through a near-touching pt, make a pointset obj
dx=0.025; g=-2.6:dx:3.8; ng=numel(g); [xx yy] = ndgrid(g,g);
sli = []; sli.N = numel(xx);
sli.x(1,:) = xx(:)'; sli.x(2,:) = yy(:)'; sli.x(3,:) = zsli*ones(1,sli.N);
slimask = ones(1,sli.N);
for k=1:K, xt = R{k}\(sli.x - t{k});  % sli trgs transformed to ellip coord sys
  slimask((xt(1,:)/E0(1)).^2+(xt(2,:)/E0(2)).^2+(xt(3,:)/E0(3)).^2<=1)=nan;
end
%figure; imagesc(g,g,reshape(slimask,ng,ng)'); axis xy equal tight; % check

Nvs = 16:8:64;   % Nv convergence, per body
%Nvs = 64;

nc = numel(Nvs);   % # of resolutions. things to save during convergence...
us = nan(trg.N,nc); ks = nan(1,nc);  % u at trg, cond #s
uslis = nan(sli.N,nc);  % slice vals
disp('N-convergence study...');
for c=1:nc;                % --------------------------------------------
  N = 2*ceil(Nvs(c)*[Nufac 1]/2);   % pick Nu and Nv (even)
  b = setupsurfquad(b,N);                  % bdry nodes for 1 ellipsoid
  n = b.N;     % n per body

  tim=tic; % QFS for 1 body Nystrom A0, and QFS eval setup...
  qfso.verb = 2;
  qfso.factor='l';   % ridiculously fast (instead filling dominates)
  q = qfs3d_create(b,interior,lpker,qfsker,tol,qfso);
  B = qfsker(b,q.s);  % bdry-from-src          % could use to fill Aij blks
  if isfield(q,'Q2'), A1 = (B*q.Q2)*q.Q1;     % Nystrom matrix (SVD way)
  else, A1 = B*q.X; end         % not bkw-stab way, for LU etc (q.X=SFB)
  if twosided
    qfso.I = q.I;        % reuse same I
    q2 = qfs3d_create(b,~interior,lpker,qfsker,tol,qfso);  % for 2-sided
    qfso = rmfield(qfso,'I'); q = rmfield(q,'I');
    B2 = qfsker(b,q2.s);    % other side I/O
    if isfield(q2,'Q2'), A2 = (B2*q2.Q2)*q2.Q1;     % I/O Nystrom matrix
    else, A2 = B2*q2.X; end
    A0 = -sgn*0.5*eye(n) + (A1+A2)/2;        % 2-sided avg Nyst + JR Id term
    clear q2 A1 A2 B2
  else, A0 = A1; clear A1       % 1-sided
  end
  fprintf('\tQFS 1-body setup %.3g s... ',toc(tim))
  tic; ks(c) = cond(A0); fprintf('cond A0 = %g in %.3g s\n',ks(c),toc);
  %figure;plot(eig(A0),'+'); axis equal  % slow: 2-sided much better
  
  % RHS...
  for k=1:K, x{k} = R{k}*b.x+t{k}; end   % set up cell arrays of bdry nodes only
  xx = horzcat(x{:});
  f = -uinc(xx)';                       % col vec of all RHS
  for k=1:K, jj=(1:n)+(k-1)*n; f(jj)=f(jj)+voltage(k); end  % add voltages to bodies
  
  if verb && c==numel(Nvs), figure(2); clf; oo.nofig=1;
    for k=1:K,showsurffunc(struct('x',x{k},'topo','s'),f((1:n)+(k-1)*n),oo); hold on; end
    caxis([min(f) max(f)]);
    plot3(trg.x(1,:),trg.x(2,:),trg.x(3,:),'k*');
    if exist('Edir','var'), plot3([0 Edir(1)],[0 Edir(2)],[0 Edir(3)],'m-'); end
    %plot3(q.s.x(1,:),q.s.x(2,:),q.s.x(3,:),'r.'); plot3(q.c.x(1,:),q.c.x(2,:),q.c.x(3,:),'g.');
    if exist('z0','var'), plot3(z0(1),z0(2),z0(3),'r*'); end  % inc src
    title('RHS, trgs, src, chk'); drawnow
  end
  
  % fill dense multibody A... (rather than apply offdiag blks on fly, or FMM)
  tim=tic;
  if precond
    tic; Q = inv(A0); fprintf('1-body precond in %.3g s\n',toc);  % all bodies
  end
  for k=1:K                   % cell array of translated+rot src objs...
    s{k}=q.s; s{k}.x=R{k}*q.s.x+t{k}; s{k}.nx=R{k}*q.s.nx; s{k}.w=q.s.w;
  end
  A = nan(K*n);
  for i=1:K, I=((1:n)+(i-1)*n);  % row index set
    for j=1:K, J=((1:n)+(j-1)*n);  % col index set
      if j==i, A(I,J) = A0;        % same diag blk
      else
        BIJ = qfsker(struct('x',x{i}), s{j});   % bdry_i from src_j
        A(I,J) = BIJ * q.X;    % X = 1-body QFS src-from-dens (unstab way)
        % *** precond version to do: eye on diag, right Q off
      end
    end
  end
  clear BIJ
  fprintf('\tA(%d)filled %.3g s\n',n*K,toc(tim))
  
  tic; if ~isfinite(gmrestol)
    dens = A\f;  %r = A*dens-f; fprintf('A rel resid nrm %g\n',norm(r)/norm(f))
  else, maxit = 1e3; tic
    [dens,flag,relres,iter] = gmres(A,f,b.N,gmrestol,min(maxit,b.N));
    fprintf('\tA GMRES %.3g sec, %d iters, rel resid nrm %.3g (requested %.3g)\n',toc,iter(2),relres,gmrestol)
  end
  clear A   % (10 GB for Nv=56, n=3604, K=10)

  for k=1:K, co{k} = q.qfsco(dens((1:n)+(k-1)*n)); end    % get all qfs src co
  % or is plain matvec SFD.dens_k ok here?
  % could then stack & send to FMM?...

  if verb && c==numel(Nvs), figure(3); clf; oo.nofig=1;
    for k=1:K,showsurffunc(struct('x',x{k},'topo','s'),dens((1:n)+(k-1)*n),oo);
      hold on; end
      caxis([min(dens) max(dens)]);
      title('(b) density $\tau$','interpreter','latex','fontsize',14);
    %figure(4); clf; oo.nofig=1; 
    %for k=1:K,showsurffunc(s{k},co{k},oo); hold on; end
    %caxis([min(vertcat(co{:})) max(vertcat(co{:}))]); title('qfs co'); drawnow
  end
  
  u = zeros(trg.N,1);         % sum QFS pots @ trg... quick if few trgs
  for k=1:K, u = u + qfsker(trg,s{k}) * co{k}; end
  uinctrg = uinc(trg.x)';     % col vec, otherwise matlab auto-bsxfuns! stupid
  utot = uinctrg + u;         % physical soln @ trgs (col vec)
  us(:,c) = utot;
  
  tic; uslis(:,c) = uinc(sli.x)' .* slimask';  % sum uinc + QFS pots @ slice trg
  for k=1:K, uslis(:,c) = uslis(:,c) + qfsker(sli,s{k}) * co{k}; end
  fprintf('\tslice eval %d pts in %.3g s\n',sli.N,toc)
  if verb && c==nc, figure(5); imagesc(g,g,reshape(uslis(:,c),ng,ng)');
    xlabel('x'); ylabel('y'); axis xy equal tight; caxis([-1 1]); title(sprintf('Nv=%d',Nvs(c))); drawnow; end
  
  % build unified bdry pointset (fields x,nx,w) for lp plain eval of dens...
  bb.x=xx; bb.nx=0*xx; bb.w = 0*xx(1,:);
  for k=1:K, jj=(1:n)+(k-1)*n; bb.nx(:,jj)=R{k}*b.nx; bb.w(jj)=b.w; end
  ufarplain = lpker(trg,bb) * dens;    % plain Nystrom eval, far only
  fprintf('BVP [%3d,%3d]: n=%d;nf=%d (per body), Ntot=%d\tut farpl,far,nr = %.9g,%.9g,%.9g\n',N(1),N(2),n,q.bf.N,n*K,ufarplain(1)+uinctrg(1),utot(1),utot(2))
  
  % to do: eval QFS u at all surf nodes, multibody
  %ub = B*co;
  %fprintf('\tsurf u err rms %.3g\t\trms dens %.3g\trms co %.3g\n',sqrt(mean((ub-f).^2)),sqrt(mean(dens.^2)),sqrt(mean(co.^2)))
   % to do: MFS multibody?
end                     % ------------------------------------------------



stop
% random stuff...

figure(6); imagesc(g,g,reshape(uslis(:,end)-uslis(:,end-1),ng,ng)'); xlabel('x'); ylabel('y'); axis xy equal tight; colorbar; title('last sli u diff');

disp('far trg dist from each ellipsoid:')
for k=1:K, dist_ellipsoid(R{k}\(trg.x(:,1)-t{k}),E0), end
% min is 1.15

disp('sli min dist to any body:')       % takes a few secs
sli.Nout = sum(~isnan(slimask));
slidist = 1e3*ones(sli.Nout,1);
for k=1:K, xt = R{k}\(sli.x(:,~isnan(slimask))-t{k});  % transf valid sli pts
  distk = nan(sli.Nout,1);
  for j=1:sli.Nout, distk(j) = dist_ellipsoid(xt(:,j),E0); end  % slow
  slidist = min(slidist,distk);
end
min(slidist)  %  2.19347149836403e-05


% PAPER FIGS...
  
figure(3); clf; oo.nofig=1; oo.siz=1.5; sc=0.3; % dens
for k=1:K,showsurffunc(struct('x',x{k},'topo','s'),sc*dens((1:n)+(k-1)*n),oo);
  hold on; end
  caxis(.5*[-1 1]); %caxis([min(dens) max(dens)]);
title('(b) density $\tau$ on $\partial\Omega$, and slice of $u_{tot}$','interpreter','latex','fontsize',16);
h=surf(g,g,zsli*ones(ng,ng),reshape(uslis(:,c),ng,ng)');
shading interp; view(5,20); hc=colorbar; set(hc,'Position',[.8 .2 0.03 0.7]) 
len=4; El=len*Edir; plot3([0 El(1)],[0 El(2)],[0 El(3)],'k-');
text(El(1)+.2,El(2),El(3),'$\mathbf{E}_{inc}$','interpreter','latex','fontsize',14);
set(gcf,'paperposition',[0 0 8 6]); print -dpng tmp.png
system('convert tmp.png -trim ellipsli.png');

figure(7); clf;
for c=1:nc-1; maxesli(c)=max(abs(uslis(:,c)-uslis(:,end))); end
maxesli(nc)=nan;
eufar=abs(us(1,:)-us(1,end)); eufar(end)=nan;
hp=semilogy(Nvs,eufar,'k+-', Nvs,abs(us(2,:)), 'k.-', Nvs, maxesli,'r.-');
axis([15 65 1e-7 0.03]);
xlabel('$N_v$ (vertical nodes per body)','interpreter','latex');
ylabel('abs error in $u$','interpreter','latex');
h=legend(hp,'far, self', 'on-surface', '$L^\infty$ slice, self');
set(h,'interpreter','latex');
%title('(c) QFS-D convergence');
text(5,0.03,'(c)','fontsize',14);
set(gcf,'paperposition',[0 0 3.3 4]); print -dpng tmp.png
system('convert tmp.png -trim ellipconv.png');
%system('convert tmp.png -trim ellipconv_zerovoltages.png');
