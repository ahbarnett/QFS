% try Laplace ext multibody BVP in 3D.
% Barnett 8/15/21
%startup(3);
clear
lpker = @(t,s) Lap3dDLPmat(t,s) + Lap3dSLPmat(t,s);       % D+S rep
verb = 1;     % plots in this driver
interior = 0;
sgn = sign_from_side(interior);
qfsker = @Lap3dSLPmat;                     % choose QFS src rep (SLP is robust)
twosided = true;   % true is better for cond#
tol = nan;       % QFS tol, dummy since not used for surfmeth='d'
gmrestol = 1e-9;     % make nan for direct solve
precond = false;

E0 = [.5,1,1.5];     % baseline semiaxes
K=10;   % build cluster of K same shape ellipsoids
dmin = 0.03;    % target dist of each to prev ellipsoids
rng(0);    % seed
[E R t] = ellipsoid_cluster(E0,K,dmin);        % grow a DLA cluster w/ dmin
figure(1); for j=1:K
  colorvec = mod(j*[.42,.1,.29],1);   % some irrationals so colors unique
  show_ellipsoid(E{j},R{j},t{j},colorvec); hold on;
  text(t{j}(1),t{j}(2),t{j}(3),sprintf('%d',j), 'fontsize',20);
end
axis equal tight vis3d; set(gca,'Clipping','off');
title('ellipsoids cluster (pts for viz only)');

Edir = [1;2;3]; Edir = Edir/norm(Edir);  % uniform field: inc dir, col vec
uinc = @(x) 0.3 * Edir'*x;               % make uinc only size 1 at cluster
%uinfty = 0.0; uinc = @(x) uinfty+0*x(1,:);   % const u_infty (not "field")
%dinc=0.1; z0 = Edir*(1+dinc);
%uinc = @(x) 2*dinc./sqrt(sum((x-z0).^2,1));  % ... or pt src @ z0
voltage = rand(K,1)-1/2; voltage(1)=0;        % body voltages for RHS

b = ellipsoid(E0(1),E0(2),E0(3)); Nufac = 4/3;  % ellipsoid (0<u<2pi, -1<v<1)
qfso.surfmeth = 'd'; qfso.param = [1.0,0.1,3,0.1];  % tuned for axes .5:1:1.5

xfar = [1;-1;2]; %[1.5;0.4;-0.3]; % good far point out to K=10
nrdist = 0*1e-6; unr = 2.0; vnr = 0.7;     % pick near targ above (u,v) surf pt
nnr = cross(b.Zu(unr,vnr),b.Zv(unr,vnr)); nnr=nnr/norm(nnr);  % unit nor
xnear = b.Z(unr,vnr) + nrdist*nnr;           % note on body 1
trg.x = [xfar, xnear];
plot3(trg.x(1,:),trg.x(2,:),trg.x(3,:),'k*'); drawnow

Nvs = 16:8:48; %64;   % Nv convergence, per body
Nvs = 40;

% *** to do set up error conv arrays, index w/ c
us = nan(size(trg.x,2),numel(Nvs));
disp('N-convergence study...');
for c=1:numel(Nvs);       % --------------------------------------------
  N = 2*ceil(Nvs(c)*[Nufac 1]/2);   % pick Nu and Nv (even)
  b = setupsurfquad(b,N);                  % bdry nodes for 1 ellipsoid
  n = b.N;     % n per body

  tim=tic; % QFS for 1 body Nystrom A0, and QFS eval setup...
  qfso.verb = 2;
  qfso.factor='l';   % ridiculously fast (instead filling dominates)
  q = qfs3d_create(b,interior,lpker,qfsker,tol,qfso);
  B = qfsker(b,q.s);  % bdry-from-src          % could use to fill Aij blks
  SFD = q.qfsco(eye(n));   % QFS src-from-dens -  avoids Q1, Q2 if stable?
  if isfield(q,'Q2'), A1 = (B*q.Q2)*q.Q1;     % Nystrom matrix
  else, A1 = B*SFD; end         % bkw-stab (?ok) way for LU
  if twosided
    qfso.I = q.I;        % reuse same I
    q2 = qfs3d_create(b,~interior,lpker,qfsker,tol,qfso);  % for 2-sided
    qfso = rmfield(qfso,'I'); q = rmfield(q,'I');
    B2 = qfsker(b,q2.s);    % other side I/O
    if isfield(q2,'Q2'), A2 = (B2*q2.Q2)*q2.Q1;     % I/O Nystrom matrix
    else, A2 = B2*q2.qfsco(eye(n)); end
    A0 = -sgn*0.5*eye(n) + (A1+A2)/2;        % 2-sided avg Nyst + JR Id term
    clear q2
  else, A0 = A1;        % 1-sided
  end
  % cond(A0), figure;plot(eig(A0),'+'); axis equal  % slow: 2-sided much better
  fprintf('\tQFS 1-body setup %.3g s... ',toc(tim))
  
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
      if j==i, A(I,J) = A0;
      else
        BIJ = qfsker(struct('x',x{i}), s{j});   % bdry_i from src_j
        A(I,J) = BIJ * SFD;        % SFD = 1-body QFS src-from-dens
        % *** precond version to do: eye on diag, right Q off
      end
    end
  end
  fprintf('\tA(%d)filled %.3g s\n',n*K,toc(tim))
  
  tic; if ~isfinite(gmrestol)
    dens = A\f;  %r = A*dens-f; fprintf('A rel resid nrm %g\n',norm(r)/norm(f))
  else, maxit = 1e3; tic
    [dens,flag,relres,iter] = gmres(A,f,b.N,gmrestol,min(maxit,b.N));
    fprintf('\tA GMRES %.3g sec, %d iters, rel resid nrm %.3g (requested %.3g)\n',toc,iter(2),relres,gmrestol)
  end

  for k=1:K, co{k} = q.qfsco(dens((1:n)+(k-1)*n)); end    % get all qfs src co
  % or is plain matvec SFD.dens_k ok here?
  % (could stack & send to FMM?)

  if verb && c==numel(Nvs), figure(3); clf; oo.nofig=1;
    for k=1:K,showsurffunc(struct('x',x{k},'topo','s'),dens((1:n)+(k-1)*n),oo);
      hold on; end
      caxis([min(dens) max(dens)]); title('dens'); drawnow
    figure(4); clf; oo.nofig=1; 
    for k=1:K,showsurffunc(s{k},co{k},oo); hold on; end
    caxis([min(vertcat(co{:})) max(vertcat(co{:}))]); title('qfs co'); drawnow
  end
  
  u = zeros(size(trg.x,2),1);        % sum QFS pots @ trg... quick if few trgs
  for k=1:K, u = u + qfsker(trg,s{k}) * co{k}; end
  uinctrg = uinc(trg.x)';     % col vec, otherwise matlab auto-bsxfuns! stupid
  utot = uinctrg + u;         % physical soln @ trgs (col vec)
  us(:,c) = utot;
  
  % build unified bdry pointset (fields x,nx,w) for lp plain eval of dens...
  bb.x=xx; bb.nx=0*xx; bb.w = 0*xx(1,:);
  for k=1:K, jj=(1:n)+(k-1)*n; bb.nx(:,jj)=R{k}*b.nx; bb.w(jj)=b.w; end
  ufarplain = lpker(trg,bb) * dens;    % plain Nystrom eval, far only
  fprintf('BVP [%3d,%3d]: n=%d;nf=%d (per body), Ntot=%d\tut farpl,far,nr = %.9g,%.9g,%.9g\n',N(1),N(2),n,q.bf.N,n*K,ufarplain(1)+uinctrg(1),utot(1),utot(2))
  
  % to do: eval QFS u at all surf nodes, multibody
  %ub = B*co;
  %fprintf('\tsurf u err rms %.3g\t\trms dens %.3g\trms co %.3g\n',sqrt(mean((ub-f).^2)),sqrt(mean(dens.^2)),sqrt(mean(co.^2)))
 
  % to do: MFS multibody
%  comfs = B\f; %r = B*co-f; fprintf('mfs resid nrm %g\n',norm(r)) % also check
%  umfs = qfsker(trg,q.s) * comfs;    % MFS pot trg
%  utotmfs = uinctrg + umfs;
%  fprintf('\t\t\t\t\t      MFS ut far,nr = %.9g,%.9g\n',utotmfs(1),utotmfs(2))
end
