% try Laplace ext BVP in 3D, "scattering" response to applied potential
% Barnett 8/14/21
%startup(3);
clear
lpker = @(t,s) Lap3dDLPmat(t,s) + Lap3dSLPmat(t,s);       % D+S rep
verb = 1;     % plots in this driver
interior = 0;
sgn = sign_from_side(interior);
qfsker = @Lap3dSLPmat;                     % choose QFS src rep (SLP is robust)
twosided = 1;   % 1 is better for cond#
tol = nan; % dummy since not used for surfmeth='d'

Edir = [1;2;3]; Edir = Edir/norm(Edir);  % inc dir, col vec
uinc = @(x) Edir'*x;                     % uniform field
%dinc=0.1; z0 = Edir*(1+dinc);
%uinc = @(x) 2*dinc./sqrt(sum((x-z0).^2,1));  % ... or pt src @ z0

% shape
%b = ellipsoid(1,1,1); Nufac=2; % sphere (0<u<2pi, -1<v<1)
%qfso.surfmeth = 'd'; qfso.param = [1.0,0.2,1.5,0.3];  % upsamp a tiny bit?
%tol=1e-9; qfso.surfmeth = 'a'; qfso.param = [1.0];    % good for sphere

b = ellipsoid(.5,1,1.5); Nufac = 4/3;  % ellipsoid (0<u<2pi, -1<v<1)
qfso.surfmeth = 'd'; qfso.param = [1.0,0.1,3,0.1];
%qfso.surfmeth = 'a'; qfso.param = [1.0];

qo = [];            % quadr options
%qo.minunodes = 16;   % makes E nonsquare always, useful for test

xfar = [1.5;0.4;-0.3];
nrdist = 0*1e-6; unr = 2.0; vnr = 0.7;     % pick near targ above (u,v) surf pt
nnr = cross(b.Zu(unr,vnr),b.Zv(unr,vnr)); nnr=nnr/norm(nnr);  % unit nor
xnear = b.Z(unr,vnr) + nrdist*nnr;
trg.x = [xfar, xnear];

Nvs = 16:8:48; %64;   % Nv convergence

eu = nan(numel(Nvs));
for i=1:numel(Nvs);       % --------- N convergence
  N = 2*ceil(Nvs(i)*[Nufac 1]/2);   % pick Nu and Nv (even)
  b = setupsurfquad(b,N,qo);                          % bdry nodes
  f = -uinc(b.x)';                       % col vec

  qfso.verb = 0;  % 2 for subtimings, err estims, etc
  qfso.factor='l';   % ridiculously fast (fills dominate)
  q = qfs3d_create(b,interior,lpker,qfsker,tol,qfso);
  B = qfsker(b,q.s);  % bdry-from-src          % could use to fill Aij blks
  if isfield(q,'Q2'), A1 = (B*q.Q2)*q.Q1;     % Nystrom matrix
  else, A1 = B*q.qfsco(eye(b.N)); end         % bkw-stab way for LU
  
  if twosided
    qfso.I = q.I;        % reuse same I
    q2 = qfs3d_create(b,~interior,lpker,qfsker,tol,qfso);  % for 2-sided
    qfso = rmfield(qfso,'I');
    B2 = qfsker(b,q2.s);    % other side I/O
    if isfield(q2,'Q2'), A2 = (B2*q2.Q2)*q2.Q1;     % I/O Nystrom matrix
    else, A2 = B2*q2.qfsco(eye(b.N)); end
    A = -sgn*0.5*eye(b.N) + (A1+A2)/2;           % 2-sided avg Nyst + JR Id term
  else, A = A1;        % 1-sided
  end
  % cond(A), figure;plot(eig(A),'+'); axis equal  % slow: 2-sided much better

  tic; iter = true;       % pick
  if ~iter
    dens = A\f;    %r = A*dens-f; fprintf('A resid nrm %g\n',norm(r))
  else, gtol = 1e-9;  maxit = 1e3; tic
    [dens,flag,relres,iter] = gmres(A,f,b.N,gtol,min(maxit,b.N));
    fprintf('\tA GMRES %g sec, %d iters, rel resid nrm %g\n',toc,iter(2),relres)
  end

  co = q.qfsco(dens);   % qfs src co

  if verb && i==numel(Nvs), figure(2); clf; oo.nofig=1;
    showsurffunc(b,dens,oo); hold on; plot3(trg.x(1,:),trg.x(2,:),trg.x(3,:),'k*');
    plot3([0 Edir(1)],[0 Edir(2)],[0 Edir(3)],'m-');
    plot3(q.s.x(1,:),q.s.x(2,:),q.s.x(3,:),'r.'); plot3(q.c.x(1,:),q.c.x(2,:),q.c.x(3,:),'g.');
    if exist('z0','var'), plot3(z0(1),z0(2),z0(3),'r*'); end  % inc src
    title('dens, trgs, src, chk, Edir');
    figure(3); clf; oo.nofig=1; showsurffunc(q.s,co,oo); title('qfs co');
  end
  
  u = qfsker(trg,q.s) * co;   % QFS pot @ trg
  uinctrg = uinc(trg.x)';     % col vec, otherwise matlab auto-bsxfuns! stupid
  utot = uinctrg + u;         % physical soln @ trgs (col vec)
  ufarplain = lpker(trg,b) * dens;    % plain Nystrom eval, far only
  fprintf('BVP [%3d,%3d]: N=%d, Nf=%d\tut farpl,far,nr = %.9g,%.9g,%.9g\n',N(1),N(2),b.N,q.bf.N,ufarplain(1)+uinctrg(1),utot(1),utot(2))
  ub = B*co;
  fprintf('\tsurf u err rms %.3g\t\trms dens %.3g\trms co %.3g\n',sqrt(mean((ub-f).^2)),sqrt(mean(dens.^2)),sqrt(mean(co.^2)))
  
  comfs = B\f; %r = B*co-f; fprintf('mfs resid nrm %g\n',norm(r)) % also check
  umfs = qfsker(trg,q.s) * comfs;    % MFS pot trg
  utotmfs = uinctrg + umfs;
  fprintf('\t\t\t\t\t      MFS ut far,nr = %.9g,%.9g\n',utotmfs(1),utotmfs(2))
end
