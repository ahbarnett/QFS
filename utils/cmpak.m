function z = cmpak(v,ncomp)          % helper:  if ncomp=2 pack pairs as C-#s
if ncomp==1, z=v;
elseif ncomp==2, z=v(1:end/2)+1i*v(end/2+1:end);
else, error('cmpak only for ncomp=1 or 2!');
end
