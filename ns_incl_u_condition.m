function u = ns_incl_u_condition(u,obs,cond,ufunc,genu)
  u=ufunc(u,obs);
  while ~cond(u,obs)
    u = genu(obs);
  end
end
