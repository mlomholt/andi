function u = wt_genu(obs,dt_time_invprior,genu_meta,invprior_meta)
  T_total=size(obs,1);
  u=genu_meta(obs);
  theta_meta=invprior_meta(u);
  u.ds_t=rand;
  t=dt_time_invprior(u.ds_t,theta_meta);
  while t<T_total
    u.ds_t=[u.ds_t rand];
    t=t+dt_time_invprior(u.ds_t(end),theta_meta);
  end
end

