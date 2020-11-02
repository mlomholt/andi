function uout = wt_adjust_u(u,obs,dt_time_invprior,adjust_meta,invprior_meta)
  T_total=size(obs,1);
  uout=adjust_meta(u,obs);
  theta_meta=invprior_meta(uout,obs);
  if ~isfield(u,'ds_t') || length(u.ds_t)<1
    u.ds_t=rand;
  end
  uout.ds_t=u.ds_t(1);
  t=dt_time_invprior(uout.ds_t,theta_meta);
  while t<T_total
    if length(u.ds_t)>length(uout.ds_t)
      uout.ds_t=[uout.ds_t u.ds_t(length(uout.ds_t)+1)];
    else
      uout.ds_t=[uout.ds_t rand];
    end
    t=t+dt_time_invprior(uout.ds_t(end),theta_meta);
  end
end

