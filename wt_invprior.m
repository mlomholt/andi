function theta = wt_invprior(u,obs,dt_dec_invprior,invprior_meta)
  theta.meta = invprior_meta(u,obs);
  theta.dts = arrayfun(@(u) dt_dec_invprior(u,theta.meta),u.ds_t); % theta.dts(i) is the length of the i'th time interval, or maybe after a floor command.
end

