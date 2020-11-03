function [alpha] = andi_get_levy_via_bgof_flex(alpha,results,obs,factor)

cond=false;
for i=1:length(results.bgof.logZmod)
  if results.bgof.logZmod(i)/log(10)>factor
    cond=true;
  end
end

if cond
  alpha.best_model=3;
  alpha.model_probs=[0 0 0 1 0];
  alpha.median=1.5;
  alpha.mean=alpha.median;
  alpha.cond_mean=alpha.median;
  alpha.cond_median=alpha.median;
end
