function [alpha] = andi_get_alpha_etc(results,get_model,get_alpha)

nsamp = length(results.samples);
alphas = NaN(1,nsamp);
model_no = NaN(1,nsamp);
probs = NaN(1,nsamp);
model_probs=zeros(1,5);

for i=1:nsamp
  probs(i) = exp(results.samples(i).logp);
  model_no(i)=get_model(results.samples(i).theta);
  model_probs(model_no(i)+1)=model_probs(model_no(i)+1)+probs(i);
  alphas(i)= get_alpha(results.samples(i).theta,model_no(i));
end
probs=probs/sum(probs);
model_probs=model_probs/sum(model_probs);
alpha.median=ns_median(alphas,probs);
alpha.mean = sum(alphas.*probs);
alpha_var = sum((alphas-alpha.mean).^2.*probs);
alpha.stddev=sqrt(alpha_var);
alpha.model_probs=model_probs;
[~,best]=max(model_probs);
alpha.best_model=best-1;
i_bests=find(model_no+1==best);
best_probs=probs(i_bests);
best_probs=best_probs/sum(best_probs);
best_alphas=alphas(i_bests);
alpha.cond_mean=sum(best_alphas.*best_probs);
alpha.cond_median=ns_median(best_alphas,best_probs);



end
%--------------
function fractile = ns_median(theta_list,posterior)
    p=0.5;
    [theta_sort, I] = sort(theta_list);
    post_sort=posterior(I);
    post_cum=[0 cumsum(post_sort) 1];
    theta_sort=[theta_sort(1) theta_sort theta_sort(end)];
    n=2;
    while p>post_cum(n)
      n=n+1;
    end
    fractile = (theta_sort(n)*(p-post_cum(n-1))+theta_sort(n-1)*(post_cum(n)-p))/(post_cum(n)-post_cum(n-1));
    end

