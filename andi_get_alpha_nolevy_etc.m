function [alpha] = andi_get_alpha_nolevy_etc(results)

nsamp = length(results.samples);
alphas = NaN(1,nsamp);
model_no = NaN(1,nsamp);
probs = NaN(1,nsamp);
model_probs=zeros(1,5);


for i=1:nsamp
  probs(i) = exp(results.samples(i).logp);
  model_no(i)=results.samples(i).theta{2}{1}-1;
  model_probs(model_no(i)+1)=model_probs(model_no(i)+1)+probs(i);
  if ismember(model_no(i),[0])
    alphas(i) = results.samples(i).theta{2}{2}{1}.meta(1);
  elseif ismember(model_no(i),[1]) %Use if ctrw is "joined"
    alphas(i) = results.samples(i).theta{2}{2}{2}{1}.meta(1);
  elseif ismember(model_no(i),[2 3])
    alphas(i)= results.samples(i).theta{2}{2}(1);
%   elseif ismember(model_no(i),[3]) % Disc
% %    alphas(i)= results.samples(i).theta{2}{2}(1); %Use when FBM included
%     alphas(i)=1.5; %use when filler model included
% %    levy_model=results.samples(i).theta{2}{2}{1};
% %    alphas(i)= results.samples(i).theta{2}{2}{2}{1}.meta(1); %When Levy included
  else
    fprintf('No model assigned for sample %i when running andi_get_alpha\n',i)
  end
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
model_probs(5)=model_probs(4);
model_probs(4)=0;
alpha.model_probs=model_probs;
[~,best]=max(model_probs);
alpha.best_model=best-1;

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

