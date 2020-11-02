function [alpha] = andi_get_tswitch_etc(results,get_model,get_alpha)

nsamp = length(results.samples);
segmented=logical(zeros(1,nsamp));

for i=1:nsamp
  probs(i) = exp(results.samples(i).logp);
  segmented(i)=(length(results.samples(i).theta{1})>0);
end
prob_seg=sum(probs(segmented));
fprintf('Probability of two segments: %.3f\n',prob_seg);
if prob_seg>=0.5
  results.samples=results.samples(segmented);
  probs=probs(segmented);
  for i=1:length(results.samples)
    for k=1:2
      seg(k).samples(i).theta=results.samples(i).theta{k+1};
      seg(k).samples(i).logp=results.samples(i).logp;
      seg(k).samples(i).post=results.samples(i).post;
    end
    tswitchs(i)=results.samples(i).theta{1};
  end
  alpha.seg1=andi_get_alpha_etc(seg(1),get_model,get_alpha);
  alpha.seg2=andi_get_alpha_etc(seg(2),get_model,get_alpha);
  alpha.tsw.median=ns_median(tswitchs,probs);
  alpha.tsw.mean = sum(tswitchs.*probs);
  alpha_var = sum((tswitchs-alpha.tsw.mean).^2.*probs);
  alpha.tsw.stddev=sqrt(alpha_var);
  %--- Code added for allowing bgof test later ---
  if exist('bgof_main')
    alpha.seg=seg;
    alpha.tswitchs=tswitchs;
  end
  %--- End of bgof preparation code ---
  out=[alpha.tsw.median,alpha.seg1.best_model,alpha.seg1.median,alpha.seg2.best_model,alpha.seg2.median];
  alpha.result_vector=out;
else
  results.samples=results.samples(~segmented);
  for i=1:length(results.samples)
    results.samples(i).theta=results.samples(i).theta{2};
  end
  alpha=andi_get_alpha_etc(results,get_model,get_alpha);
  out=[0 alpha.best_model,alpha.median,alpha.best_model,alpha.median];
  alpha.result_vector=out;
end
fprintf('Result vector: [%.2f, %i, %.3f, %i, %.3f]\n',out(1),out(2),out(3),out(4),out(5))
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

