function [alpha] = andi_get_seg_levy_via_bgof_flex(alpha,results,pos,factor)

for k=1:2
  cond=false;
  for i=1:length(alpha.seg(k).bgof.logZmod)
%    len=(3-2*k)*alpha.tsw.median+(k-1)*size(obs,1);
%    if alpha.seg(k).bgof.logZmod(i)/log(10)>20*len/200  
    if alpha.seg(k).bgof.logZmod(i)/log(10)>factor
      cond=true;
    end
  end

  if cond
    alp=struct;
    alp.best_model=3;
    alp.model_probs=[0 0 0 1 0];
    if k==1
        alp.median=tamsd(pos(1:(alpha.tsw.median+1),:));
        alp.mean=alp.median;
        alp.cond_mean=alp.median;
        alp.cond_median=alp.median;
        alp.std_dev=0;
        alpha.seg1=alp;
    else
        alp.median=tamsd(pos((alpha.tsw.median+1):end,:));
        alp.mean=alp.median;
        alp.cond_mean=alp.median;
        alp.cond_median=alp.median;
        alp.std_dev=0;
        alpha.seg1=alp;
        alpha.seg2=alp;
    end
  end
end
out=[alpha.tsw.median,alpha.seg1.best_model,alpha.seg1.median,alpha.seg2.best_model,alpha.seg2.median];
alpha.result_vector=out;
%fprintf('Updated result vector: [%.2f, %i, %.3f, %i, %.3f]\n',out(1),out(2),out(3),out(4),out(5))
end
