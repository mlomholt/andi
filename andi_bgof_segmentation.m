function alpha = andi_bgof_segmentation(obs,seg_model,alpha)
if exist('bgof_main')
  if ~isfield(seg_model,'bgof')
    seg_model.bgof=struct;
  end
  if ~isfield(seg_model.bgof,'nsamples')
    seg_model.bgof.nsamples=12;
  end
  for k=1:2
    [few_samples,chosen]=bgof_resample(alpha.seg(k).samples,seg_model.bgof.nsamples);
    tswitchs=alpha.tswitchs(chosen);
    tswitchs_un=unique(tswitchs);
    for i=1:length(tswitchs_un)
      ts=[1 tswitchs_un(i); tswitchs_un(i)+1 Inf];
      cur_samples=few_samples(tswitchs==tswitchs_un(i));
      tprob=sum([cur_samples(:).post]);
      for j=1:length(cur_samples)
        cur_samples(j).post=cur_samples(j).post/tprob;
      end
      bgof=bgof_main(obs(ts(k,1):min(ts(k,2),end),:),cur_samples,seg_model.opt.x2u,seg_model.bgof);
      if i==1
        logZmod=bgof.logZmod+log(tprob);
      else
        logZmod=bmc_logsumexpcol([logZmod; (bgof.logZmod+log(tprob))]);
      end
    end
    alpha.seg(k).bgof.logZmod=logZmod;
  end
end

