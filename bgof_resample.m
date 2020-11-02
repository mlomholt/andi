function [new_samples,varargout] = bgof_resample(samples,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  post_cum=[0 cumsum([samples(:).post])];
  post_cum=post_cum/post_cum(end);
  us=((0:(N-1))+rand)/N;
  counts = histcounts(us,post_cum);
  chosen=(counts>0);
  new_samples = samples(chosen);
  posts=num2cell(counts(chosen)/N);
  [new_samples.post]=posts{:};
  if nargout==2
    varargout{1}=chosen;
  end
end

