function an = ns_analyze(samples,model,misc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates percentiles, maximum likelihood parameters,
% mean and standard deviation for the parameters
% in samples. If model.add contains a function of the parameters, then
% percentiles etc. is also calculated for them.
%
% Contributors to the code in this file: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   n_samp = length(samples);
   left=logical(ones(1,n_samp));
   un_discs=cellfun(model.disc,{samples(:).theta},'UniformOutput',false);
   un_conts=cellfun(model.cont,{samples(:).theta},'UniformOutput',false);
   discs={};
   thetas={};
   logps={};
   logls={};
   while sum(left)>0
     next=find(left,1);
     discs{end+1}=un_discs{next};
     indices=cellfun(@(x)isequal(x,un_discs{next}),un_discs);
     thetas=horzcat(thetas,cell2mat({un_conts{indices}}')');
     left=left-indices;
     logps=horzcat(logps,cell2mat({samples(indices).logp}));
     logls=horzcat(logls,cell2mat({samples(indices).logl}));
   end 

   percentiles={};
   param_mean = {};
   param_stddev = {};
   maxLpar={};
   maxlogls={};
   log_sumps={};

   if isfield(model,'add')
     for j=1:length(model.add)
       add=NaN(1,n_samp);
       for k=1:n_samp
         add(k)=model.add{j}(samples(k).theta);
       end
       thetas=vertcat(thetas,add);
     end
   end

for i=1:length(logps)
   n_theta=length(thetas{i}(:,1));
   p_mean = zeros(1,n_theta);
   p_var = zeros(1,n_theta);
   log_sump = ns_logsumexp(logps{i});
   posterior = exp(logps{i}-log_sump);
   for j=1:n_theta;
     p_mean(j)=sum(thetas{i}(j,:).*posterior);
     p_var(j)=sum((thetas{i}(j,:)-p_mean(j)).^2.*posterior);
   end
   param_mean = horzcat(param_mean,p_mean);
   param_stddev = horzcat(param_stddev,sqrt(p_var(:)));
   percentiles=horzcat(percentiles,ns_percentiles(misc.percentiles_at,thetas{i},posterior));
   [maxlogl,jmax]=max(logls{i});
   maxLpar=horzcat(maxLpar,thetas{i}(:,jmax)');
   maxlogls=horzcat(maxlogls,maxlogl);
   log_sumps=horzcat(log_sumps,log_sump);
end
an.percentiles=percentiles;
an.param_mean=param_mean;
an.param_stddev=param_stddev;
an.maxLpar=maxLpar;
an.maxlogls=maxlogls;
an.discs=discs;
an.log_sumps=log_sumps;
end

%%% Function that calculates percentiles %%%
function percentiles = ns_percentiles(p_list,theta_list,posterior)
  ntheta=length(theta_list(:,1));
  percentiles=NaN(ntheta,length(p_list));
  for j=1:ntheta
    [theta_sort, I] = sort(theta_list(j,:));
    post_sort=posterior(I);
    post_cum=[0 cumsum(post_sort) 1];
    theta_sort=[theta_sort(1) theta_sort theta_sort(end)];
    n=2;
    for m=1:length(p_list)
      while p_list(m)>post_cum(n)
        n=n+1;
      end
      percentiles(j,m)= (theta_sort(n)*(p_list(m)-post_cum(n-1))+theta_sort(n-1)*(post_cum(n)-p_list(m)))/(post_cum(n)-post_cum(n-1));
    end
  end
end

