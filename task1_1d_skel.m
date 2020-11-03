%inference for task 1 in 1 dimension
tstart=tic;
clear models obs alpha results
midpath='task1_1d_';
datapath='/Datafiles/Task1/'; %source folder
destpath='/Submit/Task1_1/'; %destination folder
dim=1;
part_total=16618;
maxpoints=200; %cuts longer trajectories at 200 points
inv_cauchy=@(u) tan(pi*(u-1/2));
ctrw_joined=ns_join_models([andi_ctrw_model() andi_ctrw_model(0)]);
levy_joined=andi_hmm_levy_model();
models{1}=ns_join_models([andi_attm_model() ctrw_joined andi_fbm_model() levy_joined andi_sbm_model()]);
models{1}=ns_incl_scaling(models{1},@(u) 1./abs(inv_cauchy(u)));
models{1}.opt.hmc_get_u=@(u) [];
models{1}.opt.slice_head=[models{1}.opt.slice_head {'crg_'}];
options=struct;
options.nsteps=8;
options.nsamples=180;
options.ais_nwalkers=30;
options.nwalkers=24;
options.ais_betas=0.025:0.025:1;
normal_logl=@(x) -log(2*pi)*numel(x)/2-sum(sum(x.^2))/2;
options.beta_logl='sequential';
models{1}.options=options;

parfor ind=1:part_total    
    sourcename=sprintf('%s%d%s',midpath,ind,'.txt');
    pos=load(fullfile(datapath,sourcename));
    
    T_total=size(pos,1);
    if T_total>maxpoints
       obs=diff(pos(1:maxpoints,:));
    else
       obs=diff(pos);
    end

    misc=struct('mode','mc');
    fprintf('Starting inference for trajectory %i\n',ind);
    results=ns_main(obs,models,misc);
    fprintf('Finished trajectory %i\n',ind);
    alpha=andi_get_alpha_etc(results,@get_model,@get_alpha);
    outmat=[dim alpha.median];
%    outmat2=[dim alpha.model_probs];
    filename1=sprintf('%d%s',ind,'.txt');
    dlmwrite(fullfile(destpath1,filename1),outmat,'delimiter','\t');
     
end

toc(tstart)
%----------
function model_no = get_model(theta)
  model_no=theta{2}{1}-1;
end

function alpha = get_alpha(theta,model_no)
  if ismember(model_no,[0])
    alpha = theta{2}{2}{1}.meta(1);
  elseif ismember(model_no,[1]) 
    alpha = theta{2}{2}{2}{1}.meta(1);
  elseif ismember(model_no,[2 4])
    alpha= theta{2}{2}(1);
  elseif ismember(model_no,[3]) 
    alpha= theta{2}{2}(1); 
  end
end

