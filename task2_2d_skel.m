%inference for task 2 in 2 dimension

tstart=tic;
clear models obs alpha results
midpath='task2_2d_';%for task2 2dimension
datapath='/Datafiles/Task2/';
destpath='/Submit/Task2_2/';
dim=2;
part_total=13309;
maxpoints=200; %cuts longer trajectories at 200 points
inv_cauchy=@(u) tan(pi*(u-1/2));
ctrw_joined=ns_join_models([andi_ctrw_model() andi_ctrw_model(0)]);
models{1}=ns_join_models([andi_attm_model() ctrw_joined andi_fbm_model() andi_sbm_model()]);
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
    alpha=andi_get_alpha_nolevy_etc(results);
    alpha=andi_get_levy_via_bgof_flex(alpha,results,obs,3.5);
%     if(alpha.model_probs(4)==1)
%         alpha.median=tamsd(pos);
%     end 
    outmat=[dim alpha.model_probs];
    filename1=sprintf('%d%s',ind,'.txt');
    dlmwrite(fullfile(destpath,filename1),outmat,'delimiter','\t');
            
    
end

toc(tstart)