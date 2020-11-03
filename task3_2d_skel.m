%inference for task 3 in 2 dimension
tstart=tic;
clear models obs two_seg_model alpha results
midpath='task3_2d_';
datapath='/Datafiles/Task3/';
destpath='/Submit/Task3_2/';
dim=2; %dimension
part_total=10000;
ctrw_joined=ns_join_models([andi_ctrw_model() andi_ctrw_model(0)]);
models{1}=ns_join_models([andi_attm_model() ctrw_joined andi_fbm_model() andi_sbm_model()]);
inv_cauchy=@(u) tan(pi*(u-1/2));
seg_model=ns_incl_scaling(models{1},@(u) 1./abs(inv_cauchy(u)));
N_pos=200;
models{1}=ns_2_segments_model(seg_model,@(u,T_total) min(ceil((N_pos-2)*u),T_total-1));

models{1}.opt.hmc_get_u=@(u) [];
models{1}.opt.slice_head=[models{1}.opt.slice_head {'crg_'}];

options=struct;
options.nsteps=8;
options.nsamples=180;
options.ais_nwalkers=30;
options.nwalkers=24;
options.ais_betas=0.025:0.025:1;
normal_logl=@(x) -log(2*pi)*numel(x)/2-sum(sum(x.^2))/2;
models{1}.options=options;

parfor ind=partint:part_total
    sourcename=sprintf('%s%d%s',midpath,ind,'.txt');
    pos=load(fullfile(datapath,sourcename));
   
    T_total=size(pos,1);
    obs=diff(pos);
    misc=struct('mode','mc');
    fprintf('Starting inference for trajectory %i\n',ind);
    results=ns_main(obs,models,misc);
    fprintf('Finished trajectory %i\n',ind);    
    alpha=andi_get_tswitch_nolevy_etc(results);
      if isfield(alpha,'seg')
        alpha=andi_bgof_segmentation(obs,seg_model,alpha);
        alpha=andi_get_seg_levy_via_bgof_flex(alpha,results,pos,3.5);
        outmat=[dim alpha.result_vector];
      else
        alpha=andi_get_levy_via_bgof_flex(alpha,results,pos,3.5);
        if(alpha.model_probs(4)==1)
            alpha.median=tamsd(pos);
        end 
        outmat=[dim 0 alpha.best_model alpha.median alpha.best_model alpha.median];
      end
      filename2=sprintf('%d%s',ind,'.txt');
      dlmwrite(fullfile(destpath,filename2),outmat,'delimiter','\t');
       
    
end

toc(tstart)
