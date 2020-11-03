%inference for task 3 in 1 dimension

tstart=tic;
clear models obs two_seg_model alpha results
midpath='task3_1d_';
datapath='/Datafiles/Task3/';
destpath='/Submit/Task3_1/';
dim=1; %dimension
partint=1;
part_total=10000;
ctrw_joined=ns_join_models([andi_ctrw_model() andi_ctrw_model(0)]);
levy_joined=andi_hmm_levy_model();
models{1}=ns_join_models([andi_attm_model() ctrw_joined andi_fbm_model() levy_joined andi_sbm_model()]);
inv_cauchy=@(u) tan(pi*(u-1/2));
seg_model=ns_incl_scaling(models{1},@(u) 1./abs(inv_cauchy(u)));
N_pos=200;
models{1}=ns_2_segments_model(seg_model,@(u,T_total) min(ceil((N_pos-2)*u),T_total-1));

models{1}.opt.hmc_get_u=@(u) [];
models{1}.opt.slice_head=[models{1}.opt.slice_head {'crg_'}];

options=struct;
options.nsteps=10;
options.ais_nwalkers=30;
options.nwalkers=24;
options.nsamples=8*options.nwalkers;
options.ais_betas=0.025:0.025:1;
normal_logl=@(x) -log(2*pi)*numel(x)/2-sum(sum(x.^2))/2;
models{1}.options=options;




parfor ind=1:part_total
    sourcename=sprintf('%s%d%s',midpath,ind,'.txt');
    pos=load(fullfile(datapath,sourcename)); 
    T_total=size(pos,1);
    obs=diff(pos);
    misc=struct('mode','mc');
    fprintf('Starting inference for trajectory %i\n',ind);
    results=ns_main(obs,models,misc);
    fprintf('Finished trajectory %i\n',ind);    
    alpha=andi_get_tswitch_etc(results,@get_model,@get_alpha);
    outmat=[dim alpha.result_vector];    
    filename2=sprintf('%d%s',ind,'.txt');
    dlmwrite(fullfile(destpath,filename2),outmat,'delimiter','\t');
    
end
%..................................
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



