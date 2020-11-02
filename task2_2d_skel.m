tstart=tic;
clear models obs alpha results
done=[];
%done=load('undone1_task2_2d.txt');
midpath='task2_2d_';%for task2 2dimension
datapath='Datafiles/Task2/';
destpath='Submit/Task2_2/';%for task2 2 dimension
dim=2;
partint=1;
part_total=2;
maxpoints=200; %cuts longer trajectories at 200 points
%outmat1=[dim*ones(part_total,1) zeros(part_total,1)];
%outmat2=[dim*ones(part_total,1) zeros(part_total,5)];
inv_cauchy=@(u) tan(pi*(u-1/2));
ctrw_joined=ns_join_models([andi_ctrw_model() andi_ctrw_model(0)]);
models{1}=ns_join_models([andi_attm_model() ctrw_joined andi_fbm_model() andi_sbm_model()]);
models{1}=ns_incl_scaling(models{1},@(u) 1./abs(inv_cauchy(u)));
models{1}.opt.hmc_get_u=@(u) [];
models{1}.opt.slice_head=[models{1}.opt.slice_head {'crg_'}];

options=struct;
options.nsteps=10;
options.ais_nwalkers=30;
options.nwalkers=24;
options.nsamples=8*options.nwalkers;
options.ais_betas=0.025:0.025:1;
normal_logl=@(x) -log(2*pi)*numel(x)/2-sum(sum(x.^2))/2;
options.beta_logl='sequential';
models{1}.options=options;

% presumm=['../../Resultfiles/summary_' midpath 'results_'];
% premat=['../../Resultfiles/saved_' midpath 'results_'];

%if (isempty(gcp('nocreate')))
%     parpool(12)
%end 
Remmat=[];
for jj=partint:part_total
     if(ismember(jj,done))

     else
        Remmat=[Remmat;jj];
     end
end
L=length(Remmat);
parfor(kk=1:L)
    ind=Remmat(kk);
    sourcename=sprintf('%s%d%s',midpath,ind,'.txt');
    pos=load(fullfile(datapath,sourcename));
    T_total=size(pos,1);
    if T_total>maxpoints
        obs=diff(pos(1:maxpoints,:));
    else
        obs=diff(pos);
    end
    if(abs(obs(1,1))>10^9||abs(obs(1,2))>10^9)
        outmat2=[dim 0 0 0 1 0];
        filename2=sprintf('%d%s',ind,'.txt');
        dlmwrite(fullfile(destpath,filename2),outmat2,'delimiter','\t');
    else
        scale_factor=max(abs(obs));
        scale_factor(scale_factor==0)=1;
        obs=obs./scale_factor;
        %misc=struct('nssummary',[presumm sprintf('%d',i) '.txt'],'save_results',[premat sprintf('%d',i) '.mat']);
        %misc.mode='mc';
        misc=struct('mode','mc');
        fprintf('Starting inference for trajectory %i\n',ind);
        results=ns_main(obs,models,misc);
        fprintf('Finished trajectory %i\n',ind);
        alpha=andi_get_alpha_nolevy_etc(results);
        alpha=andi_get_levy_via_bgof3(alpha,results,obs,3.5);
        outmat2=[dim alpha.model_probs];
        filename2=sprintf('%d%s',ind,'.txt');
        dlmwrite(fullfile(destpath,filename2),outmat2,'delimiter','\t');
    end
end

toc(tstart)

