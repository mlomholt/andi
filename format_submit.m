function [outmat]=format_submit()

datapath=strcat(pwd,'/Submit/Task2_1/');
destpath=strcat(pwd,'/Submit/');
parttotal=10000;
dim=1;
%outmat=[dim*ones(parttotal,1) zeros(parttotal,1)]; % for task1
outmat=[dim*ones(parttotal,1) zeros(parttotal,5)]; %for task2
for jj=1:parttotal
%    jj
    sourcename=sprintf('%d%s',jj,'.txt');
    A=load(fullfile(datapath,sourcename));
    outmat(jj,:)=A;
end    

fname='submit_task2_1.txt'; %change according to task
dlmwrite(fullfile(destpath,fname),outmat,'delimiter',';');

end