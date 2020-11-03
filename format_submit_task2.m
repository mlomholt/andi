function []=format_submit_task2()
% creates result file task2.txt
datapath1=strcat(pwd,'/Submit/Task2_1/');
datapath2=strcat(pwd,'/Submit/Task2_2/');
datapath3=strcat(pwd,'/Submit/Task2_3/');
destpath=strcat(pwd,'/Submit/');
outmat=[];
for jj=1:16618
    sourcename=sprintf('%d%s',jj,'.txt');
    A=load(fullfile(datapath1,sourcename));
    outmat=[outmat;A];
end    
for jj=1:13309
    sourcename=sprintf('%d%s',jj,'.txt');
    A=load(fullfile(datapath2,sourcename));
    outmat=[outmat;A];
end    
for jj=1:10000
    sourcename=sprintf('%d%s',jj,'.txt');
    A=load(fullfile(datapath3,sourcename));
    outmat=[outmat;A];
end    

fname='task2.txt'; 
dlmwrite(fullfile(destpath,fname),outmat,'delimiter',';');

end