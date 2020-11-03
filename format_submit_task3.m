function []=format_submit_task3()
% creates result file task3.txt
datapath1=strcat(pwd,'/Submit/Task3_1/');
datapath2=strcat(pwd,'/Submit/Task3_2/');
datapath3=strcat(pwd,'/Submit/Task3_3/');
destpath=strcat(pwd,'/Submit/');
outmat=[];
for jj=1:10000
    sourcename=sprintf('%d%s',jj,'.txt');
    A=load(fullfile(datapath1,sourcename));
    outmat=[outmat;A];
end    
for jj=1:10000
    sourcename=sprintf('%d%s',jj,'.txt');
    A=load(fullfile(datapath2,sourcename));
    outmat=[outmat;A];
end    
for jj=1:10000
    sourcename=sprintf('%d%s',jj,'.txt');
    A=load(fullfile(datapath3,sourcename));
    outmat=[outmat;A];
end    

fname='task3.txt'; 
dlmwrite(fullfile(destpath,fname),outmat,'delimiter',';');

end