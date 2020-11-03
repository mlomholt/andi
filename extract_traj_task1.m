function []=extract_traj_task1()
%extracts individual trajectories for each dimension from AnDi datasets for
%task1
source_id=pwd; %source file location
dest_id = strcat(pwd,'/Datafiles/Task1/'); %destination folder
fname='task1'; 
fid=fopen(fullfile(source_id,strcat(fname,'.txt')));


jj=1;

while jj<10001
    clear A;
    st=fgetl(fid);
    A=str2num(st);
    A=A(2:end);
    filename1=sprintf('%s%s%d%s',fname,'_1d_',jj,'.txt');
    dlmwrite(fullfile(dest_id,filename1),A,'delimiter','\t','precision',12);
    jj=jj+1;
end  
kk=1;
while jj<20001
    clear A;
    st=fgetl(fid);
    A=str2num(st);
    A=A(2:end);
    rownum=length(A)/2;
    A=reshape(A,rownum,2);
    filename1=sprintf('%s%s%d%s',fname,'_2d_',kk,'.txt');
    dlmwrite(fullfile(dest_id,filename1),A,'delimiter','\t','precision',12);
    jj=jj+1;
    kk=kk+1;
end    
 kk=1;
while jj<30001
    clear A;
    st=fgetl(fid);
    A=str2num(st);
    A=A(2:end);
    rownum=length(A)/3;
    A=reshape(A,rownum,3);
    filename1=sprintf('%s%s%d%s',fname,'_3d_',kk,'.txt');
    dlmwrite(fullfile(dest_id,filename1),A,'delimiter','\t','precision',12);
    jj=jj+1;
    kk=kk+1;
end

kk=10001;
while jj<36619 
    clear A;
    st=fgetl(fid);
    A=str2num(st);
    A=A(2:end);
    filename1=sprintf('%s%s%d%s',fname,'_1d_',kk,'.txt');
    dlmwrite(fullfile(dest_id,filename1),A,'delimiter','\t','precision',12);
    jj=jj+1;
    kk=kk+1;
end    
kk=10001;
while jj<39928
    clear A;
    st=fgetl(fid);
    A=str2num(st);
    A=A(2:end);
    rownum=length(A)/2;
    A=reshape(A,rownum,2);
    filename1=sprintf('%s%s%d%s',fname,'_2d_',kk,'.txt');
    %writematrix(A,fullfile(dest_id,filename1),'Delimiter','tab');
    dlmwrite(fullfile(dest_id,filename1),A,'delimiter','\t','precision',12);
    jj=jj+1;
    kk=kk+1;
end    

fclose(fid);


end