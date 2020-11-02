function [u_fi,theta,logl,new_step_mod]=ns_slice(obs,u_fi,logl,invprior,logl_func,logLstar,step_mod,new_step_mod,nsteps,mode,nrepeats)
for m=1:nrepeats
   for n=1:length(u_fi)
     ub={u_fi,u_fi};
     ub{1}(n)=u_fi(n)-step_mod(n)*rand;
     ub{2}(n)=ub{1}(n)+step_mod(n);
     if ~isequal(mode,'ns')
       logLstar=log(rand)+logl;
     end
     for k=1:2
       while logLstar<logl_func(obs,invprior(mod(ub{k},1),obs)) && (ub{2}(n)-ub{1}(n))<1
         ub{k}(n)=ub{k}(n)+(2*k-3)*step_mod(n);
       end
     end
     cond=false;
     un_old=u_fi(n);
     while ~cond
       un_new=ub{1}(n)+(ub{2}(n)-ub{1}(n))*rand;
       u_fi(n)=mod(un_new,1);
       theta=invprior(u_fi,obs);
       logl=logl_func(obs,theta);
       cond=(logl >= logLstar);
       if ~cond
         if un_new>un_old
           ub{2}(n)=un_new;
         else
           ub{1}(n)=un_new;
         end
       end
     end
     new_step_mod(n)=min(new_step_mod(n)*exp(-sign(new_step_mod(n)-(ub{2}(n)-ub{1}(n)))/nsteps),1);
   end
end
end
