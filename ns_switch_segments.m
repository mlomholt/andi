function [walker_new1,varargout] = ns_switch_segments(obs,model,logLstar,walker1,walker2,mode)

walker_new1 = walker1;
if ~isequal(mode,'ns')
  varargout={walker2};
else
  varargout={};
end

if length(walker1.theta)>2 && length(walker2.theta)>2
  if ~isequal(mode,'ns')
    walker_new2 = walker2;
    u1=ns_get_u_i(walker1.u,'_1');
    u2=ns_get_u_i(walker2.u,'_1');
    walker1.u=ns_put_u_i(walker1.u,u2,'_1');
    walker2.u=ns_put_u_i(walker2.u,u1,'_1');
    walker1.u=model.adjust_u(walker1.u,obs);
    walker2.u=model.adjust_u(walker2.u,obs);
    walker1.theta=model.invprior(walker1.u,obs);
    walker2.theta=model.invprior(walker2.u,obs);
    walker1.logl=model.logl(obs,walker1.theta);
    walker2.logl=model.logl(obs,walker2.theta);
    logLstar=log(rand)+walker_new1.logl+walker_new2.logl;
    if walker1.logl+walker2.logl >= logLstar
      walker_new1=walker1;     % Updates walker_new1
      varargout={walker2};     % Updates walker_new2
    end
  else
    u2=ns_get_u_i(walker2.u,'_1');
    walker1.u=ns_put_u_i(walker1.u,u2,'_1');
    walker1.u=model.adjust_u(walker1.u,obs);
    walker1.theta=model.invprior(walker1.u,obs);
    walker1.logl=model.logl(obs,walker1.theta);
    if walker1.logl>=logLstar
      walker_new1=walker1;
    else
      walker1=walker_new1;
      u2=ns_get_u_i(walker2.u,'_2');
      walker1.u=ns_put_u_i(walker1.u,u2,'_2');
      walker1.u=model.adjust_u(walker1.u,obs);
      walker1.theta=model.invprior(walker1.u,obs);
      walker1.logl=model.logl(obs,walker1.theta);
      if walker1.logl>=logLstar
        walker_new1=walker1;
      end
    end
  end
end
end

