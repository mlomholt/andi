function step_mod = ns_step_mod_mean(new_step_mod)
    if all(cellfun(@(x)isequal(x,0),new_step_mod))
      step_mod=0;
    elseif any(cellfun(@(x)isstruct(x),new_step_mod))
      step_mod=struct;
      store_evolve_extra_step_mod={};
      for m=1:length(new_step_mod)
        if isfield(new_step_mod{m},'evolve_extra_step_mod')
          store_evolve_extra_step_mod=[store_evolve_extra_step_mod {new_step_mod{m}.evolve_extra_step_mod}];
          new_step_mod{m}=rmfield(new_step_mod{m},'evolve_extra_step_mod');
        end
        if ~isequal(new_step_mod{m},0)
          s=fieldnames(new_step_mod{m});
          for j=1:length(s)
            step_mod=setfield(step_mod,s{j},[]);
          end
        end
      end
      new_step_mod_fi={};
      s=fieldnames(step_mod);

      if length(store_evolve_extra_step_mod)>0
        step_mod.evolve_extra_step_mod=ns_step_mod_mean(store_evolve_extra_step_mod);
      end

      for j=1:length(s)
        for m=1:length(new_step_mod)
          if ~isequal(new_step_mod{m},0) && isfield(new_step_mod{m},s{j})
            new_step_mod_fi{m}=getfield(new_step_mod{m},s{j});
          else
            new_step_mod_fi{m}=[];
          end
        end
        ls=cellfun(@(x) length(x),new_step_mod_fi);
        step_mod_fi=zeros(1,max(ls));
        counts=zeros(1,max(ls));
        for m=1:length(new_step_mod)
          step_mod_fi(1:ls(m))=step_mod_fi(1:ls(m))+log(new_step_mod_fi{m});
          counts(1:ls(m))=counts(1:ls(m))+1;
        end
        for m=1:length(counts)
          step_mod_fi(m)=exp(step_mod_fi(m)/counts(m));
        end
        step_mod=setfield(step_mod,s{j},step_mod_fi);
      end
    elseif any(cellfun(@(x)iscell(x),new_step_mod))
      step_mod={};
      m=0;
      while length(step_mod)==m
        m=m+1;
        store={};
        for ipar=1:length(new_step_mod)
          if length(new_step_mod{ipar})>=m && ~isequal(new_step_mod{ipar},0)
            store=[store {new_step_mod{ipar}{m}}];
          end
        end
        if length(store)>0
          step_mod{m}=ns_step_mod_mean(store);
        end
      end
    else
      step_mod=new_step_mod{1};
      for ipar=2:length(new_step_mod)
        step_mod=step_mod.*new_step_mod{ipar};
      end
      step_mod=nthroot(step_mod,length(new_step_mod));
    end
end
