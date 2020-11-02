function u=ns_put_u_i(u,ui,ending)
  s=fieldnames(u);
  le=length(ending);
  for j=1:length(s)
    if length(s{j})>=le && isequal(s{j}((end-le+1):end),ending)
      u=rmfield(u,s{j});
    end
  end
  si=fieldnames(ui);
  for j=1:length(si)
    u.([si{j} ending])=ui.(si{j});
  end
end

