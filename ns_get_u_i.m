function ui=ns_get_u_i(u,ending)
  s=fieldnames(u);
  ui=struct;
  le=length(ending);
  for j=1:length(s)
    if length(s{j})>=le && isequal(s{j}((end-le+1):end),ending)
      ui.(s{j}(1:end-le))=u.(s{j});
    end
  end
end

