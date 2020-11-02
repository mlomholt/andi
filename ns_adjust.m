function uout = ns_adjust(u,uout,fields,expected)
%%%%%%%%%%%%%%%%%
% Function that adjusts the fields in 'fields' of u to have
% the lengths in 'expected' and attaches them to uout.
%
% Contributors to the code in this file: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%

  for i=1:length(fields)
    if isfield(u,fields{i}) && length(u.(fields{i}))>0
      missing=expected(i)-length(u.(fields{i}));
      if missing>0
        uout.(fields{i})=[u.(fields{i}) rand(1,missing)];
      else
        uout.(fields{i})=u.(fields{i})(1,1:expected(i));
      end
    else
      uout.(fields{i})=rand(1,expected(i));
    end
  end
end
