function[alph]=tamsd(A)
%calculates TAMSD
del_time=1; %time difference between frames in seconds.
fit_points=9;
L1 = length(A(:,1)); 
total_lags=9;
%dim=length(A(1,:));

lagtime = zeros(total_lags,1);

Tmsd =zeros(total_lags,1);

if(L1>=10)
    for i = 1:total_lags
        lagtime(i)= del_time * i ;
    end

    for jj=1:total_lags

        T_min_del_inv = 1/(L1-jj);

        for kk = 1:(L1-jj)
           Tmsd(jj)= Tmsd(jj) + T_min_del_inv *sum((A(kk+jj,:)-A(kk,:)).^2);
        end

    end 
    B=log([lagtime Tmsd]);
    p=polyfit(B(1:fit_points,1),B(1:fit_points,2),1);
    alph=p(1);
else
    alph=1.5;
end    
end
