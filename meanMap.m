function meanmap = meanMap(im,window)

[L1,L2]=size(im);
if isscalar(window), w1=window; w2=window;
else w1=window(1); w2=window(2); end

N1=floor(L1/w1); N2=floor(L2/w2);
rest1=L1-N1*w1; rest2=L2-N2*w2;
rest1l=floor(rest1/2); rest2l=floor(rest2/2);

meanmap=nan([N1 N2]);
for cnt1=1:N1
    for cnt2=1:N2
        sector=im(1+rest1l+(cnt1-1)*w1:rest1l+cnt1*w1,1+rest2l+(cnt2-1)*w2:rest2l+cnt2*w2);
        meanmap(cnt1,cnt2)=mean(sector(:));
    end
end

end