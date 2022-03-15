
function correlation = Biweight_midcorrelation(M)
    med = median(M,2);
     ma = mad(M,1,2);
     ma = 9*ma;
     fun1 = @(a,b) (a - b) ;
     v = bsxfun(fun1,M,med);
     u = bsxfun(@rdivide,v,ma);

     w = zeros(size(M));
     index = find(abs(u) <1);

     tmp = u(index);
     tmp2 = tmp.^2;
     tmp = 1 - tmp2;
     tmp2 = tmp .^2;
     w(index)= tmp2;

     fun3 = @(a,b) a - b;
     xtilda = bsxfun(fun3, M,med);
     xtil = xtilda .* w;
     xtilda = normr(xtil); 
     correlation = xtilda * xtilda' ; 
end
