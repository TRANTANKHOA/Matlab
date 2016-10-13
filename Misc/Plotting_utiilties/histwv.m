function [p,h] = histwv(v, w, bins) 
%Inputs: 
%vv - values 
%ww - weights 
%bins - number of bins (inclusive) 

%Outputs: 
%histw - weighted histogram 
%histv (optional) - histogram of values 

delta = (max(v)-min(v))/(bins-1); 
subs = round((v-min(v))/delta)+1; 
p = min(v)+(0:bins-1)*delta;
h = accumarray(subs,w,[bins,1]); 
h = h/sum(h.*(p(2)-p(1)));
end