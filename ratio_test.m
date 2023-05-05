function [mr,ir] = ratio_test(a,b)
mr = Inf;
ir = 0;
for i = 1:1: length(b)
    if a(i) > 0 
        m = b(i)/a(i);
        if m < mr
            mr = m;
            ir = i;
        end
    end
    
end
end