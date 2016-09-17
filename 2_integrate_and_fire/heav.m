function y = heav(x)
   temp = zeros(1,length(x));
   temp(x>0) = 1;
   temp(x==0) = 0.5;
   y = temp;
end