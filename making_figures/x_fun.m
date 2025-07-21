function [x,Dx]=x_fun(t)
  % x=t; Dx=1;
 % x=sin(2*pi*t)/(2*pi)+t; Dx=cos(2*pi*t)+1;
    x=-1+(sin(2*pi*t)/(2*pi)+sin(3.1*pi*t)/(3.1*pi))/2+t; Dx=(cos(2*pi*t)+cos(3.1*pi*t))/2+1;
end