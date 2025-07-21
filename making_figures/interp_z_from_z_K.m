function z=interp_z_from_z_K(t,t_K,h_K,z_K,p_K)
%
 z=zeros(size(t));
 K=length(z_K);
%
 i=1; % start from interval 1: (t_1,t_2)
%
 for j=1:length(t)
    if i<K-1
        if t(j)>t_K(i+1)
            i=i+1; % interval "i": (t_i,t_{i+1})
        end
    end
  %
     z(j)=z_K(i+1)*(t(j)-t_K(i))^2/h_K(i)^2+...
         +p_K(i)*(t(j)-t_K(i))*(t(j)-t_K(i+1))/h_K(i)^2+...
         z_K(i)*(t(j)-t_K(i+1))^2/h_K(i)^2;
 end
end