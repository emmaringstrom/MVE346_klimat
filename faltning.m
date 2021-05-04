
function M = faltning(t,M0,U,T,A,tao)

    s = 0;
    
    %for j = t
    s = s + impulssvar(T-t+1,A,tao,5).*U(t+1);
    %end
    
    M = M0 + s;
end