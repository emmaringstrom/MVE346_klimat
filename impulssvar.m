% 
% function imp = impulssvar(t,A,tao,n)
%     
%     imp = 0;
%     tao1 = zeros(length(t),n);
%     for j = t
%         tao1(j+1,:) = tao(j+1);
%         for i = 1:n
%             imp = imp + A(i) * exp(-t/tao1(j+1,i));
%         end
%     end
%     
% end

function imp = impulssvar(t,A,tao,n)
    
    imp = sum(A(1:n) * exp(-t./tao(t)));

    
end

