
%Funktion för Euler Framåt

function B = eulerforward(B0,dB)
    % Euler's Method
    % Initial conditions and setup
    h = 1;  % step size
    years = 1:h:736;  % the range of x
    B = zeros(3, length(years));  % allocate the result y
    B(:,1) = B0(:);  % the initial y value
    %n = numel(B);  % the number of y values
    % The loop to solve the DE
    for i=1:735 %n-1
        f = dB;
        B(:,i+1) = B(:,i) + h * f(B(:,i),i);
    end
end