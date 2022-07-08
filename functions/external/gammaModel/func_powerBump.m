function F = func_powerBump(x,P,f,p_exp)

% F= P - (x(1)-p_exp*f + x(2)*normpdf(f,x(3),.1));

F= P - (x(1)-p_exp*f + x(2)*sqrt(2*pi)*normpdf(f,x(3),x(4)));

% The amplitude is gain/sigma = x(2)/x(4)

