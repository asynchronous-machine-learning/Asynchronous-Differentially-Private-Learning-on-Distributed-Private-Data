function loss = lossf(D,Y,theta,sigma)

loss = sum((Y'-theta*D').^2)/length(Y)+(sigma/2)*theta*theta';


