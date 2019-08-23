function sg2 = subgradient(D,Y,theta)

sg2 = 2*(Y'-theta*D')*(-D)/length(Y);



