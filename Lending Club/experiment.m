clc
clear

%#ok<*MINV>
%#ok<*NBRAK>
%#ok<*SAGROW>

%% Reading converted data

XX=csvread('X.csv');
YY=csvread('Y.csv');

%% Sanitizing the data (removing Nan values, etc)

XX(21391,:)=[];
YY(21391,:)=[];

[I,J]=find(isnan(XX));

XX(I,J)=zeros(length(I),length(J));

%% PCA transformation

m=10;
[V,D]=eigs(XX(end-10000:end,:)'*XX(end-10000:end,:),m);                    % as described in the paper the PCA is only 
                                                                           % based on the last 10,000 entries

X=XX*V*inv(sqrt(D));
Y=YY;

%% Generating Figure 3

T=1e3;

n_attributes = size(X,2);

sigma=1e-5;

N=3;                                                                       % 3 data owners

n1=250000;
n2=n1;
n3=n1;

eps1=1;                                                                    % can be replaced with any other 
                                                                           % epsilon, e.g., 10 or 0.1
eps2=eps1;
eps3=eps1;

n_bar = n1+n2+n3;

D1 = X(1:n1,:);                                                            % The first data owner
Y1 = Y(1:n1,:);                                                            % The first data owner

D2 = X(n1+1:n1+n2,:);                                                      % The second data owner
Y2 = Y(n1+1:n1+n2,:);                                                      % The second data owner

D3 = X(n1+n2+1:n1+n2+n3,:);                                                % The third data owner
Y3 = Y(n1+n2+1:n1+n2+n3,:);                                                % The third data owner

Dtrain = X(1:n1+n2+n3,:);                                                  % All training data, i.e., D1, D2, D3
Ytrain = Y(1:n1+n2+n3,:);                                                  % All training data, i.e., Y1, Y2, Y3

no_privacy_theta=inv((2/n_bar)*(Dtrain')*Dtrain ...
                 +sigma*eye(n_attributes))*(2/n_bar)*(Dtrain'*Ytrain);     % ML model without privacy

itmax=100;                                                                 % Number of runs for generating statistics

rtrain_ave=zeros(itmax,T);
for it=1:itmax
    % initialize theta
    theta    = rand(1,n_attributes);
    theta_1  = rand(1,n_attributes);
    theta_2  = rand(1,n_attributes);
    theta_3  = rand(1,n_attributes);
    
    for k = 1:T
        rtrain_ave(it,k)=(lossf(Dtrain,Ytrain,theta,sigma) ...
                          -lossf(Dtrain,Ytrain,no_privacy_theta',sigma)) ...
                          /lossf(Dtrain,Ytrain,no_privacy_theta',sigma);
        pk = 3e3*N/T^2/sigma;                                              % step size (or training rate)
        
        switch ceil(3*rand(1,1))
            case 1
                bar_theta=(theta_1+theta)/2;
                D1_sg = subgradient(D1,Y1,bar_theta); 
                b1 = T*n_attributes/(n1*eps1);

                w1 = b1*log(rand(1,1)/rand(1,1)); 

                theta_1    = projectoperator(bar_theta ...
                                -(sigma*bar_theta)*pk/(2*N) ...
                                -n1*(D1_sg+w1)*pk/n_bar    ,1e5);

                theta      = projectoperator(bar_theta ...
                                -(sigma*bar_theta)*pk/(N*2),1e5);
            case 2
                bar_theta=(theta_2+theta)/2;

                D2_sg = subgradient(D2,Y2,bar_theta);
                b2 = T*n_attributes/(n2*eps2);
                w2 = b2*log(rand(1,1)/rand(1,1)); 

                theta_2    = projectoperator(bar_theta ...
                                -(sigma*bar_theta)*pk/(2*N) ...
                                -n1*(D2_sg+w2)*pk/n_bar    ,1e5);

                theta      = projectoperator(bar_theta ...
                                -(sigma*bar_theta)*pk/(N*2),1e5);               
            case 3
                bar_theta=(theta_3+theta)/2;

                D3_sg = subgradient(D3,Y3,bar_theta);
                b3 = T*n_attributes/(n3*eps3);
                w3 = b3*log(rand(1,1)/rand(1,1));

                theta_3    = projectoperator(bar_theta ...
                                -(sigma*bar_theta)*pk/(2*N) ...
                                -n1*(D3_sg+w3)*pk/n_bar    ,1e5);                    

                theta      = projectoperator(bar_theta ...
                                -(sigma*bar_theta)*pk/(N*2),1e5);
        end
    end
end


figure('Renderer', 'painters', 'Position', [600 300 500 200])
boxplot(rtrain_ave,'symbol','','Widths',0.3,'BoxStyle','filled','Whisker',0)
set(gca,'YScale','log')
set(gca,'XScale','log')
set(gca,'xtick',[1 10 100 1000]);
grid on
hold on
ylim([1e-3 1e1])
yticks(10.^[-3:1:1]) 
xlim([1 T])
semilogy(median(rtrain_ave),'-k','linewidth',2)
set(gca, 'XTickLabels', {'$10^0$','$10^1$','$10^2$','$10^3$'},'TickLabelInterpreter','latex')
a = get(get(gca,'children'),'children');
set(a{2}, 'Color', [0.6627    0.6627    0.6627]);
set(gca,'FontSize', 17)

%% Generating Figures 5 and 6

N=ceil(logspace(3,5,10));
E=logspace(-1,1,10);

NN=3;

index_n=0;
for n1=N
    index_n=index_n+1;
    
    n2=n1;
    n3=n1;
    
    n_bar = n1+n2+n3;

    D1 = X(1:n1,:);                                                        % The first data owner
    Y1 = Y(1:n1,:);                                                        % The first data owner

    D2 = X(n1+1:n1+n2,:);                                                  % The second data owner
    Y2 = Y(n1+1:n1+n2,:);                                                  % The second data owner

    D3 = X(n1+n2+1:n1+n2+n3,:);                                            % The third data owner
    Y3 = Y(n1+n2+1:n1+n2+n3,:);                                            % The third data owner

    Dtrain = X(1:n1+n2+n3,:);                                              % All training data, i.e., D1, D2, D3
    Ytrain = Y(1:n1+n2+n3,:);                                              % All training data, i.e., Y1, Y2, Y3

    no_privacy_theta=inv((2/n_bar)*(Dtrain')*Dtrain ...
                     +sigma*eye(n_attributes))*(2/n_bar)*(Dtrain'*Ytrain); % ML model without privacy
    
    index_eps=0;
    for eps1=E
        index_eps=index_eps+1;
        
        eps2=eps1;
        eps3=eps1;
        
        rtrain_ave=0;
        for it=1:itmax
            % initialize theta
            theta    = rand(1,n_attributes);
            theta_1  = rand(1,n_attributes);
            theta_2  = rand(1,n_attributes);
            theta_3  = rand(1,n_attributes);

            for k = 1:T
                pk = 3e3*NN/T^2/sigma;

                switch ceil(3*rand(1,1))
                    case 1
                        bar_theta=(theta_1+theta)/2;

                        D1_sg = subgradient(D1,Y1,bar_theta); 
                        b1 = T*n_attributes/(n1*eps1);

                        w1 = b1*log(rand(1,1)/rand(1,1)); 

                        theta_1    = projectoperator(bar_theta ...
                                        -(sigma*bar_theta)*pk/(2*NN) ...
                                        -n1*(D1_sg+w1)*pk/n_bar,1e5);

                        theta      = projectoperator(bar_theta ...
                                        -(sigma*bar_theta)*pk/(NN*2),1e5);
                    case 2
                        bar_theta=(theta_2+theta)/2;

                        D2_sg = subgradient(D2,Y2,bar_theta);
                        b2 = T*n_attributes/(n2*eps2);

                        w2 = b2*log(rand(1,1)/rand(1,1)); 

                        theta_2    = projectoperator(bar_theta ...
                                        -(sigma*bar_theta)*pk/(2*NN) ...
                                        -n1*(D2_sg+w2)*pk/n_bar,1e5);

                        theta      = projectoperator(bar_theta ...
                                        -(sigma*bar_theta)*pk/(NN*2),1e5);               
                    case 3
                        bar_theta=(theta_3+theta)/2;

                        D3_sg = subgradient(D3,Y3,bar_theta);
                        b3 = T*n_attributes/(n3*eps3);

                        w3 = b3*log(rand(1,1)/rand(1,1));

                        theta_3    = projectoperator(bar_theta ...
                                        -(sigma*bar_theta)*pk/(2*NN) ...
                                        -n1*(D3_sg+w3)*pk/n_bar,1e5);                    

                        theta      = projectoperator(bar_theta ...
                                        -(sigma*bar_theta)*pk/(NN*2),1e5);
                end
            end
            rtrain_ave=rtrain_ave+((lossf(Dtrain,Ytrain,theta,sigma) ...
                      -lossf(Dtrain,Ytrain,no_privacy_theta',sigma)) ...
                      /lossf(Dtrain,Ytrain,no_privacy_theta',sigma))/itmax;
        end
        Z(index_n,index_eps)=rtrain_ave; 
    end
end

N=ceil(logspace(3,5,10));
E=logspace(-1,1,10);

index=0;
for i=1:length(N)
    for j=1:length(E)
        index=index+1;
        xdata(1,index)=(N(i));
        xdata(2,index)=(E(j));
        zdata(1,index)=(Z(i,j));
    end
end

scatter3(xdata(1,:),xdata(2,:),zdata,'ok','filled')

hold on

NN=ceil(logspace(3,5,50));
EE=logspace(-1,1,50);
for i=1:length(NN)
    for j=1:length(EE)
        ZZ(i,j)=1./(EE(j)^2*NN(i)^2);
    end
end
s=surf(NN,EE,2.1536e+09*ZZ,'EdgeColor',[1 1 1]/2,'FaceColor','none');

set(gca,'ZScale','log')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'ColorScale','log')
colormap(gray);
caxis([1e-3 1e6])
view(135,30);
zlim([1e-3 1e6])

zticks(10.^[-3 0 3 6])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize', 17)

figure('Renderer', 'painters', 'Position', [600 300 500 200])

loglog(E,Z(1,:) ,'sk','linewidth',1,'MarkerSize',10,'MarkerFaceColor','k')
hold on
loglog(E,Z(6,:) ,'dk','linewidth',1,'MarkerSize',10,'MarkerFaceColor','k')
loglog(E,Z(10,:),'ok','linewidth',1,'MarkerSize',10,'MarkerFaceColor','k')
loglog(E,10^(mean(log10(Z(1 ,:))'+2*log10(E)'))./E.^2,'-k','linewidth',2)
loglog(E,10^(mean(log10(Z(6 ,:))'+2*log10(E)'))./E.^2,'-k','linewidth',2)
loglog(E,10^(mean(log10(Z(10,:))'+2*log10(E)'))./E.^2,'-k','linewidth',2)

grid on

ylim([1e-3 1e9])
yticks(10.^(-3:3:9));

h=legend('$n=10^3$','$n=10^4$','$n=10^5$');
set(h,'fontsize',17,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize', 17)

figure('Renderer', 'painters', 'Position', [600 300 500 200])
loglog(N,Z(:,1) ,'sk','linewidth',1,'MarkerSize',10,'MarkerFaceColor','k')
hold on
loglog(N,Z(:,5) ,'dk','linewidth',1,'MarkerSize',10,'MarkerFaceColor','k')
loglog(N,Z(:,10),'ok','linewidth',1,'MarkerSize',10,'MarkerFaceColor','k')
loglog(N,10^(mean(log10(Z(:,1 ))+2*log10(N)'))./N.^2,'-k','linewidth',2)
loglog(N,10^(mean(log10(Z(:,5 ))+2*log10(N)'))./N.^2,'-k','linewidth',2)
loglog(N,10^(mean(log10(Z(:,10))+2*log10(N)'))./N.^2,'-k','linewidth',2)

grid on

ylim([1e-3 1e9])
yticks(10.^(-3:3:9));

h=legend('$\epsilon=0.1$','$\epsilon=1$','$\epsilon=10$');
set(h,'fontsize',17,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize', 17)

%% Generating Figure 7

clear theta

n = 10000;                                                                 % size of each dataset

NN=[2 4 8 16 32 64];                                                       % number of datasets
EE=logspace(-1,1,10);

index_N=0;
for N=NN
    index_N=index_N+1;
    
    n_bar = N*n;

    for i=1:N
        D_train{i}=X((i-1)*n+1:i*n,:);                                     % Each dataset
        Y_train{i}=Y((i-1)*n+1:i*n,:);                                     % Each dataset
    end
    
    Dtrain = X(1:n*N,:);                                                   % All training datasets
    Ytrain = Y(1:n*N,:);                                                   % All training datasets
    
    no_privacy_theta_1=inv((2/n)*D_train{1}'*D_train{1} ...
                  +sigma*eye(n_attributes))*(2/n)*(D_train{1}'*Y_train{1});% ML model without privacy based on 
                                                                           % data of the first data owner
    no_privacy_theta  =inv((2/n_bar)*(Dtrain')*Dtrain ...
                      +sigma*eye(n_attributes))*(2/n_bar)*(Dtrain'*Ytrain);% ML model without privacy based on all data
        
    index_e=0;
    for E=EE
        index_e=index_e+1;
        
        rtrain_ave=0;
        for it=1:itmax
            % initialize theta
            theta_0 = zeros(1,n_attributes);
            for i=1:N
                theta{i} = zeros(1,n_attributes); 
            end

            for k = 1:T
                pk = 3e3*N/T^2/sigma;

                i=ceil(N*rand(1,1));

                bar_theta=(theta{i}+theta_0)/2;

                D_sg = subgradient(D_train{i},Y_train{i},bar_theta); 
                b = T*n_attributes/(n*E);

                w = b*log(rand(1,1)/rand(1,1)); 

                theta{i}   = projectoperator(bar_theta ...
                                -(sigma*bar_theta)*pk/(2*N) ...
                                -n*(D_sg+w)*pk/n_bar       ,1e8);

                theta_0    = projectoperator(bar_theta ...
                                -(sigma*bar_theta)*pk/(2*N),1e8);

            end
            rtrain_ave=rtrain_ave+((lossf(Dtrain,Ytrain,theta_0,sigma) ...
                      -lossf(Dtrain,Ytrain,no_privacy_theta',sigma))   ...
                      /lossf(Dtrain,Ytrain,no_privacy_theta',sigma))/itmax;
        end
        U (index_N,index_e)=rtrain_ave;
        UU(index_N,index_e)=lossf(Dtrain,Ytrain,no_privacy_theta_1',sigma) ...
                            /lossf(Dtrain,Ytrain,no_privacy_theta',sigma)-1;
    end
end

NN=[2 4 8 16 32 64];
EE=logspace(-1,1,10);

surf(NN,EE,U','EdgeColor',[1 1 1]/2,'FaceColor','none')

hold on
surf(NN,EE,UU','EdgeColor','none','FaceColor',[1 1 1]/2)

index=0;
for i=1:length(NN)
    for j=1:length(EE)
        index=index+1;
        xdata(1,index)=(NN(i));
        xdata(2,index)=(EE(j));
        zdata(1,index)=(U(i,j));
    end
end

scatter3(xdata(1,:),xdata(2,:),zdata,'ok','filled')

view(150,20);
set(gca,'yscale','log')
set(gca,'xscale','log')
set(gca,'zscale','log')
xlim([1e0 1e2])
zlim([1e-2 1e4])

set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize', 17)