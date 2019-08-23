clc
clear

%#ok<*AGROW>
%#ok<*SAGROW>
%#ok<*MINV>

%% Reading converted data

XX=csvread('X.csv');
YY=csvread('Y.csv');

%% Sanitizing the data (removing Nan values, removing useless parts, etc)

XX(:,13)=[];

[I,J]=find(isnan(XX));

XX(I,:)=[];

Unique=unique(XX(:,1));                                                    % All hospitals

m=9;
[V,D]=eigs(XX(end-50000:end,[2:5 7:11])'*XX(end-50000:end,[2:5 7:11]),m);  % PCA based on the last 50000 entries
XX(:,[2:5 7:11])=XX(:,[2:5 7:11])*V*inv(sqrt(D));

index_data=[];                                                             % Selecting hospitals with more than 10000 records
for i=1:length(Unique)
    if sum(XX(:,1)==Unique(i))>=10000
        index_data=[index_data i]; 
    end
end

N=length(index_data);                                                      % Number of hospitals with more than 10000 records


for i=1:N                                                                  % Data of hospitals with more than 10000 records
    index=find(XX(:,1)==Unique(index_data(i)));
    D_train{i}=XX(index,[2:5 7:11]);
    Y_train{i}=XX(index,6);
    num(i)=length(D_train{i});
end

%% Generating Figure 11

sigma=1e-5;

n_attributes = size(D_train{1},2);

T=1e3;

itmax=100;                                                                 % Number of runs for generating statistics

EE=[0.1,1,10];

n_bar=0;
for i=1:N
    n_bar=n_bar+length(D_train{i});
end


Dtrain = XX(:,[2:5 7:11]);                                                 % All training data, i.e., D1, D2, D3
Ytrain = XX(:,6);                                                          % All training data, i.e., Y1, Y2, Y3

no_privacy_theta=inv((2/n_bar)*(Dtrain')*Dtrain ...
                 +sigma*eye(n_attributes))*(2/n_bar)*(Dtrain'*Ytrain);     % ML model without privacy
    
index_e=0;
for E=EE
    index_e=index_e+1;

    rtrain_ave=0;
    for it=1:itmax
        % initialize theta
        theta_0 = ones(1,n_attributes);
        for i=1:N
            theta{i} = ones(1,n_attributes); 
        end

        for k = 1:T
            pk = 1e4*N/T^2/sigma;

            i=ceil(N*rand(1,1));

            bar_theta=(theta{i}+theta_0)/2;

            D_sg = subgradient(D_train{i},Y_train{i},bar_theta); 
            b = T*n_attributes/(length(D_train{i})*E);

            w = b*log(rand(1,1)/rand(1,1)); 

            theta{i}   = projectoperator(bar_theta ...
                            -(sigma*bar_theta)*pk/(2*N) ...
                            -length(D_train{i})*(D_sg+w)*pk/n_bar,1e5);

            theta_0    = projectoperator(bar_theta ...
                            -(sigma*bar_theta)*pk/(2*N)           ,1e5);
        end
        rtrain_ave=rtrain_ave+((lossf(Dtrain,Ytrain,theta_0,sigma) ...
                  -lossf(Dtrain,Ytrain,no_privacy_theta',sigma)) ...
                  /lossf(Dtrain,Ytrain,no_privacy_theta',sigma))/itmax;
    end
    Z(index_e)=rtrain_ave;
end

for i=1:N
    no_privacy_theta_1=inv((2/length(D_train{i}))*D_train{i}'*D_train{i} ...
                  +sigma*eye(n_attributes))*(2/length(D_train{i}))*(D_train{i}'*Y_train{i});% ML model without privacy based on 
                                                                           % data of the first data owner
    L(i)=length(D_train{i});
    Q(i)=lossf(Dtrain,Ytrain,no_privacy_theta_1',sigma) ...
         /lossf(Dtrain,Ytrain,no_privacy_theta',sigma)-1;
end

loglog(L,Q,'ok','MarkerSize',7,'MarkerFaceColor','k')
hold on
loglog([1e4 1e5],[Z(1) Z(1)],'-k','linewidth',2)
loglog([1e4 1e5],[Z(2) Z(2)],'--k','linewidth',2)
loglog([1e4 1e5],[Z(3) Z(3)],'-.k','linewidth',2)
ylim([1e-3 1e3])
grid on

set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize', 17)