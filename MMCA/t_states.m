clc
clear
%求解不同参数对密度的影响
N=2000;
lamda=0.3;
delta=0.6;
mu=0.8;
ter=40;
gama=1.5;
m=0.1;
stp=30;%
PUS_con = zeros(1,ter);         %存放每个时刻t的变量的密度
PAS_con = zeros(1,ter);
PAG_con= zeros(1,ter);
PUS_con(1) = 0.8;         %存放每个时刻t的变量的密度
PAS_con(1) =0.15;
PAG_con(1)= 0.05;
load ba 
load ER 
% 每次仿真的结果预分配
       parfor xun = 1:ter
              beta_R = xun / 80
            PUS=0.8*ones(1,N);%每个节点的初始状态的值
            PAG=0.1*ones(1,N);
            PAS=0.1*ones(1,N);

            PUS_UPDATE=zeros(1,N);
            PAS_UPDATE=zeros(1,N);
            PAG_UPDATE=zeros(1,N);
      
            rA=zeros(1,N);
            qA=zeros(1,N);
            qU=zeros(1,N);

            RA=zeros(N,N);
            QA=zeros(N,N);
            QU=zeros(N,N);
   
            for t = 1:stp
                for i =1:N
                    for j =1:N
                      RA(j,i)=1-A(i,j)*(PAG(1,j)+PAS(1,j))*lamda;
                      QA(j,i)=1-B(i,j)*PAG(1,j)*beta_R*gama;
                      QU(j,i)=1-B(i,j)*PAG(1,j)*beta_R;
                    end
                    tempprodRA=cumprod(RA(:,i));
                    rA(1,i)=tempprodRA(N);
                    tempprodQA=cumprod(QA(:,i));
                    qA(1,i)=tempprodQA(N);
                    tempprodQU=cumprod(QU(:,i));
                    qU(1,i)=tempprodQU(N);
 PUS_UPDATE(1,i)=PUS(1,i)* rA(1,i)*qU(1,i)*(1-m)+PAS(1,i)*delta*qU(1,i)*(1-m)+PAG(1,i)*delta*mu*(1-m);
 PAS_UPDATE(1,i)=PUS(1,i)*(1-rA(1,i))*qA(1,i)*(1-m)+PAS(1,i)*(1-delta)*qA(1,i)*(1-m)+PAG(1,i)*(1-delta)*mu*(1-m);  
 PAG_UPDATE(1,i)=PUS(1,i)*rA(1,i)*(1-qU(1,i)+qU(1,i)*m)+PUS(1,i)*(1-rA(1,i))*(1-qA(1,i)+qA(1,i)*m)+PAS(1,i)*(1-delta)*((1-qA(1,i))+qA(1,i)*m)...
                 +PAS(1,i)*delta*((1-qU(1,i))+qU(1,i)*m)+PAG(1,i)*(1-mu+mu*m);
                end
                PUS=PUS_UPDATE;
                PAG=PAG_UPDATE;
                PAS=PAS_UPDATE;
            end
           PAG_con(1,xun)=sum(PAG)/N; 
           PAS_con(1,xun)=sum(PAS)/N; 
           PUS_con(1,xun)=sum(PUS)/N;       
       end
 
xzhou=(1:ter)/80;
hold on;
box on;
grid off;
set(gca,'Fontsize',15);
set(gca, 'TickLabelInterpreter', 'latex');
plot(xzhou,PAG_con,'-o','color',[46/256 106/256 89/256],'MarkerFaceColor',[46/256 106/256 89/256]);
plot(xzhou,PAS_con,'-^','color',[176/256 76/256 72/256],'MarkerFaceColor',[176/256 76/256 72/256]);
plot(xzhou,PUS_con,'-v','color',[220/256 178/256 105/256],'MarkerFaceColor',[220/256 178/256 105/256]);
set(gcf,'DefaultTextInterpreter','latex');
xlabel('$\beta$','FontSize',15);ylabel('$proportion $','FontSize',15);   
h=legend({'$\rho^{AG}$', '$\rho^{AS}$', '$\rho^{US}$'}, 'Interpreter', 'latex', 'FontSize', 15);
set(h,'Interpreter','latex','FontSize',15)%,
save('PAG_con.mat','PAG_con')
save('PAS_con.mat','PAS_con')
save('PUS_con.mat','PUS_con') 
saveas(gcf, 'mmc.fig')