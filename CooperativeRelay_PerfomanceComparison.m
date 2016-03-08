%Cooperative Relay Communication
%Performance Comparison Between Optimal Scheme and the Random Scheme
clear all
T = 0.2e-3;     %duration of time slot
u = 20;         %Number of intervals of finite state space SNR received
avgSNR = 30;    %average SNR of Rayleigh fading channel = 30dB
W = 1e6;        %Bandwidth = 1MHz
t = 3e-6;       %Duation of obervation
M = 15:5:60;    %Number of SU Relay Nodes
lamda = avgSNR;
U = logspace(0,2,u);
Qu = zeros(1,20);
rk = zeros(1,20);
pk = zeros(1,20);

Qu_rk = zeros(20,3);

transData = zeros(length(M),4,100);
transData_Optimal = zeros(length(M),4,100);
for trial = 1:100
    %Initialization
    for L = 1:length(M)      %for all the number of nodes in the system

        gammaK = exprnd(lamda,1,M(L));
        gammaK_lin = 10.^(gammaK /10);    
        alpha = 0.5 + 0.5 .* rand(1,1);
        beta_low = (1-alpha)/alpha;
        Z = zeros(1,M(L));
        rk_N = zeros(M(L),5);
        initOrder = zeros(M(L),5);
        %Qu, Rk pk matrix calculation
        for s= 1:u
            if s == 1
                Qu(1) = 1 - exp(-U(s)/avgSNR);
                rk(1) = W *log2(1 + U(s));
                Qu_rk(1,1)= s;                
                Qu_rk(1,2)= Qu(1);
                Qu_rk(1,3)= rk(1);
            elseif s == 20
                Qu(20) = exp(-U(s)/avgSNR);
                rk(s) = W *log2(1 + U(s));
                Qu_rk(s,1)= s;                
                Qu_rk(s,2)= Qu(s);
                Qu_rk(s,3)= rk(s);
            else
                Qu(s) = exp(-U(s)/avgSNR) - exp(-U(s+1)/avgSNR);
                rk(s) = W *log2(1 + U(s));
                Qu_rk(s,1)= s;                
                Qu_rk(s,2)= Qu(s);
                Qu_rk(s,3)= rk(s);
            end 

        end     
      
        observationStep = 0;
        %% Intutive Order
         nodeSelect = 0;
         
        %for each number of links
         for node = 1:1:M(L)             
            theta = rand(1);      %probability of a node availability
            r = rand(1);
            beta = sum(r >= cumsum([0,theta])) - 1;                  
            H = beta > beta_low;        
            Rk = W * log2(1 + gammaK_lin(node));        
            Xk = Rk * H; 
            rk_N(node,1) = node;
            rk_N(node,2) = H;
            rk_N(node,3) = Rk;
            rk_N(node,4) = Xk;
            rk_N(node,5) = theta;
         
             for s=1:20
               pk(s) = Qu(s) * theta; 
             end
             leftSum = 0;
             rightSum = 0;
             if node == 1
                Cm = 1 - (M(L)*t)/T;
                Z(1) = Cm * rk*pk';
             else
                 Cmj = 1 - (M(L)-node)*t/T;
                 for k = 1:u
                        if (Cmj * rk(k)) > Z(node - 1)
                            leftSum = leftSum + Cmj * rk(k) * pk(k);
                        else
                            rightSum = rightSum + Z(node - 1) * pk(k);
                        end
                 end
                 Z(node) = leftSum + rightSum;
             end
         end
         
        initOrder = sortrows(rk_N,4);
        k_count = 0;
        for node = 1:+1:M(L)-1
            k_count = k_count +1;      
            Ck = 1 - (node*t)/T;
            Xk = initOrder(M(L) - k_count + 1,4);
            Yk = Ck * Xk;

            if Yk < Z(M(L)- node)
               continue;
            else
                %select the kth SU Node as cooperative relay
                nodeSelect = initOrder(M(L)-k_count+1,1);
                break;
            end      

        end
        %select the node as the SU node
        transData_Optimal(L,1,trial) = M(L);
        transData_Optimal(L,2,trial) = nodeSelect;
        transData_Optimal(L,3,trial) = (Xk*(T-(k_count*t)))/1000;
        transData_Optimal(L,4,trial) = k_count;
        
        
        %% Sequential Order
        Z = zeros(1,M(L));
        nodeSelect = 0;
        k_count = 0;
        %for each number of links
         for node = 1:+1:M(L)-1
            k_count = k_count +1;      
            Ck = 1 - (node*t)/T;
            Xk = rk_N(node,4);
            Yk = Ck * Xk;

            if Yk < Z(M(L)- node)
               continue;
            else
                %select the kth SU Node as cooperative relay
                nodeSelect = rk_N(node,1);
                break;
            end      

        end
        %select the node as the SU node
        transData(L,1,trial) = M(L);
        transData(L,2,trial) = nodeSelect;
        transData(L,3,trial) = (Xk*(T-k_count*t))/1000;
        transData(L,4,trial) = k_count;


    end
end
figure
data = mean(transData,3);
plot(M,data(:,3),'-sr');
hold on
xlim([15,60]);
xlabel('Number of SU Candidate Relays')
ylabel('Transmitted data in KBits')
grid on
data_optimal = mean(transData_Optimal,3);
plot(M,data_optimal(:,3),'-*g');


