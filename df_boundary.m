function [flag_ol,IEQ] = df_boundary(u_G_MP,P_G_MP,t,flag_SP,flag_df)
% 当第n次迭代MP求解结果在SP校核中发生频差指标(60s频差、最大频差)越限时，计算第n+1次迭代MP优化割约束的边界
% u_G_MP：火电机组运行状态
% P_G_MP：火电机组出力有名值
% t：频差指标越限时段
% flag_SP：取1为SP1，取2为SP2
% flag_df：取1为60s频差，取2为最大频差
% flag_ol：t时段机组组合条件下，系统能否通过调整火电机组出力避免频差指标越限，取1为无法避免，取0则相反
% IEQ：flag_ol为0时，使得频差指标不越限的系统最小积分电量总和；flag_ol为1时取0
%% 数据导入
load data.mat
%% 确定积分电量参数及功率扰动值
flag = flag_SP*10+flag_df;
switch flag
    case 11     % SP1(功率突增扰动)，60s频差
        I = G.I60u;
        d = G.dotu;
        k = G.k60u;
        b = G.b60u;
        P_D = step_W*sum(P_Wp(t,:))+step_L*P_L(t);
    case 12     % SP1(功率突增扰动)，最大频差
        I = G.Imu;
        d = G.dotu;
        k = G.kmu;
        b = G.bmu;
        P_D = step_W*sum(P_Wp(t,:))+step_L*P_L(t);
    case 21     % SP2(功率突减扰动)，60s频差
        I = G.I60d;
        d = G.dotd;
        k = G.k60d;
        b = G.b60d;
        P_D = -step_W*sum(P_Wp(t,:))-step_L*P_L(t);
    case 22     % SP2(功率突减扰动)，最大频差
        I = G.Imd;
        d = G.dotd;
        k = G.kmd;
        b = G.bmd;
        P_D = -step_W*sum(P_Wp(t,:))-step_L*P_L(t);
end
%% 利用最大积分电量计算flag_ol与IEQ
[~,idx] = max(I,[],2);
for n = 1:G.N
    Pmax(1,n) = d(n,idx(n));  % 各机组取最大积分电量时，对应的出力标幺值
end
[~,dfend,dfpeak] = frequency(u_G_MP,Pmax.*G.Pn,P_L(t),P_D);
if flag_df == 1
    if abs(dfend)>df60s
        flag_ol = 1;
        IEQ = 0;
        return
    end
elseif flag_df == 2
    if abs(dfpeak)>dfmax
        flag_ol = 1;
        IEQ = 0;
        return
    end
end
%% 计算flag_ol为0时IEQ的值
% 计算I_set初始值，I_G_MP为各机组在P_G_MP出力条件下的积分电量
for n = 1:G.N
    if P_G_MP(n)/G.Pn(n)<=d(n,2)
        I_G_MP(n) = u_G_MP(n)*(k(n,1)*P_G_MP(n)/G.Pn(n)+b(n,1));
    elseif P_G_MP(n)/G.Pn(n)>d(n,2) && P_G_MP(n)/G.Pn(n)<=d(n,3)
        I_G_MP(n) = u_G_MP(n)*(k(n,2)*P_G_MP(n)/G.Pn(n)+b(n,2));
    elseif P_G_MP(n)/G.Pn(n)>d(n,3) && P_G_MP(n)/G.Pn(n)<=d(n,4)
        I_G_MP(n) = u_G_MP(n)*(k(n,3)*P_G_MP(n)/G.Pn(n)+b(n,3));
    elseif P_G_MP(n)/G.Pn(n)>d(n,4)
        I_G_MP(n) = u_G_MP(n)*(k(n,4)*P_G_MP(n)/G.Pn(n)+b(n,4));
    end
end
I_set = sum(I_G_MP);
% 优化求解IEQ
while 1
    I_set = I_set+step_I*sum(max(I,[],2)'.*u_G_MP); % 更新I_set
    if I_set>sum(max(I,[],2)'.*u_G_MP)              % 如果循环求解IEQ发现，只有当I_set取各机组的最大积分电量组合时，
        flag_ol = 1;                                % 才能保证频率指标不越限，则认为flag_ol=1。
        IEQ = 0;
        return
    end
    % 变量定义
    P_G = sdpvar(1,G.N,'full');     % 火电机组出力
    u_rpr = binvar(1,G.N,'full');   % 常规调峰状态
    u_dpr = binvar(1,G.N,'full');   % 不投油深度调峰状态
    u_dpro = binvar(1,G.N,'full');  % 投油深度调峰状态
    u_ieq = binvar(4,G.N,'full');   % 积分电量四分段直线状态
    I_G = sdpvar(1,G.N,'full');     % 火电机组积分电量
    loss = sdpvar(1,G.N,'full');    % 各机组转子损耗成本
    % 目标函数
    C_coal = sum(G.coal2.*(P_G.^2)+G.coal1.*P_G+G.coal0.*u_G_MP);   % 煤耗成本
    C_carbon = sum(G.carbon1.*P_G+G.carbon0.*u_G_MP)*Ccb;           % 碳交易成本
    C_loss = sum(loss);                                             % 转子损耗成本
    C_oil = sum(G.oil.*u_dpro)*Coil;                                % 投油成本
    obj = C_coal+C_carbon+C_loss+C_oil;                             % 总成本
    % 约束条件
    constraints = [];
    constraints = [constraints,u_rpr+u_dpr+u_dpro==u_G_MP];                         % 调峰状态约束
    constraints = [constraints,-u_G_MP*M<=P_G<=u_G_MP*M];                           % 火电机组出力约束
    constraints = [constraints,G.Prpr-(1-u_rpr)*M<=P_G<=G.Pn+(1-u_rpr)*M];          % 常规调峰出力约束
    constraints = [constraints,G.Pdpr-(1-u_dpr)*M<=P_G<=G.Prpr-m+(1-u_dpr)*M];      % 不投油深度调峰出力约束
    constraints = [constraints,G.Pdpro-(1-u_dpro)*M<=P_G<=G.Pdpr-m+(1-u_dpro)*M];   % 投油深度调峰出力约束
    constraints = [constraints,G.buy*G.Pn.*(G.loss1*P_G./G.Pn+G.loss0)-(1-u_dpr-u_dpro)*M<=loss<=...
        G.buy*G.Pn.*(G.loss1*P_G./G.Pn+G.loss0)+(1-u_dpr-u_dpro)*M];                % 转子损耗约束
    constraints = [constraints,-(u_dpr+u_dpro)*M<=loss<=(u_dpr+u_dpro)*M];          % 转子损耗约束
    constraints = [constraints,sum(u_ieq)==u_G_MP];                                 % 积分电量状态约束
    for i = 1:4                                                                     % 按积分电量分段的出力约束
        constraints = [constraints,d(:,i)'.*G.Pn-(1-u_ieq(i,:))*M<=P_G<=...
            d(:,i+1)'.*G.Pn+(1-u_ieq(i,:))*M];
    end
    constraints = [constraints,-u_G_MP*M<=I_G<=u_G_MP*M];                           % 积分电量约束
    for i = 1:4                                                                     % 积分电量约束
        constraints = [constraints,k(:,i)'.*P_G./G.Pn+b(:,i)'-(1-u_ieq(i,:))*M<=I_G<=...
        k(:,i)'.*P_G./G.Pn+b(:,i)'+(1-u_ieq(i,:))*M];
    end
    constraints = [constraints,sum(I_G)==I_set];                                    % 积分电量总和等于循环设定值
    % Gurobi求解
    ops = sdpsettings('solver','gurobi','verbose',0);
%     ops.gurobi.MIPGap = 0.001;
    ops.gurobi.TuneTimeLimit = Inf;
    optimize(constraints,obj,ops);
    % 判断频差指标是否越限
    [~,dfend,dfpeak] = frequency(u_G_MP,value(P_G),P_L(t),P_D);
    if flag_df == 1
        if abs(dfend)>df60s
            continue        % 若越限则继续循环
        else
            flag_ol = 0;
            IEQ = I_set;    % 若未越限，则输出IEQ
            return
        end
    elseif flag_df == 2
        if abs(dfpeak)>dfmax
            continue        % 若越限则继续循环
        else
            flag_ol = 0;
            IEQ = I_set;    % 若未越限，则输出IEQ
            return
        end
    end
end
end