% 调度子程序(常规日前调度MP)
function MP_Dispatch(Hvf_SP1,Hvf_SP2,Hdf_SP1,Hdf_SP2,ImAve_SP1,ImAve_SP2,I60Ave_SP1,I60Ave_SP2,Im_SP1,Im_SP2,I60_SP1,I60_SP2)
%% 数据及参数导入
load data.mat
%% 变量定义
% 火电机组变量
u_G = binvar(T,G.N,'full');     % 启停机状态
u_rpr = binvar(T,G.N,'full');   % 常规调峰状态
u_dpr = binvar(T,G.N,'full');   % 不投油深度调峰状态
u_dpro = binvar(T,G.N,'full');  % 投油深度调峰状态
u_vpr = binvar(T,G.N,'full');   % 上爬坡分区状态，取0为1区，取1为2区，前提是开机状态
u_ieq = binvar(T,G.N,5,'full'); % 积分电量分段状态变量
P_G = sdpvar(T,G.N,'full');     % 火电机组有功出力
loss = sdpvar(T,G.N,'full');    % 火电机组单一时段转子损耗成本
I_60u = sdpvar(T,G.N,'full');   % 功率突增扰动下的火电机组60s频差积分电量
I_60d = sdpvar(T,G.N,'full');   % 功率突减扰动下的火电机组60s频差积分电量
I_mu = sdpvar(T,G.N,'full');    % 功率突增扰动下的火电机组最大频差积分电量
I_md = sdpvar(T,G.N,'full');    % 功率突减扰动下的火电机组最大频差积分电量
% 电池储能变量
u_B = binvar(T,B.N,'full');     % 电池储能状态，放电为1，充电为0
P_B = sdpvar(T,B.N,'full');     % 电池储能功率，放电为正，充电为负
E_B = sdpvar(T,B.N,'full');     % 电池储能剩余能量
% 风电场变量
P_W = sdpvar(T,N_W,'full');     % 风电场有功出力
% 系统变量
pha = sdpvar(T,N_nd,'full');    % 节点相位
P_tl = sdpvar(T,N_tl,'full');   % 线路有功功率(有名值，方向：小节点指向大节点)
%% 目标函数
C_coal = sum(sum(repmat(G.coal2,T,1).*(P_G.^2)+repmat(G.coal1,T,1).*P_G+repmat(G.coal0,T,1).*u_G)); % 火电煤耗成本
C_on = sum(sum(repmat(G.on,T-1,1).*(u_G(2:T,:)-u_G(1:T-1,:)+abs(u_G(2:T,:)-u_G(1:T-1,:)))/2));      % 火电开机成本
C_loss = sum(sum(loss));                                                                            % 火电深度调峰转子损耗成本
C_oil = sum(sum(Coil*repmat(G.oil,T,1).*u_dpro));                                                   % 火电深度调峰油耗成本
C_carbon = sum(sum(repmat(G.carbon1,T,1).*P_G+repmat(G.carbon0,T,1).*u_G))*Ccb;                     % 火电碳交易成本
C_wind = sum(sum(P_Wp-P_W))*penalty;                                                                % 弃风惩罚成本
obj = C_coal+C_on+C_loss+C_oil+C_carbon+C_wind;
%% 火电机组运行约束
constraints = [];
% 火电机组状态约束
constraints = [constraints,u_rpr+u_dpr+u_dpro==u_G];
% 火电机组转子损耗约束
constraints = [constraints,repmat(G.buy*G.Pn,T,1).*(G.loss1*P_G./repmat(G.Pn,T,1)+G.loss0)-(1-u_dpr-u_dpro)*M<=loss<=...
    repmat(G.buy*G.Pn,T,1).*(G.loss1*P_G./repmat(G.Pn,T,1)+G.loss0)+(1-u_dpr-u_dpro)*M];
constraints = [constraints,-(u_dpr+u_dpro)*M<=loss<=(u_dpr+u_dpro)*M];
% 火电机组出力约束(按调峰状态确定)
constraints = [constraints,-u_G*M<=P_G<=u_G*M];
constraints = [constraints,repmat(G.Prpr,T,1)-(1-u_rpr)*M<=P_G<=repmat(G.Pn,T,1)+(1-u_rpr)*M];
constraints = [constraints,repmat(G.Pdpr,T,1)-(1-u_dpr)*M<=P_G<=repmat(G.Prpr-m,T,1)+(1-u_dpr)*M];
constraints = [constraints,repmat(G.Pdpro,T,1)-(1-u_dpro)*M<=P_G<=repmat(G.Pdpr-m,T,1)+(1-u_dpro)*M];
% 火电机组出力约束(按阶段性爬坡确定)
constraints = [constraints,P_G>=repmat(G.vdot,T,1)-(2-u_G-u_vpr)*M];
constraints = [constraints,P_G<=repmat(G.vdot-m,T,1)+(1-u_G+u_vpr)*M];
% 火电机组阶段性爬坡约束
constraints = [constraints,P_G(2:T,:)<=repmat(G.Prpr,T-1,1)+(dt-repmat(G.Tdpr,T-1,1)-(repmat(G.Pdpr,T-1,1)-P_G(1:T-1,:))./...
    repmat(G.vdpro,T-1,1)).*repmat(G.vrpr,T-1,1)+(1-u_G(1:T-1,:)+u_vpr(1:T-1,:))*M];    % 上爬坡1区
constraints = [constraints,P_G(1:T-1,:)<=repmat(G.Prpr,T-1,1)+(dt-repmat(G.Tdpr,T-1,1)-(repmat(G.Pdpr,T-1,1)-P_G(2:T,:))./...
    repmat(G.vdpro,T-1,1)).*repmat(G.vrpr,T-1,1)+(1-u_G(2:T,:)+u_vpr(2:T,:))*M];        % 下爬坡1区
% 火电机组启停时间约束
for n = 1:G.N
    for t = 1:T-G.Tud(n)
        constraints = [constraints,sum(u_G(t+1:t+G.Tud(n),n))>=G.Tud(n)*(u_G(t+1,n)-u_G(t,n))];
        constraints = [constraints,sum(1-u_G(t+1:t+G.Tud(n),n))>=G.Tud(n)*(u_G(t,n)-u_G(t+1,n))];
    end
    for t = T-G.Tud(n)+1:T-1
        constraints = [constraints,sum(u_G(t+1:T,n))>=(T-t)*(u_G(t+1,n)-u_G(t,n))];
        constraints = [constraints,sum(1-u_G(t+1:T,n))>=(T-t)*(u_G(t,n)-u_G(t+1,n))];
    end
end
% 火电机组持续调峰约束
for n = 1:G.N
    for t = 1:T-G.Tpr
        constraints = [constraints,sum(u_rpr(t+1:t+G.Tpr,n))>=G.Tpr*(u_rpr(t+1,n)-u_rpr(t,n))];
        constraints = [constraints,sum(u_dpr(t+1:t+G.Tpr,n))>=G.Tpr*(u_dpr(t+1,n)-u_dpr(t,n))];
        constraints = [constraints,sum(u_dpro(t+1:t+G.Tpr,n))>=G.Tpr*(u_dpro(t+1,n)-u_dpro(t,n))];
    end
    for t = T-G.Tpr+1:T-1
        constraints = [constraints,sum(u_rpr(t+1:T,n))>=(T-t)*(u_rpr(t+1,n)-u_rpr(t,n))];
        constraints = [constraints,sum(u_dpr(t+1:T,n))>=(T-t)*(u_dpr(t+1,n)-u_dpr(t,n))];
        constraints = [constraints,sum(u_dpro(t+1:T,n))>=(T-t)*(u_dpro(t+1,n)-u_dpro(t,n))];
    end
end
% 火电机组启停-调峰衔接约束
constraints = [constraints,u_rpr(2:T,:)>=u_G(2:T,:)-u_G(1:T-1,:)];
%% 电池储能运行约束
% 电池储能功率约束
constraints = [constraints,-(1-u_B)*M<=P_B<=repmat(B.Pmax,T,1)+(1-u_B)*M];
constraints = [constraints,-repmat(B.Pmax,T,1)-u_B*M<=P_B<=u_B*M];
% 电池储能能量约束
constraints = [constraints,E_B(1,:)==B.Eini];
constraints = [constraints,repmat(B.Emin,T,1)<=E_B<=repmat(B.Emax,T,1)];
constraints = [constraints,P_B(1:T-1,:)/B.eff-(1-u_B(1:T-1,:))*M<=...
    E_B(1:T-1,:)-E_B(2:T,:)<=P_B(1:T-1,:)/B.eff+(1-u_B(1:T-1,:))*M];
constraints = [constraints,P_B(1:T-1,:)*B.eff-u_B(1:T-1,:)*M<=...
    E_B(1:T-1,:)-E_B(2:T,:)<=P_B(1:T-1,:)*B.eff+u_B(1:T-1,:)*M];
constraints = [constraints,B.Eini-(1-u_B(T,:))*M<=E_B(T,:)-P_B(T,:)/B.eff<=B.Emax+(1-u_B(T,:))*M];
constraints = [constraints,B.Eini-u_B(T,:)*M<=E_B(T,:)-P_B(T,:)*B.eff<=B.Emax+u_B(T,:)*M];
%% 系统运行约束
% 风电消纳约束
constraints = [constraints,0<=P_W<=P_Wp];
constraints = [constraints,sum(sum(P_W))/sum(sum(P_Wp))>=uti_Sys];
% 直流潮流约束(有名值)
for n = 1:N_tl
    constraints = [constraints,P_tl(:,n)==Ub^2*(pha(:,hn_tl(n))-pha(:,en_tl(n)))/X_tl(n)];
end
% 平衡节点(26节点)约束
constraints = [constraints,pha(:,26)==0];
% 线路功率约束
constraints = [constraints,-P_tlmax<=P_tl<=P_tlmax];
% 节点功率平衡约束(流入为正，流出为负，有名值)
for n = 1:N_nd
    -P_L*P_nd(n);           % 添加有功负荷
    if detp(n)=="N"
    elseif detp(n)=="G"
        ans+P_G(:,den(n));  % 添加火电机组有功出力
    elseif detp(n)=="W"
        ans+P_W(:,den(n));  % 添加风电场有功出力
    elseif detp(n)=="B"
        ans+P_B(:,den(n));  % 添加电池储能有功出力
    end
    for i = 1:length(tln{n})
        if tltp{n}(i)==1
            ans-P_tl(:,tln{n}(i));  % 添加线路有功功率(本节点为线路首节点)
        elseif tltp{n}(i)==0
            ans+P_tl(:,tln{n}(i));  % 添加线路有功功率(本节点为线路末节点)
        end
    end
    constraints = [constraints,ans==0];
end
% 系统备用约束
constraints = [constraints,sum(u_G.*repmat(G.Pn,T,1)-P_G,2)+sum(repmat(B.Pmax,T,1)-P_B,2)>=sum(P_Wer,2)+P_Ler];
constraints = [constraints,sum(u_G.*repmat(G.Pn,T,1)-P_G,2)+sum(B.eff*(E_B-repmat(B.Emin,T,1))-P_B,2)>=sum(P_Wer,2)+P_Ler];
constraints = [constraints,sum(P_G-u_G.*repmat(G.Pdpro,T,1),2)+sum(P_B+repmat(B.Pmax,T,1),2)>=sum(P_Wer,2)+P_Ler];
constraints = [constraints,sum(P_G-u_G.*repmat(G.Pdpro,T,1),2)+sum(P_B-(E_B-repmat(B.Emax,T,1))/B.eff,2)>=sum(P_Wer,2)+P_Ler];
%% SP反馈约束
% 火电机组开机数量约束
for t = 1:T
    if Hvf_SP1(t)>0
        constraints = [constraints,sum(u_G(t,:).*G.Pn.*G.H)>=Hvf_SP1(t)];
    end
    if Hvf_SP2(t)>0
        constraints = [constraints,sum(u_G(t,:).*G.Pn.*G.H)>=Hvf_SP2(t)];
    end
    if Hdf_SP1(t)>0
        constraints = [constraints,sum(u_G(t,:).*G.Pn.*G.H)>=Hdf_SP1(t)];
    end
    if Hdf_SP2(t)>0
        constraints = [constraints,sum(u_G(t,:).*G.Pn.*G.H)>=Hdf_SP2(t)];
    end
    if ImAve_SP1(t)>0
        constraints = [constraints,sum(u_G(t,:).*G.ImuAve)>=ImAve_SP1(t)];
    end
    if ImAve_SP2(t)>0
        constraints = [constraints,sum(u_G(t,:).*G.ImdAve)>=ImAve_SP2(t)];
    end
    if I60Ave_SP1(t)>0
        constraints = [constraints,sum(u_G(t,:).*G.I60uAve)>=I60Ave_SP1(t)];
    end
    if I60Ave_SP2(t)>0
        constraints = [constraints,sum(u_G(t,:).*G.I60dAve)>=I60Ave_SP2(t)];
    end
end
% 积分电量计算式约束
for t = 1:T
    for n = 1:G.N
        if Im_SP1(t)>0 || Im_SP2(t)>0 || I60_SP1(t)>0 || I60_SP2(t)>0
            % 状态约束
            constraints = [constraints,sum(u_ieq(t,n,:))==u_G(t,n)];
            % 出力约束
            for i = 1:5
                constraints = [constraints,G.dotI(n,i)*G.Pn(n)-(1-u_ieq(t,n,i))*M<=P_G(t,n)<=...
                    G.dotI(n,i+1)*G.Pn(n)+(1-u_ieq(t,n,i))*M];
            end
        end
        if Im_SP1(t)>0
            % 功率突增扰动下的最大频差积分电量约束
            constraints = [constraints,-u_G(t,n)*M<=I_mu(t,n)<=u_G(t,n)*M];
            for i = 1:4
                if i == 1
                    constraints = [constraints,G.kmu(n,i)/G.Pn(n)*P_G(t,n)+G.bmu(n,i)-(1-u_ieq(t,n,i)-u_ieq(t,n,i+1))*M<=...
                        I_mu(t,n)<=G.kmu(n,i)/G.Pn(n)*P_G(t,n)+G.bmu(n,i)+(1-u_ieq(t,n,i)-u_ieq(t,n,i+1))*M];
                else
                    constraints = [constraints,G.kmu(n,i)/G.Pn(n)*P_G(t,n)+G.bmu(n,i)-(1-u_ieq(t,n,i+1))*M<=...
                        I_mu(t,n)<=G.kmu(n,i)/G.Pn(n)*P_G(t,n)+G.bmu(n,i)+(1-u_ieq(t,n,i+1))*M];
                end
            end
        end
        if Im_SP2(t)>0
            % 功率突减扰动下的最大频差积分电量约束
            constraints = [constraints,-u_G(t,n)*M<=I_md(t,n)<=u_G(t,n)*M];
            for i = 1:4
                if i == 4
                    constraints = [constraints,G.kmd(n,i)/G.Pn(n)*P_G(t,n)+G.bmd(n,i)-(1-u_ieq(t,n,i)-u_ieq(t,n,i+1))*M<=...
                        I_md(t,n)<=G.kmd(n,i)/G.Pn(n)*P_G(t,n)+G.bmd(n,i)+(1-u_ieq(t,n,i)-u_ieq(t,n,i+1))*M];
                else
                    constraints = [constraints,G.kmd(n,i)/G.Pn(n)*P_G(t,n)+G.bmd(n,i)-(1-u_ieq(t,n,i))*M<=...
                        I_md(t,n)<=G.kmd(n,i)/G.Pn(n)*P_G(t,n)+G.bmd(n,i)+(1-u_ieq(t,n,i))*M];
                end
            end
        end
        if I60_SP1(t)>0
            % 功率突增扰动下的60s频差积分电量约束
            constraints = [constraints,-u_G(t,n)*M<=I_60u(t,n)<=u_G(t,n)*M];
            for i = 1:4
                if i == 1
                    constraints = [constraints,G.k60u(n,i)/G.Pn(n)*P_G(t,n)+G.b60u(n,i)-(1-u_ieq(t,n,i)-u_ieq(t,n,i+1))*M<=...
                        I_60u(t,n)<=G.k60u(n,i)/G.Pn(n)*P_G(t,n)+G.b60u(n,i)+(1-u_ieq(t,n,i)-u_ieq(t,n,i+1))*M];
                else
                    constraints = [constraints,G.k60u(n,i)/G.Pn(n)*P_G(t,n)+G.b60u(n,i)-(1-u_ieq(t,n,i+1))*M<=...
                        I_60u(t,n)<=G.k60u(n,i)/G.Pn(n)*P_G(t,n)+G.b60u(n,i)+(1-u_ieq(t,n,i+1))*M];
                end
            end
        end
        if I60_SP2(t)>0
            % 功率突减扰动下的60s频差积分电量约束
            constraints = [constraints,-u_G(t,n)*M<=I_60d(t,n)<=u_G(t,n)*M];
            for i = 1:4
                if i == 4
                    constraints = [constraints,G.k60d(n,i)/G.Pn(n)*P_G(t,n)+G.b60d(n,i)-(1-u_ieq(t,n,i)-u_ieq(t,n,i+1))*M<=...
                        I_60d(t,n)<=G.k60d(n,i)/G.Pn(n)*P_G(t,n)+G.b60d(n,i)+(1-u_ieq(t,n,i)-u_ieq(t,n,i+1))*M];
                else
                    constraints = [constraints,G.k60d(n,i)/G.Pn(n)*P_G(t,n)+G.b60d(n,i)-(1-u_ieq(t,n,i))*M<=...
                        I_60d(t,n)<=G.k60d(n,i)/G.Pn(n)*P_G(t,n)+G.b60d(n,i)+(1-u_ieq(t,n,i))*M];
                end
            end
        end
    end
end
% 积分电量边界约束
for t = 1:T
    if Im_SP1(t)>0
        constraints = [constraints,sum(I_mu(t,:))>=Im_SP1(t)];
    end
    if Im_SP2(t)>0
        constraints = [constraints,sum(I_md(t,:))>=Im_SP2(t)];
    end
    if I60_SP1(t)>0
        constraints = [constraints,sum(I_60u(t,:))>=I60_SP1(t)];
    end
    if I60_SP2(t)>0
        constraints = [constraints,sum(I_60d(t,:))>=I60_SP2(t)];
    end
end
%% 求解
ops = sdpsettings('solver','gurobi','verbose',1);
% ops.gurobi.TimeLimit = 200;
ops.gurobi.MIPGap = 0.001;
ops.gurobi.TuneTimeLimit = Inf;
solveinfo = optimize(constraints,obj,ops)
%% 结果保存
% 火电机组
u_G_MP = value(u_G);
u_rpr_MP = value(u_rpr);
u_dpr_MP = value(u_dpr);
u_dpro_MP = value(u_dpro);
u_vpr_MP = value(u_vpr);
P_G_MP = value(P_G);
% 电池储能
u_B_MP = value(u_B);
P_B_MP = value(P_B);
E_B_MP = value(E_B);
% 风电场
P_W_MP = value(P_W);
% 系统
pha_MP = value(pha);
P_tl_MP = value(P_tl);
% 成本
C_coal_MP = value(C_coal);
C_on_MP = value(C_on);
C_loss_MP = value(C_loss);
C_oil_MP = value(C_oil);
C_carbon_MP = value(C_carbon);
C_wind_MP = value(C_wind);
C_total = value(obj);
% 保存mat文件
save MP_result u_G_MP u_rpr_MP u_dpr_MP u_dpro_MP u_vpr_MP P_G_MP u_B_MP P_B_MP E_B_MP P_W_MP pha_MP ...
    P_tl_MP C_coal_MP C_on_MP C_loss_MP C_oil_MP C_carbon_MP C_wind_MP C_total solveinfo
end