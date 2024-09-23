function [vfini,dfend,dfpeak] = frequency(u_G,P_G,P_ld,P_D)
% u_G：火电机组运行状态
% P_G：火电机组出力有名值
% P_ld：额定频率下的系统负荷有名值
% P_D：功率扰动有名值
% vfini：初始频率变化率标幺值
% dfend：60s频差标幺值
% dfpeak：最大频差标幺值
%% 数据导入
load data.mat
%% 关键参数确定
pT = zeros(1,G.N);          % 主汽压力标幺值
CV = zeros(1,G.N);          % 阀门开度标幺值
Tch = zeros(1,G.N);         % 高压缸前汽室容积时间常数
Trh = zeros(1,G.N);         % 再热器容积时间常数
os = zeros(1,G.N);          % 高压缸功率自然过调系数
Fhp = zeros(1,G.N);         % 高压缸输出功率占比
pg = P_G./G.Pn;             % 火电出力标幺值
kG = 2*sum(u_G.*G.Pn.*G.H); % 系统惯性响应系数
for n = 1:G.N
    if pg(n)>=G.Pc1s(n)
        pT(n) = G.prc1;
        Tch(n) = G.Tchc1(n)/pg(n);
    elseif (pg(n)<G.Pc1s(n))&&(pg(n)>=G.Psc2(n))
        pT(n) = G.ks(n)*pg(n);
        Tch(n) = G.Tchs(n);
    elseif (pg(n)<G.Psc2(n))&&(P_G(n)>=G.Pdpro(n))
        pT(n) = G.prc2(n);
        Tch(n) = G.Tchc2(n)/pg(n);
    end
end
for n = 1:G.N
    if pg(n)>=G.Prh
        Trh(n) = G.Trhc(n)/pg(n);
    elseif (pg(n)<G.Prh)&&(P_G(n)>=G.Pdpro(n))
        Trh(n) = G.Trhs(n);
    end
end
CV = pg./pT;
os = G.osn.*(G.osk3*(pg.^3)+G.osk2*(pg.^2)+G.osk1*pg+G.osk0);
Fhp = G.Fhpn.*(G.Fhp4*(pg.^4)+G.Fhp3*(pg.^3)+G.Fhp2*(pg.^2)+G.Fhp1*pg+G.Fhp0);
%% 差分变量定义
df = zeros(N_ni,1);     % 系统频差(标幺值)
dfg = zeros(N_ni,G.N);  % 火电机组调频死区后频差(标幺值)
dCV = zeros(N_ni,G.N);  % 火电机组阀门开度变化量(标幺值)
dQ = zeros(N_ni,G.N);   % 火电机组主汽流量变化量(标幺值)
dpT = zeros(N_ni,G.N);  % 火电机组主汽压力变化量(标幺值)
dQsh = zeros(N_ni,G.N); % 火电机组过热器蒸汽流量变化量(标幺值)
dpD = zeros(N_ni,G.N);  % 火电机组汽包压力变化量(标幺值)
dPg = zeros(N_ni,G.N);  % 火电机组功率限幅前输出功率变化量(标幺值)
vPg = zeros(N_ni,G.N);  % 火电机组功率限幅前输出功率变化率(标幺值)
dP_G = zeros(N_ni,G.N); % 火电机组PFR功率变化量(有名值)
dfb = zeros(N_ni,B.N);  % 电池储能调频死区后频差(标幺值)
dPb = zeros(N_ni,B.N);  % 电池储能功率限幅前输出功率变化量(标幺值)
dP_B = zeros(N_ni,B.N); % 电池储能PFR功率变化量(有名值)
%% 数值积分
for t = 2:N_ni
    df(t) = dt_ni/kG*(sum(dP_G(t-1,:))+sum(dP_B(t-1,:))-P_D)+(1-dt_ni*kL*P_ld/kG)*df(t-1);
    for n = 1:G.N
        dCV(t,n) = (1-dt_ni/G.Tse(n))*dCV(t-1,n)-dt_ni/G.R(n)/G.Tse(n)*dfg(t-1,n);
        dpD(t,n) = dpD(t-1,n)-dt_ni/G.Td(n)*dQsh(t-1,n);
        dpT(t,n) = dpT(t-1,n)+dt_ni/G.Tsh(n)*(dQsh(t-1,n)-dQ(t-1,n));
        dPg(t,n) = dt_ni*vPg(t-1,n)+dPg(t-1,n);
        dfg(t,n) = max(u_G(n)*(df(t)-G.dz),0)*(df(t)>=0)+min(u_G(n)*(df(t)+G.dz),0)*(df(t)<0);
        dQ(t,n) = dCV(t,n)*dpT(t,n)+pT(n)*dCV(t,n)+CV(n)*dpT(t,n);
        dQsh(t,n) = G.ksh(n)*sqrt(abs(dpD(t,n)-dpT(t,n)))*sign(dpD(t,n)-dpT(t,n));
        dP_G(t,n) = min([G.Pn(n)*dPg(t,n),G.al*G.Pn(n),G.Pn(n)-P_G(n)])*(dPg(t,n)>=0)+...
                    max([G.Pn(n)*dPg(t,n),-G.al*G.Pn(n),G.Pdpro(n)-P_G(n)])*(dPg(t,n)<0);
        vPg(t,n) = Fhp(n)*(1+os(n))/Tch(n)*dQ(t,n)+(1-dt_ni*(Trh(n)+Tch(n))/Trh(n)/Tch(n))*vPg(t-1,n)+...
                   (dt_ni/Trh(n)/Tch(n)-Fhp(n)*(1+os(n))/Tch(n))*dQ(t-1,n)-dt_ni/Trh(n)/Tch(n)*dPg(t-1,n);
    end
    for n = 1:B.N
        dfb(t,n) = max(df(t)-B.dz,0)*(df(t)>=0)+min(df(t)+B.dz,0)*(df(t)<0);
        dPb(t,n) = (1-dt_ni/B.Tc(n))*dPb(t-1,n)+(B.kv(n)/B.Tc(n)-dt_ni*B.kd(n)/B.Tc(n))*dfb(t-1,n)-B.kv(n)/B.Tc(n)*dfb(t,n);
        dP_B(t,n) = min(B.Pn(n)*dPb(t,n),B.al*B.Pn(n))*(dPb(t,n)>=0)+max(B.Pn(n)*dPb(t,n),-B.al*B.Pn(n))*(dPb(t,n)<0);
    end
end
%% 结果保存
vfini = (df(2)-df(1))/dt_ni;
dfend = df(end);
dfpeak = max(abs(df))*sign(dfend);
save NIresult df dP_G dP_B
end