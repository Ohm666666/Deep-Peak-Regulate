% 数据导入
%% 火电机组
G.N = 8;                                                                    % 数量
G.Pn = readmatrix('parameters.xlsx','Sheet','unit','Range','F3:M3');        % 额定功率
G.Prpr = readmatrix('parameters.xlsx','Sheet','unit','Range','F4:M4');      % RPR出力下限
G.Pdpr = readmatrix('parameters.xlsx','Sheet','unit','Range','F5:M5');      % DPR出力下限
G.Pdpro = readmatrix('parameters.xlsx','Sheet','unit','Range','F6:M6');     % DPRO出力下限
G.vrpr = readmatrix('parameters.xlsx','Sheet','unit','Range','F7:M7');      % RPR爬坡率
G.vdpr = readmatrix('parameters.xlsx','Sheet','unit','Range','F8:M8');      % DPR爬坡率
G.vdpro = readmatrix('parameters.xlsx','Sheet','unit','Range','F9:M9');     % DPRO爬坡率
G.coal2 = readmatrix('parameters.xlsx','Sheet','unit','Range','F10:M10');   % 煤耗成本二次项系数
G.coal1 = readmatrix('parameters.xlsx','Sheet','unit','Range','F11:M11');   % 煤耗成本一次项系数
G.coal0 = readmatrix('parameters.xlsx','Sheet','unit','Range','F12:M12');   % 煤耗成本常数项
G.carbon1 = readmatrix('parameters.xlsx','Sheet','unit','Range','F13:M13'); % 碳排放特性一次项系数
G.carbon0 = readmatrix('parameters.xlsx','Sheet','unit','Range','F14:M14'); % 碳排放特性常数项
G.Pc1s = readmatrix('parameters.xlsx','Sheet','unit','Range','F15:M15');    % 定压1段到滑压段的临界功率标幺值
G.prc1 = readmatrix('parameters.xlsx','Sheet','unit','Range','F16:F16');    % 定压1段主汽压力标幺值(所有机组相同)
G.Psc2 = readmatrix('parameters.xlsx','Sheet','unit','Range','F17:M17');    % 滑压段到定压2段的临界功率标幺值
G.ks = readmatrix('parameters.xlsx','Sheet','unit','Range','F18:M18');      % 滑压段斜率
G.prc2 = readmatrix('parameters.xlsx','Sheet','unit','Range','F19:M19');    % 定压2段主汽压力标幺值
G.Prh = readmatrix('parameters.xlsx','Sheet','unit','Range','F20:F20');     % 再热定压段到滑压段的临界功率标幺值(所有机组相同)
G.on = readmatrix('parameters.xlsx','Sheet','unit','Range','F21:M21');      % 单次开机成本
G.Tud = readmatrix('parameters.xlsx','Sheet','unit','Range','F22:M22');     % 启停机最小持续时间
G.Tpr = readmatrix('parameters.xlsx','Sheet','unit','Range','F23:F23');     % 各调峰状态最小持续时间(所有机组相同)
G.Tchc1 = readmatrix('parameters.xlsx','Sheet','unit','Range','F26:M26');   % 定压1段高压缸前汽室容积时间常数反比例系数
G.Tchs = readmatrix('parameters.xlsx','Sheet','unit','Range','F27:M27');    % 滑压段高压缸前汽室容积时间常数值
G.Tchc2 = readmatrix('parameters.xlsx','Sheet','unit','Range','F28:M28');   % 定压2段高压缸前汽室容积时间常数反比例系数
G.Trhc = readmatrix('parameters.xlsx','Sheet','unit','Range','F29:M29');    % 再热定压段再热器容积时间常数反比例系数
G.Trhs = readmatrix('parameters.xlsx','Sheet','unit','Range','F30:M30');    % 再热滑压段再热器容积时间常数值
G.osn = readmatrix('parameters.xlsx','Sheet','unit','Range','F31:M31');     % 高压缸功率自然过调系数额定值
G.osk3 = readmatrix('parameters.xlsx','Sheet','unit','Range','F32:F32');    % 高压缸功率自然过调系数三次项系数(所有机组相同)
G.osk2 = readmatrix('parameters.xlsx','Sheet','unit','Range','F33:F33');    % 高压缸功率自然过调系数二次项系数(所有机组相同)
G.osk1 = readmatrix('parameters.xlsx','Sheet','unit','Range','F34:F34');    % 高压缸功率自然过调系数一次项系数(所有机组相同)
G.osk0 = readmatrix('parameters.xlsx','Sheet','unit','Range','F35:F35');    % 高压缸功率自然过调系数常数项(所有机组相同)
G.Fhpn = readmatrix('parameters.xlsx','Sheet','unit','Range','F36:M36');    % 高压缸输出功率占比额定值
G.Fhp4 = readmatrix('parameters.xlsx','Sheet','unit','Range','F37:F37');    % 高压缸输出功率占比四次项系数(所有机组相同)
G.Fhp3 = readmatrix('parameters.xlsx','Sheet','unit','Range','F38:F38');    % 高压缸输出功率占比三次项系数(所有机组相同)
G.Fhp2 = readmatrix('parameters.xlsx','Sheet','unit','Range','F39:F39');    % 高压缸输出功率占比二次项系数(所有机组相同)
G.Fhp1 = readmatrix('parameters.xlsx','Sheet','unit','Range','F40:F40');    % 高压缸输出功率占比一次项系数(所有机组相同)
G.Fhp0 = readmatrix('parameters.xlsx','Sheet','unit','Range','F41:F41');    % 高压缸输出功率占比常数项(所有机组相同)
G.R = readmatrix('parameters.xlsx','Sheet','unit','Range','F42:M42');       % 调差系数
G.Tse = readmatrix('parameters.xlsx','Sheet','unit','Range','F43:M43');     % 伺服机构时间常数
G.Td = readmatrix('parameters.xlsx','Sheet','unit','Range','F44:M44');      % 汽包容积时间常数
G.Tsh = readmatrix('parameters.xlsx','Sheet','unit','Range','F45:M45');     % 过热器容积时间常数
G.ksh = readmatrix('parameters.xlsx','Sheet','unit','Range','F46:M46');     % 过热器及主汽管道流量系数
G.al = readmatrix('parameters.xlsx','Sheet','unit','Range','F47:F47');      % 功率限幅比例(所有机组相同)
G.dz = readmatrix('parameters.xlsx','Sheet','unit','Range','F48:F48');      % 调频死区标幺值(所有机组相同)
G.H = readmatrix('parameters.xlsx','Sheet','unit','Range','F49:M49');       % 惯性时间常数
G.buy = readmatrix('parameters.xlsx','Sheet','unit','Range','F50:F50');     % 机组单位造价(所有机组相同)
G.loss1 = readmatrix('parameters.xlsx','Sheet','unit','Range','F51:F51');   % 寿命损耗率一次项系数(所有机组相同)
G.loss0 = readmatrix('parameters.xlsx','Sheet','unit','Range','F52:F52');   % 寿命损耗率常数项(所有机组相同)
G.oil = readmatrix('parameters.xlsx','Sheet','unit','Range','F53:M53');     % DPRO燃料油消耗速度
G.Tdpr = readmatrix('parameters.xlsx','Sheet','unit','Range','F55:M55');    % DPR区域爬坡用时
G.vdot = readmatrix('parameters.xlsx','Sheet','unit','Range','F57:M57');    % 功率分区点有名值
G.dotu = readmatrix('parameters.xlsx','Sheet','IEQ_dfm','Range','B11:F18'); % 功率突增扰动下的积分电量分段点标幺值
G.dotd = readmatrix('parameters.xlsx','Sheet','IEQ_dfm','Range','B2:F9');   % 功率突减扰动下的积分电量分段点标幺值
G.dotI = [G.dotd(:,1:4) G.dotu(:,4:5)];                                     % 积分电量分段点标幺值
G.Imu = readmatrix('parameters.xlsx','Sheet','IEQ_dfm','Range','G11:K18');  % 功率突增扰动下的分段点处的最大频差积分电量
G.Imd = readmatrix('parameters.xlsx','Sheet','IEQ_dfm','Range','G2:K9');    % 功率突减扰动下的分段点处的最大频差积分电量
G.kmu = readmatrix('parameters.xlsx','Sheet','IEQ_dfm','Range','L11:O18');  % 功率突增扰动下的最大频差积分电量斜率
G.kmd = readmatrix('parameters.xlsx','Sheet','IEQ_dfm','Range','L2:O9');    % 功率突减扰动下的最大频差积分电量斜率
G.bmu = readmatrix('parameters.xlsx','Sheet','IEQ_dfm','Range','P11:S18');  % 功率突增扰动下的最大频差积分电量截距
G.bmd = readmatrix('parameters.xlsx','Sheet','IEQ_dfm','Range','P2:S9');    % 功率突减扰动下的最大频差积分电量截距
G.I60u = readmatrix('parameters.xlsx','Sheet','IEQ_df60','Range','G11:K18');% 功率突增扰动下的分段点处的60s频差积分电量
G.I60d = readmatrix('parameters.xlsx','Sheet','IEQ_df60','Range','G2:K9');  % 功率突减扰动下的分段点处的60s频差积分电量
G.k60u = readmatrix('parameters.xlsx','Sheet','IEQ_df60','Range','L11:O18');% 功率突增扰动下的60s频差积分电量斜率
G.k60d = readmatrix('parameters.xlsx','Sheet','IEQ_df60','Range','L2:O9');  % 功率突减扰动下的60s频差积分电量斜率
G.b60u = readmatrix('parameters.xlsx','Sheet','IEQ_df60','Range','P11:S18');% 功率突增扰动下的60s频差积分电量截距
G.b60d = readmatrix('parameters.xlsx','Sheet','IEQ_df60','Range','P2:S9');  % 功率突减扰动下的60s频差积分电量截距
G.ImuAve = zeros(1,G.N);                                                    % 功率突增扰动下的最大频差平均积分电量
G.ImdAve = zeros(1,G.N);                                                    % 功率突减扰动下的最大频差平均积分电量
G.I60uAve = zeros(1,G.N);                                                   % 功率突增扰动下的60s频差平均积分电量
G.I60dAve = zeros(1,G.N);                                                   % 功率突减扰动下的60s频差平均积分电量
for n = 1:G.N
    %% 计算G.ImuAve
    I = 0;
    for i = 1:10000
        P = 0.3+i/10000*0.7;
        if P>=G.dotu(n,1) && P<=G.dotu(n,2)
            I = I+G.kmu(n,1)*P+G.bmu(n,1);
        elseif P>G.dotu(n,2) && P<=G.dotu(n,3)
            I = I+G.kmu(n,2)*P+G.bmu(n,2);
        elseif P>G.dotu(n,3) && P<=G.dotu(n,4)
            I = I+G.kmu(n,3)*P+G.bmu(n,3);
        elseif P>G.dotu(n,4) && P<=G.dotu(n,5)
            I = I+G.kmu(n,4)*P+G.bmu(n,4);
        end
    end
    G.ImuAve(n) = I/(10000+1);
    %% 计算G.ImdAve
    I = 0;
    for i = 1:10000
        P = 0.3+i/10000*0.7;
        if P>=G.dotd(n,1) && P<=G.dotd(n,2)
            I = I+G.kmd(n,1)*P+G.bmd(n,1);
        elseif P>G.dotd(n,2) && P<=G.dotd(n,3)
            I = I+G.kmd(n,2)*P+G.bmd(n,2);
        elseif P>G.dotd(n,3) && P<=G.dotd(n,4)
            I = I+G.kmd(n,3)*P+G.bmd(n,3);
        elseif P>G.dotd(n,4) && P<=G.dotd(n,5)
            I = I+G.kmd(n,4)*P+G.bmd(n,4);
        end
    end
    G.ImdAve(n) = I/(10000+1);
    %% 计算G.I60uAve
    I = 0;
    for i = 1:10000
        P = 0.3+i/10000*0.7;
        if P>=G.dotu(n,1) && P<=G.dotu(n,2)
            I = I+G.k60u(n,1)*P+G.b60u(n,1);
        elseif P>G.dotu(n,2) && P<=G.dotu(n,3)
            I = I+G.k60u(n,2)*P+G.b60u(n,2);
        elseif P>G.dotu(n,3) && P<=G.dotu(n,4)
            I = I+G.k60u(n,3)*P+G.b60u(n,3);
        elseif P>G.dotu(n,4) && P<=G.dotu(n,5)
            I = I+G.k60u(n,4)*P+G.b60u(n,4);
        end
    end
    G.I60uAve(n) = I/(10000+1);
    %% 计算G.I60dAve
    I = 0;
    for i = 1:10000
        P = 0.3+i/10000*0.7;
        if P>=G.dotd(n,1) && P<=G.dotd(n,2)
            I = I+G.k60d(n,1)*P+G.b60d(n,1);
        elseif P>G.dotd(n,2) && P<=G.dotd(n,3)
            I = I+G.k60d(n,2)*P+G.b60d(n,2);
        elseif P>G.dotd(n,3) && P<=G.dotd(n,4)
            I = I+G.k60d(n,3)*P+G.b60d(n,3);
        elseif P>G.dotd(n,4) && P<=G.dotd(n,5)
            I = I+G.k60d(n,4)*P+G.b60d(n,4);
        end
    end
    G.I60dAve(n) = I/(10000+1);
end
%% 电池储能
B.N = 2;                                                                    % 数量
B.Pn = readmatrix('parameters.xlsx','Sheet','unit','Range','P7:Q7');        % 额定功率
B.Pmax = readmatrix('parameters.xlsx','Sheet','unit','Range','P9:Q9');      % 功率上限
B.Emax = readmatrix('parameters.xlsx','Sheet','unit','Range','P11:Q11');    % 能量上限
B.Emin = readmatrix('parameters.xlsx','Sheet','unit','Range','P12:Q12');    % 能量下限
B.Eini = readmatrix('parameters.xlsx','Sheet','unit','Range','P13:Q13');    % 初始能量
B.eff = readmatrix('parameters.xlsx','Sheet','unit','Range','P14:P14');     % 效率(所有机组相同)
B.kv = readmatrix('parameters.xlsx','Sheet','unit','Range','P15:Q15');      % 虚拟惯量系数
B.kd = readmatrix('parameters.xlsx','Sheet','unit','Range','P16:Q16');      % 频率下垂系数
B.Tc = readmatrix('parameters.xlsx','Sheet','unit','Range','P17:Q17');      % 变流器时间常数
B.dz = readmatrix('parameters.xlsx','Sheet','unit','Range','P18:P18');      % 调频死区标幺值(所有机组相同)
B.al = readmatrix('parameters.xlsx','Sheet','unit','Range','P19:P19');      % 功率限幅比例(所有机组相同)
%% 风电及负荷预测值
N_W = 2;                                                                    % 风电场数量
P_Wp = readmatrix('parameters.xlsx','Sheet','wind','Range','B3:C26');       % 风电预测功率
P_Wer = readmatrix('parameters.xlsx','Sheet','wind','Range','D3:E26');      % 风电预测误差
P_L = readmatrix('parameters.xlsx','Sheet','load','Range','B2:B25');        % 负荷预测功率
P_Ler = readmatrix('parameters.xlsx','Sheet','load','Range','C2:C25');      % 负荷预测误差
%% 线路参数
N_tl = 46;                                                                  % 线路数量
hn_tl = readmatrix('parameters.xlsx','Sheet','line','Range','B2:B47');      % 线路首节点编号
en_tl = readmatrix('parameters.xlsx','Sheet','line','Range','C2:C47');      % 线路末节点编号
X_tl = readmatrix('parameters.xlsx','Sheet','line','Range','D2:D47');       % 线路电抗(有名值)
%% 节点参数
N_nd = 39;                                                                  % 节点数量
P_nd = readmatrix('parameters.xlsx','Sheet','node','Range','B3:B41');       % 节点有功负荷比例
detp = readmatrix('parameters.xlsx','Sheet','node','Range','C3:C41','OutputType','string');  % 节点所接设备类型
den = readmatrix('parameters.xlsx','Sheet','node','Range','D3:D41');        % 节点所接设备编号
tln = cell(N_nd,1);                                                         % 节点相连线路编号
tltp = cell(N_nd,1);                                                        % 节点相连线路类型，节点为线路首节点时取1，为线路末节点时取0
for i = 1:N_nd
    for n = 1:N_tl
        if hn_tl(n)==i
            tln{i} = [tln{i},n];
            tltp{i} = [tltp{i},1];
        elseif en_tl(n)==i
            tln{i} = [tln{i},n];
            tltp{i} = [tltp{i},0];
        end
    end
end
%% 其他参数
dt = readmatrix('parameters.xlsx','Sheet','other','Range','B2:B2');         % 单位调度时长
T = readmatrix('parameters.xlsx','Sheet','other','Range','B3:B3');          % 日前调度时段数
Coil = readmatrix('parameters.xlsx','Sheet','other','Range','B4:B4');       % 油价
kL = readmatrix('parameters.xlsx','Sheet','other','Range','B5:B5');         % 负荷频率响应系数
fn = readmatrix('parameters.xlsx','Sheet','other','Range','B6:B6');         % 额定频率
vfmax = readmatrix('parameters.xlsx','Sheet','other','Range','B7:B7');      % 初始频率变化率限值标幺值
df60s = readmatrix('parameters.xlsx','Sheet','other','Range','B8:B8');      % 60s频差限值标幺值
dfmax = readmatrix('parameters.xlsx','Sheet','other','Range','B9:B9');      % 最大频差限值标幺值
Ccb = readmatrix('parameters.xlsx','Sheet','other','Range','B10:B10');      % 碳交易价格
dt_ni = readmatrix('parameters.xlsx','Sheet','other','Range','B11:B11');    % 数值积分步长
T_ni = readmatrix('parameters.xlsx','Sheet','other','Range','B12:B12');     % 数值积分总时长
N_ni = readmatrix('parameters.xlsx','Sheet','other','Range','B13:B13');     % 数值积分时步数
P_tlmax = readmatrix('parameters.xlsx','Sheet','other','Range','B14:B14');  % 线路容量上限
Ub = readmatrix('parameters.xlsx','Sheet','other','Range','B15:B15');       % 电压等级
penalty = readmatrix('parameters.xlsx','Sheet','other','Range','B16:B16');  % 弃风惩罚成本系数
uti_Sys = readmatrix('parameters.xlsx','Sheet','other','Range','B17:B17');  % 系统消纳率下限
step_W = readmatrix('parameters.xlsx','Sheet','other','Range','B18:B18');   % 风电阶跃扰动比例
step_L = readmatrix('parameters.xlsx','Sheet','other','Range','B19:B19');   % 负荷阶跃扰动比例
step_I = readmatrix('parameters.xlsx','Sheet','other','Range','B20:B20');   % 积分电量迭代调整比例
M = 1000000;                                                                % 足够大的正数
m = 0.00001;                                                                % 足够小的正数
%% 结果输出
clear i n I P
save data.mat