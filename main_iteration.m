% 迭代控制主程序(频率越限校验SP，SP1为功率突增扰动，SP2为功率突减扰动)
clear;clc;tic
%% 数据及参数
Data
Its = 0;    % 迭代次数
icm_H = min(abs(G.Pn.*G.H-(G.Pn.*G.H)')+diag(NaN*ones(1,G.N)),[],'all','omitnan');          % 惯量总和边界迭代增量
icm_ImAve_SP1 = min(abs(G.ImuAve-G.ImuAve')+diag(NaN*ones(1,G.N)),[],'all','omitnan');      % 功率突增扰动下的最大频差平均积分电量总和边界迭代增量
icm_ImAve_SP2 = min(abs(G.ImdAve-G.ImdAve')+diag(NaN*ones(1,G.N)),[],'all','omitnan');      % 功率突减扰动下的最大频差平均积分电量总和边界迭代增量
icm_I60Ave_SP1 = min(abs(G.I60uAve-G.I60uAve')+diag(NaN*ones(1,G.N)),[],'all','omitnan');   % 功率突增扰动下的60s频差平均积分电量总和边界迭代增量
icm_I60Ave_SP2 = min(abs(G.I60dAve-G.I60dAve')+diag(NaN*ones(1,G.N)),[],'all','omitnan');   % 功率突减扰动下的60s频差平均积分电量总和边界迭代增量
flag_vfini_SP1 = zeros(T,1);    % 功率突增扰动下初始频率变化率越限标志，越限为1，未越限为0
flag_vfini_SP2 = zeros(T,1);    % 功率突减扰动下初始频率变化率越限标志，越限为1，未越限为0
flag_dfend_SP1 = zeros(T,1);    % 功率突增扰动下60s频差越限标志，越限为1，未越限为0
flag_dfend_SP2 = zeros(T,1);    % 功率突减扰动下60s频差越限标志，越限为1，未越限为0
flag_dfpeak_SP1 = zeros(T,1);   % 功率突增扰动下最大频差越限标志，越限为1，未越限为0
flag_dfpeak_SP2 = zeros(T,1);   % 功率突减扰动下最大频差越限标志，越限为1，未越限为0
%% 可行性优化割边界初值
% ————————————————————————手动迭代————————————————————————
load main_result.mat
Its = size(Iteration,2);
Hvf_SP1 = Iteration(Its).boundary.Hvf_SP1;
Hvf_SP2 = Iteration(Its).boundary.Hvf_SP2;
Hdf_SP1 = Iteration(Its).boundary.Hdf_SP1;
Hdf_SP2 = Iteration(Its).boundary.Hdf_SP2;
Im_SP1 = Iteration(Its).boundary.Im_SP1;
Im_SP2 = Iteration(Its).boundary.Im_SP2;
I60_SP1 = Iteration(Its).boundary.I60_SP1;
I60_SP2 = Iteration(Its).boundary.I60_SP2;
ImAve_SP1 = Iteration(Its).boundary.ImAve_SP1;
ImAve_SP2 = Iteration(Its).boundary.ImAve_SP2;
I60Ave_SP1 = Iteration(Its).boundary.I60Ave_SP1;
I60Ave_SP2 = Iteration(Its).boundary.I60Ave_SP2;
% ————————————————————————————————————————————————————————

% Hvf_SP1 = zeros(T,1);       % 功率突增扰动下由初始频率变化率决定的系统的惯量总和边界(状态边界)
% Hvf_SP2 = zeros(T,1);       % 功率突减扰动下由初始频率变化率决定的系统的惯量总和边界(状态边界)
% Hdf_SP1 = zeros(T,1);       % 功率突增扰动下由最大频差决定的系统的惯量总和边界(状态边界)
% Hdf_SP2 = zeros(T,1);       % 功率突减扰动下由最大频差决定的系统的惯量总和边界(状态边界)
% Im_SP1 = zeros(T,1);        % 功率突增扰动下给定机组组合形式的最大频差积分电量总和边界(出力边界)
% Im_SP2 = zeros(T,1);        % 功率突减扰动下给定机组组合形式的最大频差积分电量总和边界(出力边界)
% I60_SP1 = zeros(T,1);       % 功率突增扰动下给定机组组合形式的60s频差积分电量总和边界(出力边界)
% I60_SP2 = zeros(T,1);       % 功率突减扰动下给定机组组合形式的60s频差积分电量总和边界(出力边界)
% ImAve_SP1 = zeros(T,1);     % 功率突增扰动下系统的最大频差平均积分电量总和边界(状态边界)
% ImAve_SP2 = zeros(T,1);     % 功率突减扰动下系统的最大频差平均积分电量总和边界(状态边界)
% I60Ave_SP1 = zeros(T,1);    % 功率突增扰动下系统的60s频差平均积分电量总和边界(状态边界)
% I60Ave_SP2 = zeros(T,1);    % 功率突减扰动下系统的60s频差平均积分电量总和边界(状态边界)
% border.SP1.df60.u_G = cell(T,1);    % border记录了迭代过程中各时段出现过的、有能力避免频差指标(最大频差、
% border.SP1.df60.Imin = cell(T,1);   % 60s频差)越限的DPR机组组合形式及对应频差的积分电量总和的最小边界。
% border.SP1.dfm.u_G = cell(T,1);
% border.SP1.dfm.Imin = cell(T,1);
% border.SP2.df60.u_G = cell(T,1);
% border.SP2.df60.Imin = cell(T,1);
% border.SP2.dfm.u_G = cell(T,1);
% border.SP2.dfm.Imin = cell(T,1);
%% 迭代
while 1
    %% MP-SP迭代
    Its = Its+1;
    MP_Dispatch(Hvf_SP1,Hvf_SP2,Hdf_SP1,Hdf_SP2,ImAve_SP1,ImAve_SP2,I60Ave_SP1,I60Ave_SP2,Im_SP1,Im_SP2,I60_SP1,I60_SP2)
    load MP_result.mat
    for t = 1:T
        %% SP1(功率突增扰动)校核
        [vfini_SP1(t),dfend_SP1(t),dfpeak_SP1(t)] = frequency(u_G_MP(t,:),P_G_MP(t,:),P_L(t),step_W*sum(P_Wp(t,:))+step_L*P_L(t));
        % 判断初始频率变化率是否越限，并更新边界
        if abs(vfini_SP1(t))>vfmax
            Hvf_SP1(t) = sum(u_G_MP(t,:).*G.Pn.*G.H)+icm_H;
            flag_vfini_SP1(t) = 1;
        end
        % 判断60s频差是否越限，并更新边界
        record = 0;     % 记录标志位，有记录则置1
        if abs(dfend_SP1(t))>df60s
            for j = 1:length(border.SP1.df60.Imin{t})      % 查找border.SP1.df60记录，若有则直接调用记录结果
                if sum(round(u_G_MP(t,:))==round(border.SP1.df60.u_G{t}(j,:)))==G.N
                    [~,idx] = max(border.SP1.df60.Imin{t});
                    border.SP1.df60.Imin{t}(idx) = (1+step_I)*border.SP1.df60.Imin{t}(idx);
                    I60_SP1(t) = max(border.SP1.df60.Imin{t});
                    record = 1;
                    break
                end
            end
            if record == 0    % 未在border.SP1.df60中找到相同的记录，需要调用df_boundary函数计算边界
                [flag_ol,Imin] = df_boundary(u_G_MP(t,:),P_G_MP(t,:),t,1,1);
                if flag_ol == 1
                    I60Ave_SP1(t) = sum(u_G_MP(t,:).*G.I60uAve)+icm_I60Ave_SP1;
                else
                    border.SP1.df60.u_G{t} = [border.SP1.df60.u_G{t};u_G_MP(t,:)];
                    border.SP1.df60.Imin{t} = [border.SP1.df60.Imin{t},Imin];
                    I60_SP1(t) = max(border.SP1.df60.Imin{t});
                end
            end
            flag_dfend_SP1(t) = 1;
        end
        % 判断最大频差是否越限，并更新边界
        record = 0;
        if abs(dfpeak_SP1(t))>dfmax && abs(dfpeak_SP1(t))>abs(dfend_SP1(t))
            for j = 1:length(border.SP1.dfm.Imin{t})       % 查找border.SP1.dfm记录，若有则直接调用记录结果
                if sum(round(u_G_MP(t,:))==round(border.SP1.dfm.u_G{t}(j,:)))==G.N
                    [~,idx] = max(border.SP1.dfm.Imin{t});
                    border.SP1.dfm.Imin{t}(idx) = (1+step_I)*border.SP1.dfm.Imin{t}(idx);
                    Im_SP1(t) = max(border.SP1.dfm.Imin{t});
                    record = 1;
                    break
                end
            end
            if record == 0    % 未在border.SP1.dfm中找到相同的记录，需要调用df_boundary函数计算边界
                [flag_ol,Imin] = df_boundary(u_G_MP(t,:),P_G_MP(t,:),t,1,2);
                if flag_ol == 1
                    ImAve_SP1(t) = sum(u_G_MP(t,:).*G.ImuAve)+icm_ImAve_SP1;
                    Hdf_SP1(t) = sum(u_G_MP(t,:).*G.Pn.*G.H)+icm_H;
                else
                    border.SP1.dfm.u_G{t} = [border.SP1.dfm.u_G{t};u_G_MP(t,:)];
                    border.SP1.dfm.Imin{t} = [border.SP1.dfm.Imin{t},Imin];
                    Im_SP1(t) = max(border.SP1.dfm.Imin{t});
                end
            end
            flag_dfpeak_SP1(t) = 1;
        end
        %% SP2(功率突减扰动)校核
        [vfini_SP2(t),dfend_SP2(t),dfpeak_SP2(t)] = frequency(u_G_MP(t,:),P_G_MP(t,:),P_L(t),-step_W*sum(P_Wp(t,:))-step_L*P_L(t));
        % 判断初始频率变化率是否越限，并更新边界
        if abs(vfini_SP2(t))>vfmax
            Hvf_SP2(t) = sum(u_G_MP(t,:).*G.Pn.*G.H)+icm_H;
            flag_vfini_SP2(t) = 1;
        end
        % 判断60s频差是否越限，并更新边界
        record = 0;
        if abs(dfend_SP2(t))>df60s
            for j = 1:length(border.SP2.df60.Imin{t})      % 查找border.SP2.df60记录，若有则直接调用记录结果
                if sum(round(u_G_MP(t,:))==round(border.SP2.df60.u_G{t}(j,:)))==G.N
                    [~,idx] = max(border.SP2.df60.Imin{t});
                    border.SP2.df60.Imin{t}(idx) = (1+step_I)*border.SP2.df60.Imin{t}(idx);
                    I60_SP2(t) = max(border.SP2.df60.Imin{t});
                    record = 1;
                    break
                end
            end
            if record == 0    % 未在border.SP2.df60中找到相同的记录，需要调用df_boundary函数计算边界
                [flag_ol,Imin] = df_boundary(u_G_MP(t,:),P_G_MP(t,:),t,2,1);
                if flag_ol == 1
                    I60Ave_SP2(t) = sum(u_G_MP(t,:).*G.I60dAve)+icm_I60Ave_SP2;
                else
                    border.SP2.df60.u_G{t} = [border.SP2.df60.u_G{t};u_G_MP(t,:)];
                    border.SP2.df60.Imin{t} = [border.SP2.df60.Imin{t},Imin];
                    I60_SP2(t) = max(border.SP2.df60.Imin{t});
                end
            end
            flag_dfend_SP2(t) = 1;
        end
        % 判断最大频差是否越限，并更新边界
        record = 0;
        if abs(dfpeak_SP2(t))>dfmax && abs(dfpeak_SP2(t))>abs(dfend_SP2(t))
            for j = 1:length(border.SP2.dfm.Imin{t})       % 查找border.SP2.dfm记录，若有则直接调用记录结果
                if sum(round(u_G_MP(t,:))==round(border.SP2.dfm.u_G{t}(j,:)))==G.N
                    [~,idx] = max(border.SP2.dfm.Imin{t});
                    border.SP2.dfm.Imin{t}(idx) = (1+step_I)*border.SP2.dfm.Imin{t}(idx);
                    Im_SP2(t) = max(border.SP2.dfm.Imin{t});
                    record = 1;
                    break
                end
            end
            if record == 0    % 未在border.SP2.dfm中找到相同的记录，需要调用df_boundary函数计算边界
                [flag_ol,Imin] = df_boundary(u_G_MP(t,:),P_G_MP(t,:),t,2,2);
                if flag_ol == 1
                    ImAve_SP2(t) = sum(u_G_MP(t,:).*G.ImdAve)+icm_ImAve_SP2;
                    Hdf_SP2(t) = sum(u_G_MP(t,:).*G.Pn.*G.H)+icm_H;
                else
                    border.SP2.dfm.u_G{t} = [border.SP2.dfm.u_G{t};u_G_MP(t,:)];
                    border.SP2.dfm.Imin{t} = [border.SP2.dfm.Imin{t},Imin];
                    Im_SP2(t) = max(border.SP2.dfm.Imin{t});
                end
            end
            flag_dfpeak_SP2(t) = 1;
        end
    end
    %% 迭代结果记录
    disp(['已完成',num2str(Its),'次迭代，当前用时',num2str(toc),'秒。'])
    disp(['功率突增扰动下，初始频率变化率越限时段',num2str(sum(flag_vfini_SP1)),'个，最大频差越限时段',...
        num2str(sum(flag_dfpeak_SP1)),'个，60s频差越限时段',num2str(sum(flag_dfend_SP1)),'个。'])
    disp(['功率突减扰动下，初始频率变化率越限时段',num2str(sum(flag_vfini_SP2)),'个，最大频差越限时段',...
        num2str(sum(flag_dfpeak_SP2)),'个，60s频差越限时段',num2str(sum(flag_dfend_SP2)),'个。'])
    Iteration(Its).MP_result = load('MP_result.mat');   % 当前迭代轮次MP求解结果
    Iteration(Its).time = toc;                          % 截至当前迭代轮次的总用时
    Iteration(Its).findex.SP1.vfini = vfini_SP1';       % 当前迭代轮次的频率指标
    Iteration(Its).findex.SP1.dfend = dfend_SP1';
    Iteration(Its).findex.SP1.dfpeak = dfpeak_SP1';
    Iteration(Its).findex.SP2.vfini = vfini_SP2';
    Iteration(Its).findex.SP2.dfend = dfend_SP2';
    Iteration(Its).findex.SP2.dfpeak = dfpeak_SP2';
    Iteration(Its).boundary.Hvf_SP1 = Hvf_SP1;          % 当前迭代轮次的优化割约束边界
    Iteration(Its).boundary.Hvf_SP2 = Hvf_SP2;
    Iteration(Its).boundary.Hdf_SP1 = Hdf_SP1;
    Iteration(Its).boundary.Hdf_SP2 = Hdf_SP2;
    Iteration(Its).boundary.Im_SP1 = Im_SP1;
    Iteration(Its).boundary.Im_SP2 = Im_SP2;
    Iteration(Its).boundary.I60_SP1 = I60_SP1;
    Iteration(Its).boundary.I60_SP2 = I60_SP2;
    Iteration(Its).boundary.ImAve_SP1 = ImAve_SP1;
    Iteration(Its).boundary.ImAve_SP2 = ImAve_SP2;
    Iteration(Its).boundary.I60Ave_SP1 = I60Ave_SP1;
    Iteration(Its).boundary.I60Ave_SP2 = I60Ave_SP2;
    %% 迭代结束条件
    if solveinfo.problem ~= 0
        disp('MP不收敛。')
        break
    end
    if sum(flag_vfini_SP1+flag_vfini_SP2+flag_dfend_SP1+flag_dfend_SP2+flag_dfpeak_SP1+flag_dfpeak_SP2)==0
        break
    else
        flag_vfini_SP1 = zeros(T,1);
        flag_vfini_SP2 = zeros(T,1);
        flag_dfend_SP1 = zeros(T,1);
        flag_dfend_SP2 = zeros(T,1);
        flag_dfpeak_SP1 = zeros(T,1);
        flag_dfpeak_SP2 = zeros(T,1);
    end
end
%% 结果保存
save main_result Iteration border