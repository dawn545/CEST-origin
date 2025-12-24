% run_z_to_B1interp_simple.m
% 简化版脚本：从 zspect_dataavecorrect.mat 构建 Z_stack，然后做 B1 插值校正

clear; clc; close all;                                     % 清空工作区、命令行、关闭图窗
proj_root = 'e:/20240613_Rat_dead_muscle';                 % 项目根目录
series = 'PCr';                                            % 选择代谢物系列：PCr 或 Cr
fit_type = 'smoothingspline';                              % 拟合方式（平滑样条）
num_test_offsets = 1;                                      % 若>0 仅处理前几个频移点（调试用）
B1_input_index = [];                                       % 如果为空则默认使用全部B1输入

cd(proj_root);                                             % 切换到项目根目录
addpath(proj_root);                                        % 添加路径以加载函数与 MAT 文件

% 依据 series 设置 B1_input 及对应文件夹名称
if strcmp(series, 'PCr')
    B1_input = [0.35, 0.65, 0.95, 1.25, 1.55];              % PCr B1 激励列表
    folder_names = {...                                     % 对应扫描序列文件夹
        '053_CESTEPI_dead_PCr_0.35_1500_2000_20240613', ...
        '054_CESTEPI_dead_PCr_0.65_1500_2000_20240613', ...
        '055_CESTEPI_dead_PCr_0.95_1500_2000_20240613', ...
        '056_CESTEPI_dead_PCr_1.25_1500_2000_20240613', ...
        '057_CESTEPI_dead_PCr_1.55_1500_2000_20240613'};
elseif strcmp(series, 'Cr')
    B1_input = [0.35, 0.65, 0.95, 1.25, 1.55, 2.05, 2.55];   % Cr B1 激励列表
    folder_names = {...                                     % 对应扫描序列文件夹
        '059_CESTEPI_dead_Cr_0.35_1500_2000_20240613', ...
        '060_CESTEPI_dead_Cr_0.65_1500_2000_20240613', ...
        '061_CESTEPI_dead_Cr_0.95_1500_2000_20240613', ...
        '062_CESTEPI_dead_Cr_1.25_1500_2000_20240613', ...
        '063_CESTEPI_dead_Cr_1.55_1500_2000_20240613', ...
        '064_CESTEPI_dead_Cr_2.05_1500_2000_20240613', ...
        '065_CESTEPI_dead_Cr_2.55_1500_2000_20240613'};
else
    error('Unknown series');                                % series 输入错误时报错
end

% 读取 B1-map
b1map_file = fullfile(proj_root,'B1-map.mat');
if ~exist(b1map_file,'file')
    error('B1-map.mat not found');
end
S = load(b1map_file);                                      % 加载文件
fn = fieldnames(S);                                        % 读取变量名（不关心名字是什么）
rel_B1map = S.(fn{1});                                      % 提取 B1 map

% 初始化图像大小参数（若未找到会被后续覆盖）
ny = 96; nx = 70; nz = 1; num_offset = 0;
num_B1 = numel(B1_input);                                   % B1 序列数量
found_zs = false;

fprintf('Searching for zspect_dataavecorrect.mat in poly-fit-water subdirectories...\n');
% 遍历每个 B1 文件夹查找 zspect_dataavecorrect.mat 并读取大小
for b1i = 1:num_B1
    folder = fullfile(proj_root, folder_names{b1i}, 'poly-fit-water');
    zsfile_local = fullfile(folder,'zspect_dataavecorrect.mat');
    if exist(zsfile_local,'file')
        fprintf('Found: %s\n', zsfile_local);
        tmp = load(zsfile_local);
        vars = fieldnames(tmp);
        for i=1:numel(vars)
            v = tmp.(vars{i});
            if isnumeric(v) && ndims(v)==3                         % 找到 3D 数据（y,x,offset）
                zs = v; found_zs = true; break;
            end
        end
        if found_zs
            [ny,nx,num_offset] = size(zs);                         % 读取大小
            nz = 1;
            fprintf('  -> Shape: %d x %d x %d (y x x x offset)\n', ny,nx,num_offset);
            break;                                                % 找到一个就停止
        end
    end
end

if ~found_zs
    error('zspect_dataavecorrect.mat not found');
end

% 创建 Z_stack (y,x,z,offset,B1)
Z_stack = zeros(ny,nx,nz,num_offset,num_B1,'single');

% 把每个 B1 文件夹的 zspect 数据填入 Z_stack
for b1i = 1:num_B1
    folder = fullfile(proj_root, folder_names{b1i}, 'poly-fit-water');
    zsfile_local = fullfile(folder,'zspect_dataavecorrect.mat');
    if exist(zsfile_local,'file')
        tmp = load(zsfile_local);
        vars = fieldnames(tmp);
        found = false;
        for vi=1:numel(vars)
            v = tmp.(vars{vi});
            if isnumeric(v) && ndims(v)==3 && all(size(v,1:2)==[ny,nx])
                Z_stack(:,:,:,:,b1i) = single(v);                  % 写入当前 B1 层
                found = true; break;
            end
        end
        if ~found
            error('No zspect_dataavecorrect in folder');
        end
    end
end

% 若只做调试可限制 offset 数量
if num_test_offsets>0
    num_offset_proc = min(num_test_offsets, size(Z_stack,4));
    Z_stack_proc = Z_stack(:,:,:,1:num_offset_proc,:);
else
    Z_stack_proc = Z_stack;
end

% 默认输出 B1 与输入一致
B1_output = B1_input;
if isempty(B1_input_index), B1_input_index = 1:num_B1; end
SEGMENT = ones(ny,nx,nz);                                  % mask 全1（不做分区）

% B1 插值矫正
Z_stack_corr = Z_B1_correction(Z_stack_proc, rel_B1map, B1_input, B1_output, SEGMENT, fit_type, B1_input_index);

% 保存结果
out_file = fullfile(proj_root, sprintf('Z_stack_corr_%s_simple.mat', series));
save(out_file,'Z_stack_corr','Z_stack','B1_input','B1_output','rel_B1map','-v7.3');
fprintf('Saved: %s\n', out_file);

% 构建重建后的 Z_recon（按 B1 输出组织）
Z_recon = Z_stack_corr;
if size(Z_recon,4)>1
    Z_recon_mean = squeeze(mean(Z_recon,4));               % 多频移 → 按 offset 求均值
else
    Z_recon_mean = squeeze(Z_recon);                       % 无 offset → 直接 squeeze
end

save(fullfile(proj_root,sprintf('Z_recon_mean_%s.mat',series)),'Z_recon_mean','-v7.3');
fprintf('Saved reconstructed Z to Z_recon_mean_%s.mat\n',series);

% 绘制并保存一个像素校正前后对比（默认使用图像中心）
% 可选：将插值类型设置为 'pchip' 或 'smoothingspline'
interp_plot_type = 'pchip';
cy = round(ny/2); cx = round(nx/2); cz = 1; coff = 1;
orig_vals = squeeze(double(Z_stack_proc(cy,cx,cz,coff,:)));    % 原数据
corr_vals = squeeze(double(Z_stack_corr(cy,cx,cz,coff,:)));    % 校正数据
rel_pixel = rel_B1map(cy,cx);                                  % 该像素的相对B1

% Use nominal B1 (as provided in B1_input) for plotting by default
b1_nominal = double(B1_input(:));                              % 名义 B1 值
b1_x = b1_nominal;                                              % 横坐标：名义 B1（不乘 rel_pixel）
absB1_output_pixel = double(B1_output(:)) ;                      % 校正后点横坐标（名义 B1）

% 1) 计算两种拟合/插值：smoothing spline（与原脚本一致）和 pchip（穿点插值）
%    并在原始点位置上求值以方便对比
try
    % smoothing spline（与原脚本兼容，用作对照）
    fo_ = fitoptions('method','SmoothingSpline','SmoothingParam',0.998);
    ft_ = fittype('SmoothingSpline');
    cf = fit(b1_x, orig_vals(:), ft_, fo_);
    fit_at_orig = cf(b1_x);
catch
    fit_at_orig = interp1(b1_x, orig_vals, b1_x, 'pchip'); % 兜底
end

% pchip（穿点插值，可保证拟合曲线通过原始采样点）
pchip_at_orig = interp1(b1_x, orig_vals, b1_x, 'pchip');

% 保存诊断表格（便于检查）
% include rel_pixel in diagnostic table so user can see scaling
rel_val = double(rel_pixel);
T = table((1:numel(b1_x))', B1_input(:), b1_x, repmat(rel_val,numel(b1_x),1), orig_vals, fit_at_orig, pchip_at_orig, corr_vals, ...
    'VariableNames', {'idx','B1_nominal','B1_plot_x','rel_pixel','orig','fit_at_orig','pchip_at_orig','corr'});
diag_matfile = fullfile(proj_root, sprintf('B1_interp_diag_pixel_%d_%d_offset%d_%s.mat', cy, cx, coff, series));
save(diag_matfile, 'T');
try
    writetable(T, fullfile(proj_root, sprintf('B1_interp_diag_pixel_%d_%d_offset%d_%s.csv', cy, cx, coff, series)));
catch
    % 如果 writetable 不可用则忽略
end

% 计算误差指标
fit_error = fit_at_orig - double(orig_vals);
pchip_error = pchip_at_orig - double(orig_vals);
rmse_fit = sqrt(mean(fit_error.^2));
rmse_pchip = sqrt(mean(pchip_error.^2));
fprintf('Diagnostic saved to %s (CSV and MAT).\n', diag_matfile);
fprintf('RMSE fit (smoothing spline) = %.6g, RMSE pchip = %.6g\n', rmse_fit, rmse_pchip);

% 2) 绘图：使用 pchip 作为默认可视化（保证通过点），同时绘制校正后点
% Interpolate over the nominal B1 range so x-axis stays within provided B1_input
% ---- x 轴定义 ----
b1_nominal = double(B1_input(:));              % 蓝线：名义 B1
b1_actual  = b1_nominal * rel_pixel./2;       % 黄线：实际 B1

% ---- 拟合原始曲线（名义 B1）----
B1_fine_nom = linspace(min(b1_nominal), max(b1_nominal), 200);
pchip_nom   = interp1(b1_nominal, orig_vals, B1_fine_nom, 'pchip');

% ---- 绘图 ----
fig = figure('Visible','off'); hold on; grid on;

% 蓝线：原始 Z-B1（名义）
plot(b1_nominal, orig_vals, 'ro', ...
    'MarkerSize',8, 'DisplayName','orig (nominal B1)');
plot(B1_fine_nom, pchip_nom, '-b', ...
    'LineWidth',1.8, 'DisplayName','orig fit (nominal B1)');

% 黄线：校正后 Z（实际 B1）
plot(b1_actual, corr_vals, 's-', ...
    'Color',[0.85 0.65 0.13], ...
    'LineWidth',2, ...
    'MarkerSize',7, ...
    'DisplayName','B1-correct (actual B1)');

xlabel('B1 (\muT)');
ylabel('Z (a.u.)');
title(sprintf('%s pixel(%d,%d) offset#%d  relB1=%.3g', ...
    series, cy, cx, coff, rel_pixel));
legend('Location','best');

out_fig = fullfile(proj_root, ...
    sprintf('B1_interp_example_%s_simple_pixel_%d_%d.png',series,cy,cx));
saveas(fig, out_fig);
close(fig);


fprintf('Example plot saved to %s\n', out_fig);
fprintf('Done.\n');                                            % 结束
