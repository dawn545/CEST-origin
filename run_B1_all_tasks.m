% run_B1_all_tasks.m
% 整合版：
%1) 构建 Z_stack
%2) B1 校正
%3) 三个验证实验
%4) B1 数量优化探索
%5) 导出 B1_nom / B1_real（供别组使用）

clear; clc; close all;

proj_root = 'E:/20240613_Rat_dead_muscle';
series = 'PCr'; % 'PCr' or 'Cr'
fit_type = 'smoothingspline'; % 'smoothingspline' / 'linear' / 'spline' / poly2~poly5
run_experiments = true;
run_optimization = true;

cd(proj_root);
addpath(proj_root);

%% ===================== B1 设置 =====================
if strcmp(series,'PCr')
 B1_input = [0.35,0.65,0.95,1.25,1.55];
 folder_names = {
 '053_CESTEPI_dead_PCr_0.35_1500_2000_20240613'
 '054_CESTEPI_dead_PCr_0.65_1500_2000_20240613'
 '055_CESTEPI_dead_PCr_0.95_1500_2000_20240613'
 '056_CESTEPI_dead_PCr_1.25_1500_2000_20240613'
 '057_CESTEPI_dead_PCr_1.55_1500_2000_20240613'};
elseif strcmp(series,'Cr')
 B1_input = [0.35,0.65,0.95,1.25,1.55,2.05,2.55];
 folder_names = {
 '059_CESTEPI_dead_Cr_0.35_1500_2000_20240613'
 '060_CESTEPI_dead_Cr_0.65_1500_2000_20240613'
 '061_CESTEPI_dead_Cr_0.95_1500_2000_20240613'
 '062_CESTEPI_dead_Cr_1.25_1500_2000_20240613'
 '063_CESTEPI_dead_Cr_1.55_1500_2000_20240613'
 '064_CESTEPI_dead_Cr_2.05_1500_2000_20240613'
 '065_CESTEPI_dead_Cr_2.55_1500_2000_20240613'};
else
 error('Unknown series');
end
num_B1 = numel(B1_input);

%% =====================选择数据子目录 =====================
candidate_subdirs = {'poly-fit-40-water','poly-fit-water'};
data_subdir = '';
for i =1:numel(candidate_subdirs)
 f0 = fullfile(proj_root, folder_names{1}, candidate_subdirs{i}, 'zspect_dataavecorrect.mat');
 if exist(f0, 'file')
 data_subdir = candidate_subdirs{i};
 break;
 end
end
if isempty(data_subdir)
 error('No zspect_dataavecorrect.mat found in poly-fit-40-water or poly-fit-water');
end
fprintf('Using data subdir: %s\n', data_subdir);

%% =====================读取 B1 map =====================
S = load(fullfile(proj_root,'B1-map.mat'));
fn = fieldnames(S);
B1_map_raw = S.(fn{1});
rel_B1map = B1_map_raw /2; % txt要求：先除以2得到rel-map

%% =====================读取 offset =====================
offset_file = fullfile(proj_root,folder_names{1},data_subdir,'offset.mat');
if ~exist(offset_file,'file')
 error('offset.mat not found: %s', offset_file);
end
tmp_off = load(offset_file);
fn_offset = fieldnames(tmp_off);
offset = tmp_off.(fn_offset{1});

%% ===================== 探测并读取 Z_stack =====================
found = false;
for b1i=1:num_B1
 f = fullfile(proj_root,folder_names{b1i},data_subdir,'zspect_dataavecorrect.mat');
 if exist(f,'file')
 tmp = load(f);
 vars = fieldnames(tmp);
 for k=1:numel(vars)
 v = tmp.(vars{k});
 if isnumeric(v) && ndims(v)==3
 [ny,nx,num_offset] = size(v);
 nz =1;
 found = true;
 break;
 end
 end
 end
 if found, break; end
end
if ~found, error('No zspect data found'); end

Z_stack = zeros(ny,nx,nz,num_offset,num_B1,'single');
for b1i=1:num_B1
 f = fullfile(proj_root,folder_names{b1i},data_subdir,'zspect_dataavecorrect.mat');
 if ~exist(f,'file')
 error('Missing zspect file: %s', f);
 end
 tmp = load(f);
 vars = fieldnames(tmp);
 loaded = false;
 for k=1:numel(vars)
 v = tmp.(vars{k});
 if isnumeric(v) && ndims(v)==3
 Z_stack(:,:,:,:,b1i) = single(v);
 loaded = true;
 break;
 end
 end
 if ~loaded
 error('No3D numeric var in %s', f);
 end
end

SEGMENT = ones(ny,nx,nz);
B1_output = B1_input;

%% ===================== 基准B1校正 =====================
% 注意：Z_B1_correction内部会再/2，因此这里传B1_map_raw，等效于使用(B1_map_raw/2)
Z_stack_corr_full = Z_B1_correction(...
 Z_stack, B1_map_raw, B1_input, B1_output, SEGMENT, fit_type,1:num_B1);

save(fullfile(proj_root,sprintf('allinone_baseline_%s.mat',series)), ...
 'Z_stack','Z_stack_corr_full','B1_input','B1_map_raw','rel_B1map','offset','-v7.3');

%% ===================== 像素诊断图（实际B1） =====================
cy = round(ny/2);
cx = round(nx/2);
cz =1;
coff = round(num_offset/2);

orig_vals = squeeze(double(Z_stack(cy,cx,cz,coff,:)));
corr_vals = squeeze(double(Z_stack_corr_full(cy,cx,cz,coff,:)));
rel_pixel = double(rel_B1map(cy,cx));
B1_real_pixel = double(B1_input(:))' * rel_pixel;

fig_diag = figure('Visible','off'); hold on; grid on;
plot(B1_real_pixel, orig_vals, 'ro', 'MarkerSize',8, 'DisplayName', 'orig');
plot(B1_real_pixel, corr_vals, 's-', 'LineWidth',2, 'MarkerSize',7, 'DisplayName', 'B1-corrected');
xlabel('B1_{actual} (\muT)');
ylabel('Z (a.u.)');
title(sprintf('%s pixel(%d,%d) offset#%d', series, cy, cx, coff));
legend('Location','best');
saveas(fig_diag, fullfile(proj_root, sprintf('allinone_diag_%s.png',series)));
close(fig_diag);

%% ===================== 三个验证实验 =====================
if run_experiments
 fprintf('\n========== 开始运行验证实验 ==========\n');

 % 实验1：删除某个偏频值的B1数据
 offset_to_remove = round(num_offset/2);
 Z_stack_exp1 = Z_stack;
 Z_stack_exp1(:,:,:,offset_to_remove,:) = [];
 offset_indices_exp1 = [1:offset_to_remove-1, offset_to_remove+1:num_offset];

 Z_stack_corr_exp1 = Z_B1_correction(...
 Z_stack_exp1, B1_map_raw, B1_input, B1_output, SEGMENT, fit_type,1:num_B1);

 diff_exp1 = Z_stack_corr_full(:,:,:,offset_indices_exp1,:) - Z_stack_corr_exp1;
 rmse_exp1 = sqrt(mean(diff_exp1(:).^2, 'omitnan'));

 save(fullfile(proj_root,sprintf('exp1_results_%s.mat',series)), ...
 'Z_stack_corr_exp1','diff_exp1','rmse_exp1','offset_to_remove','offset_indices_exp1','-v7.3');

 % 实验2：删除某个B1，插值重建
 B1_to_remove = round(num_B1/2);
 B1_indices_keep = setdiff(1:num_B1, B1_to_remove);
 B1_input_reduced = B1_input(B1_indices_keep);

 Z_stack_reduced = Z_stack(:,:,:,:,B1_indices_keep);
 B1_output_exp2 = B1_input(B1_to_remove);

 Z_stack_interp = Z_B1_correction(...
 Z_stack_reduced, B1_map_raw, B1_input_reduced, B1_output_exp2, SEGMENT, fit_type,1:numel(B1_indices_keep));

 Z_original = Z_stack(:,:,:,:,B1_to_remove);
 diff_exp2 = Z_original - Z_stack_interp;
 rmse_exp2 = sqrt(mean(diff_exp2(:).^2, 'omitnan'));

 Z_stack_reconstructed = Z_stack;
 Z_stack_reconstructed(:,:,:,:,B1_to_remove) = Z_stack_interp;

 Z_stack_corr_exp2 = Z_B1_correction(...
 Z_stack_reconstructed, B1_map_raw, B1_input, B1_output, SEGMENT, fit_type,1:num_B1);

 diff_corr_exp2 = Z_stack_corr_full - Z_stack_corr_exp2;
 rmse_corr_exp2 = sqrt(mean(diff_corr_exp2(:).^2, 'omitnan'));

 B1_nom_missing = B1_output_exp2;
 B1_real_missing = B1_output_exp2 .* rel_B1map;

 save(fullfile(proj_root,sprintf('exp2_results_%s.mat',series)), ...
 'Z_stack_interp','Z_stack_reconstructed','diff_exp2','rmse_exp2','rmse_corr_exp2', ...
 'B1_to_remove','B1_nom_missing','B1_real_missing','-v7.3');

 % 实验3：大误差点用原始值替代
 threshold = mean(abs(diff_exp2(:)), 'omitnan') +2*std(abs(diff_exp2(:)), 'omitnan');
 large_error_mask = abs(diff_exp2) > threshold;
 num_replaced = sum(large_error_mask(:));

 Z_stack_hybrid = Z_stack_interp;
 Z_stack_hybrid(large_error_mask) = Z_original(large_error_mask);

 Z_stack_reconstructed_hybrid = Z_stack;
 Z_stack_reconstructed_hybrid(:,:,:,:,B1_to_remove) = Z_stack_hybrid;

 Z_stack_corr_exp3 = Z_B1_correction(...
 Z_stack_reconstructed_hybrid, B1_map_raw, B1_input, B1_output, SEGMENT, fit_type,1:num_B1);

 diff_corr_exp3 = Z_stack_corr_full - Z_stack_corr_exp3;
 rmse_corr_exp3 = sqrt(mean(diff_corr_exp3(:).^2, 'omitnan'));

 B1_nom_hybrid = B1_input(B1_to_remove);
 B1_real_hybrid = B1_nom_hybrid .* rel_B1map;

 save(fullfile(proj_root,sprintf('exp3_results_%s.mat',series)), ...
 'Z_stack_hybrid','Z_stack_reconstructed_hybrid','diff_corr_exp3','rmse_corr_exp3', ...
 'threshold','num_replaced','B1_nom_hybrid','B1_real_hybrid','-v7.3');

 % 验证图（全部实际B1横轴）
 test_offset_full = round(num_offset/2);
 if test_offset_full == offset_to_remove
 if test_offset_full < num_offset
 test_offset_full = test_offset_full +1;
 else
 test_offset_full = test_offset_full -1;
 end
 end
 test_offset_exp1 = find(offset_indices_exp1 == test_offset_full,1);
 B1_plot = double(B1_input(:))' * rel_pixel;

 fig_val = figure('Position',[100,100,1200,800], 'Visible','off');

 subplot(2,2,1);
 plot(B1_plot, squeeze(Z_stack_corr_full(cy,cx,cz,test_offset_full,:)), 'o-', 'LineWidth',2, 'DisplayName', '完整数据'); hold on;
 plot(B1_plot, squeeze(Z_stack_corr_exp1(cy,cx,cz,test_offset_exp1,:)), 's-', 'LineWidth',2, 'DisplayName', '删除偏频后');
 xlabel('B1_{actual} (\muT)'); ylabel('Z值'); title(sprintf('实验1: 删除偏频#%d', offset_to_remove)); legend('Location','best'); grid on;

 subplot(2,2,2);
 plot(B1_plot, squeeze(Z_stack(cy,cx,cz,test_offset_full,:)), 'o-', 'LineWidth',2, 'DisplayName', '原始数据'); hold on;
 plot(B1_plot(B1_indices_keep), squeeze(Z_stack_reduced(cy,cx,cz,test_offset_full,:)), 'x', 'MarkerSize',10, 'LineWidth',2, 'DisplayName', '保留的B1');
 plot(B1_plot(B1_to_remove), squeeze(Z_stack_interp(cy,cx,cz,test_offset_full)), 'd', 'MarkerSize',10, 'LineWidth',2, 'DisplayName', '插值结果');
 xlabel('B1_{actual} (\muT)'); ylabel('Z值'); title(sprintf('实验2: 插值B1#%d', B1_to_remove)); legend('Location','best'); grid on;

 subplot(2,2,3);
 plot(B1_plot, squeeze(Z_stack(cy,cx,cz,test_offset_full,:)), 'o-', 'LineWidth',2, 'DisplayName', '原始数据'); hold on;
 plot(B1_plot(B1_to_remove), squeeze(Z_stack_hybrid(cy,cx,cz,test_offset_full)), 's', 'MarkerSize',10, 'LineWidth',2, 'DisplayName', '混合结果');
 xlabel('B1_{actual} (\muT)'); ylabel('Z值'); title(sprintf('实验3: 混合策略 (替换%d点)', num_replaced)); legend('Location','best'); grid on;

 subplot(2,2,4);
 bar([rmse_exp1, rmse_corr_exp2, rmse_corr_exp3]);
 set(gca, 'XTickLabel', {'实验1', '实验2', '实验3'});
 ylabel('RMSE'); title('三个实验的RMSE比较'); grid on;

 saveas(fig_val, fullfile(proj_root, sprintf('validation_experiments_%s.png', series)));
 close(fig_val);

 fprintf('实验完成: RMSE1=%.6f, RMSE2=%.6f, RMSE3=%.6f\n', rmse_exp1, rmse_corr_exp2, rmse_corr_exp3);
end

%% ===================== B1数量优化探索 =====================
if run_optimization
 fprintf('\n========== 开始B1数量优化探索 ==========\n');
 optimization_results = struct();

 % 策略1:4个B1插值1个
 if num_B1 >=5
 B1_remove_idx = round(num_B1/2);
 B1_keep_idx = setdiff(1:num_B1, B1_remove_idx);
 B1_keep_idx = B1_keep_idx(1:4);

 Z_reduced_s1 = Z_stack(:,:,:,:,B1_keep_idx);
 B1_input_s1 = B1_input(B1_keep_idx);
 B1_output_s1 = B1_input(B1_remove_idx);

 Z_interp_s1 = Z_B1_correction(Z_reduced_s1, B1_map_raw, B1_input_s1, B1_output_s1, SEGMENT, fit_type,1:4);
 diff_s1 = Z_stack(:,:,:,:,B1_remove_idx) - Z_interp_s1;
 optimization_results.strategy1.rmse = sqrt(mean(diff_s1(:).^2, 'omitnan'));
 optimization_results.strategy1.B1_keep = B1_keep_idx;
 optimization_results.strategy1.B1_remove = B1_remove_idx;
 end

 % 策略2:3个B1插值2个
 if num_B1 >=5
 B1_remove_idx = [round(num_B1/3), round(2*num_B1/3)];
 B1_keep_idx = setdiff(1:num_B1, B1_remove_idx);
 B1_keep_idx = B1_keep_idx(1:min(3,length(B1_keep_idx)));

 Z_reduced_s2 = Z_stack(:,:,:,:,B1_keep_idx);
 B1_input_s2 = B1_input(B1_keep_idx);
 B1_output_s2 = B1_input(B1_remove_idx);

 Z_interp_s2 = Z_B1_correction(Z_reduced_s2, B1_map_raw, B1_input_s2, B1_output_s2, SEGMENT, fit_type,1:length(B1_keep_idx));
 diff_s2 = Z_stack(:,:,:,:,B1_remove_idx) - Z_interp_s2;
 optimization_results.strategy2.rmse = sqrt(mean(diff_s2(:).^2, 'omitnan'));
 optimization_results.strategy2.B1_keep = B1_keep_idx;
 optimization_results.strategy2.B1_remove = B1_remove_idx;
 end

 % 策略3:3个B1插值1个
 if num_B1 >=4
 B1_remove_idx = round(num_B1/2);
 B1_keep_idx = setdiff(1:num_B1, B1_remove_idx);
 B1_keep_idx = B1_keep_idx(1:min(3,length(B1_keep_idx)));

 Z_reduced_s3 = Z_stack(:,:,:,:,B1_keep_idx);
 B1_input_s3 = B1_input(B1_keep_idx);
 B1_output_s3 = B1_input(B1_remove_idx);

 Z_interp_s3 = Z_B1_correction(Z_reduced_s3, B1_map_raw, B1_input_s3, B1_output_s3, SEGMENT, fit_type,1:length(B1_keep_idx));
 diff_s3 = Z_stack(:,:,:,:,B1_remove_idx) - Z_interp_s3;
 optimization_results.strategy3.rmse = sqrt(mean(diff_s3(:).^2, 'omitnan'));
 optimization_results.strategy3.B1_keep = B1_keep_idx;
 optimization_results.strategy3.B1_remove = B1_remove_idx;
 end

 save(fullfile(proj_root,sprintf('optimization_results_%s.mat',series)), 'optimization_results', '-v7.3');
end

fprintf('\n========== 所有任务完成 ==========' );
fprintf('\n输出文件：allinone_baseline / exp1-3 / optimization / validation图\n');
