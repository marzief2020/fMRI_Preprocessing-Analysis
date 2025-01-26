% run_names should be changed according to the session: ses-1 or ses-2

% this codes get the rp-*.txt files and draw the FD and MP based on them.
% check the order of rp files for each ses and make  sure the ref_run is chosen correctly.

function PreProc_MP_plot(task_dir, iEcho, filter_rp, ref_run)

% spm fmri
% switch ses
%     case 'ses-1'
%         run_names = {'rest', 'LPI-1', 'LPI-2', 'LPI-3', 'LPI-4',  'SL', 'ML' }; %ref_run = 3;
%     case 'ses-2'
%         run_names = {'rest', 'GL-1', 'GL-2', 'GL-3', 'GC-1', 'GC-2'}; %ref_run = ;
% end

rp_mat=[]; range_data = []; range_temp = []; nSCans_temp=[]; nScans_runs=[];
rp_files = spm_select('FPListRec', task_dir, filter_rp); % check the ref runs and choose the ref_run accordingly
nRuns = size(rp_files, 1);

func_fn_cell = cell(1, nRuns);
iTask=1;
for iRun=1:nRuns
    [d_rp_fn, e_rp_fn, f_rp_fn] =  fileparts(rp_files(iRun, :));
    func_fn = spm_select('FPlist', d_rp_fn, ['^sub.*nii$']);
    func_fn_cell(1, iRun) = {mat2cell(func_fn, ones(size(func_fn,1),1))};

    nScans_runs(iRun) = size(func_fn_cell{iRun} , 1) + sum(nSCans_temp);
    nSCans_temp(iRun) =  size(func_fn_cell{iRun}, 1);

    rp_temp = spm_load(rp_files(iRun,:));
    rp_mat = [rp_mat; rp_temp];

    task_name = e_rp_fn(22:end);
    if regexp(task_name, '^rest')
        run_names{iRun}= 'rest';

    elseif  regexp(task_name, '^LinearPI')
        run_names{iRun}= ['LPI-', num2str(iTask)];
        task_id= ['task-LinearPI'];
        iTask=iTask+1;

    elseif  regexp(task_name, '^SceneLocalizer')
        run_names{iRun}= 'SL';
        

    elseif  regexp(task_name, '^MotionLocalizer')
        run_names{iRun}= 'ML';
        

    elseif  regexp(task_name, '^GridLocalizer')
        run_names{iRun}= ['GL-', num2str(iTask)];
        task_id= ['task-GridLocalizer'];
        iTask=iTask+1;

    elseif  regexp(task_name, '^DelcodeGridCell')
        run_names{iRun}= ['GD-', num2str(iTask)];
        task_id= ['task-DelcodeGridCell'];
        iTask=iTask+1;
        
            
    end
    % d_parts = split(d_rp_fn, '\');
end


figure();
sgtitle (['MP\_runwise\_', e_rp_fn(4:9)])

dist_mm=(rp_mat(:,1:3));
dist_mm_range=abs(ceil(max(dist_mm, [], 'all')));
subplot(2,1,1);plot(dist_mm);
set(gca,'xlim',[0 size(rp_mat,1)+1]); yticks([-dist_mm_range :1: dist_mm_range])
ylabel('mm');

set(gca, 'xtick', nScans_runs, 'xticklabel',run_names);

angle_deg=(rp_mat(:,4:6).*180/pi);
angle_deg_range=abs(ceil(max(angle_deg,[], 'all')));
subplot(2,1,2);plot(angle_deg);
set(gca,'xlim',[0 size(rp_mat,1)+1]); yticks([-angle_deg_range :1: angle_deg_range])
ylabel('degrees')
set(gca, 'xtick', nScans_runs, 'xticklabel', run_names);

MP_rp_fn = ['MP_run_', e_rp_fn(4:9), '.png'];
% saveas (gcf, fullfile (func_dir, MP_rp_fn)) ;
close

%% generate rp files and plot based on the movements compared to the first scan

V= spm_vol(func_fn_cell);
Q=[]; Q_run=[]; FD=[]; FD_run=[]; MP=[]; MP_run=[];

for iRun=1:size(V, 2)

    [d_fn, e_fn, f_fn] = fileparts(V{iRun}{1}.fname);

    for iScan=1:size(V{iRun}, 1)
        P_temp(iScan, :) = spm_imatrix(V{iRun}{iScan}.mat/V{ref_run}{1}.mat);
        P_temp_run(iScan, :) = spm_imatrix(V{iRun}{iScan}.mat/V{iRun}{1}.mat);
    end

    Q_temp = P_temp(:,1:6);
    Q = [Q; Q_temp];

    Q_temp_run = P_temp_run(:,1:6);
    Q_run = [Q_run; Q_temp_run];

    save(fullfile(d_fn, ['rp_total_', e_fn(1:end-11), '.txt']),'Q_temp','-ascii');

    % Compute the absolute derivative of the six realignment parameters
    abs_diff_params = [zeros(1,6); abs(diff(Q_temp))];
    abs_diff_params_run = [zeros(1,6); abs(diff(Q_temp_run))];

    % Compute the framewise displacement (FD) as the sum of the absolute derivative parameters
    FD_temp = sum(abs_diff_params, 2);
    save(fullfile(d_fn, ['FD_', e_fn(1:end-11), '.txt']),'FD_temp','-ascii');

    FD_temp_run = sum(abs_diff_params_run, 2);
    save(fullfile(d_fn, ['FD_run_', e_fn(1:end-11), '.txt']),'FD_temp_run','-ascii');

    % Calculate the MP_abs
    trans = Q_temp(:, 1:3); % Extract the translational parameters
    rot = Q_temp(:, 4:6); % Extract the rotational parameters
    trans_run = Q_temp_run(:, 1:3); % Extract the translational parameters
    rot_run = Q_temp_run(:, 4:6); % Extract the rotational parameters
    % diffTrans = [0, 0, 0; abs(trans)]; % Calculate the difference in translational parameters
    % diffRot = [0, 0, 0; diff(rot)]; % Calculate the difference in rotational parameters
    MP_abs = sum(abs(trans), 2) + sum(abs(rot), 2); % Calculate the MP
    MP_abs_run = sum(abs(trans_run-trans_run(1,:)), 2) + sum(abs(rot_run-rot_run(1,:)), 2);
    save(fullfile(d_fn, ['MP_ref_', e_fn(1:end-11), '.txt']),'MP_abs','-ascii');
    save(fullfile(d_fn, ['MP_run_', e_fn(1:end-11), '.txt']),'MP_abs_run','-ascii');

    FD = [FD; FD_temp];
    FD_run = [FD_run; FD_temp_run];
    MP = [MP; MP_abs];
    MP_run = [MP_run; MP_abs_run];
    clear P_temp Q_temp FD_temp MP_abs MP_abs_run P_temp_run Q_temp_run

end

% save(fullfile(func_dir, ['movement_params.txt']),'Q','-ascii');

fontsizeL = 16; fontsizeM = 18;

fg = spm_figure('FindWin','Graphics');
if isempty(fg), return; end

spm_figure('Clear','Graphics');
ax = axes('Position',[0.1 0.90 0.8 0.1],'Parent',fg,'Visible','off');
set(get(ax,'Title'),'String','Image realignment',...
    'FontSize',16,'FontWeight','Bold','Visible','on');
x     =  0.1;
y     =  0.9;

nV= sum(cellfun(@numel, V));

ax = axes('Position',[0.1 0.35 0.8 0.15],'Parent',fg,'XGrid','on','YGrid','on',...
    'NextPlot','replacechildren','ColorOrder',[0 0 1;0 0.5 0;1 0 0]);
plot(Q(:,1:3),'Parent',ax)
s  = {'x translation','y translation','z translation'};
%text([2 2 2], Params(2, 1:3), s, 'Fontsize',10,'Parent',ax)
legend(ax, s, 'Location','Best')
set(get(ax,'Title'),'String','translation','FontSize',fontsizeL,'FontWeight','Bold');
set(get(ax,'Xlabel'),'String','image');
set(get(ax,'Ylabel'),'String','mm');
set(ax, 'XTick', nScans_runs)
set(ax, 'XTickLabel', run_names)
axis tight;

ax = axes('Position',[0.1 0.1 0.8 0.15],'Parent',fg,'XGrid','on','YGrid','on',...
    'NextPlot','replacechildren','ColorOrder',[0 0 1;0 0.5 0;1 0 0]);
plot(Q(:,4:6)*180/pi,'Parent',ax)
s  = {'pitch','roll','yaw'};
%text([2 2 2], Params(2, 4:6)*180/pi, s, 'Fontsize',10,'Parent',ax)
legend(ax, s, 'Location','Best')
set(get(ax,'Title'),'String','rotation','FontSize',fontsizeL,'FontWeight','Bold');
set(get(ax,'Xlabel'),'String','image');
set(get(ax,'Ylabel'),'String','degrees');
set(gcf, 'Resize', 'on')
axis tight;
set(ax, 'XTick', nScans_runs)
set(ax, 'XTickLabel', run_names)

ax = axes('Position',[0.1 0.60 0.8 0.1],'Parent',fg,'XGrid','on','YGrid','on',...
    'NextPlot','replacechildren','ColorOrder',[0 0 1;0 0.5 0;1 0 0]);
plot(ax, MP_run', 'LineWidth', 2); ylabel(ax, 'mm','fontsize',fontsizeL)
set(get(ax,'Title'),'String','MP_abs_run','FontSize',16,'FontWeight','Bold'); axis tight;

ax = axes('Position',[0.1 0.80 0.8 0.1],'Parent',fg,'XGrid','on','YGrid','on',...
    'NextPlot','replacechildren','ColorOrder',[0 0 1;0 0.5 0;1 0 0]);
plot(ax, MP', 'LineWidth', 2); ylabel(ax, 'mm','fontsize',fontsizeL)
set(get(ax,'Title'),'String','MP_abs_total','FontSize',16,'FontWeight','Bold'); axis tight;

sgtitle([task_id, ' - ', e_fn(1:6)], 'fontsize',fontsizeM)

MP_fn = ['MP_' task_id, '_', e_fn(1:6), '_','echo-', num2str(iEcho), '.png'];
saveas (gcf, fullfile (task_dir, MP_fn)) ;

end