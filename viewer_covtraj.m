function viewer_covtraj(params)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot realtime covariance matrix as trajectories of
% coordinates (autocov_A, autocov_B, crosscov_AB).
% AB is any unique combination of sources from an 
% input signals matrix of size samples-by-sources.
%
% stefano.orsolini@gmail.com
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

close all

addpath(genpath(fullfile(userpath,'FilterM')))

% total number of samples and sources
[Nsamp, Nsour] = size(params.signals);

% number of trajectories
Ntraj = ((Nsour^2)-Nsour)/2;

% generate unique pairwise enumeration
S_pairs = combnk(1:Nsour,2);

% buffer size in samples
b_span = ceil(params.t_show * params.Fs);

% raw signals buffer
b_raw = zeros(b_span,Nsour);
% preprocessed signals buffer
b_m = zeros(b_span,Nsour);
% space components buffer
b_ry = zeros(Nsour,Nsour,b_span);

% instance reverse time axis in seconds
t_b = linspace(-params.t_show,0,b_span);

% instance time plot figure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if params.show_timeseries
    f1_handle = figure('Color','w');
    
    for j=1:Nsour
        s1_handle(j) = subplot(Nsour,1,j);
        % instance buffer with initial value zero
        p1_handle(j) = plot(t_b,b_m(:,j),'LineWidth',1.5);
    end
end

% instance plot3 figure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f2_handle = figure('Color','w');
s2_handle = subplot(1,1,1);

for j=1:Ntraj
    p2_handle(j) = plot3(zeros(b_span,1), zeros(b_span,1), zeros(b_span,1),'LineWidth',1.5);
    hold on
    if params.show_details
        % details on last sample
        tx2_handle(j) = text(0,0,0,num2str(S_pairs(j,:)));
    else
        % spot on last sample
        sc2_handle(j) = scatter3(0,0,0,25,'filled','b');
    end
end

xlabel('auto-cov A')
ylabel('auto-cov B')
zlabel('cross-cov AB')

if params.show_plane
    % show boundaries of crosscovariance
    [Xp,Yp] = meshgrid(linspace(0,params.ac_lim,50));

    Zp = sqrt(Xp.*Yp);
    mesh(Xp,Yp,Zp,'FaceColor','k','FaceAlpha',0.05,'EdgeColor','k','EdgeAlpha',0.1)
    Zp = - sqrt(Xp.*Yp);
    mesh(Xp,Yp,Zp,'FaceColor','k','FaceAlpha',0.05,'EdgeColor','k','EdgeAlpha',0.1)
    
    if params.fixed_axes
        % reset axes limits
        reset_plot3_lims(s2_handle,params);
    end
end

% %Add slider: max power colorbar
% sl2_handle = uicontrol('Parent',f2_handle,'Style','slider','value',params.timescale, 'min',1, 'max',1000);
% % update slider position
% place_slider(f2_handle,sl2_handle);

view(80,10)
axis square
hold off
drawnow

% % complete colormaps list in:
% % fullfile(matlabroot,'toolbox/matlab/graph3d/Contents.m')
% cmap_list = {'parula','hsv','hot','gray','bone','copper','pink', ...
%              'jet','prism','cool','autumn','spring','winter','summer'};
% % extend list for all necessary trajectories
% cmap_list = repmat(cmap_list,ceil(Ntraj/length(cmap_list)),1);
%
% for j=1:Ntraj
%     custom_cmap = eval([cmap_list{j},'(',num2str(b_span),')']);
%     custom_cmap = [flipud(uint8(custom_cmap*255)) uint8(ones(b_span,1))].';
%     set(p2_handle(j).Edge,'ColorBinding','interpolated','ColorData',custom_cmap)
% end

if params.show_box
    % text box for experimental paradigm label
    b2_handle(1) = uicontrol(f2_handle,'style','text','HorizontalAlignment','left','FontWeight','bold','BackgroundColor','w');
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if params.do_global
    % buffer of the OLS model components for global component removal
    EV = [ones(b_span,1), zeros(b_span,1)];
end

if params.do_filter && ~params.do_global
    % transfer funtion coefficients (remove first from a)
    b = params.b.';
    a = params.a(2:end).';
    % input buffer of the filter (Moving-Average side)
    B_fi = zeros(Nsour,length(b));
    % output buffer of the filter (Auto-Regressive side)
    B_fo = zeros(Nsour,length(a));
end

%signal counter
i = 0;
%frame counter
n = 0;
n_rate = floor(params.Fs/params.FPS);
% fake sampling time
tic
while i ~= Nsamp
    % increment signal counter
    i = i+1;
    % increment video frame counter
    n = n+1;
    
    % update buffer: remove oldest sample, concatenate newest
    s_new = params.signals(i,:);
    b_raw = [b_raw(2:end,:); s_new];
    
    if ~params.do_global && ~params.do_filter
        % use original signal
        b_m = b_raw;
    end
    
    if params.do_global
        % the global component is the spatial average value
        EV = [EV(2:end,:); 1, mean(s_new)];
        % calculate Ordinary Least Squares coefficients
        betas = (EV'*EV)\(EV'*b_raw);
        % use OLS residual
        b_m = b_raw - (EV(:,1)*betas(1,:) + EV(:,2)*betas(2,:));
        
        if params.do_filter
            % cascade bandpass filter
            b_m = FiltFiltM(params.b,params.a,b_m);
        end
    end
    
    if params.do_filter && ~params.do_global
        % update filter input buffer with newest sample from signal buffer
        B_fi = [s_new.', B_fi(:,1:end-1)];
        
        % result of Moving-Average side (scalar product!)
        sMA = B_fi * b;
        % result of Auto-Regressive side (scalar product!)
        sAR = B_fo * a;
        
        % filter output value
        filt_out = sMA - sAR;
        
        % output buffer of the filter
        B_fo = [filt_out, B_fo(:,1:end-1)];
        
        % use filtered signal
        b_m = [b_m(2:end,:); filt_out.'];
    end
    
    % covariance matrix of the signals buffer
    C = cov(b_m);

    % update components buffer
    b_ry = cat(3, b_ry(:,:,2:end), C);
    
    % check out flag
    get_out_1 = exist('f1_handle','var') && ~ishandle(f1_handle);
    get_out_2 = exist('f2_handle','var') && ~ishandle(f2_handle);
    if get_out_1 || get_out_2
        break
    elseif n == n_rate
        % reset video frame counter
        n = 0;
        
        if params.show_box
            % update experimental paradigm label
            label = params.events(find(params.events(:,1)<i,1,'last'),3);
            if isempty(label)
                label = 0;
            end
            % refresh textbox position and size
            place_textbox(f2_handle,b2_handle,label)
        end

        if params.show_timeseries
            % update timeseries graph
            for j=1:Nsour
                set(p1_handle(j), 'YData', b_m(:,j));
                set(s1_handle(j), 'Ylim', [-200 200]);
            end
        end
        
        if params.show_local
            % reset local limits
            X_min=Inf; X_max=-Inf; Y_min=Inf; Y_max=-Inf; Z_min=Inf; Z_max=-Inf;
        end
        
        % update 3d graph
        for j=1:Ntraj
            % 'auto-cov A'
            ac_A = b_ry(S_pairs(j,1),S_pairs(j,1),:);
            set(p2_handle(j),'XData',ac_A);
            % 'auto-cov B'
            ac_B = b_ry(S_pairs(j,2),S_pairs(j,2),:);
            set(p2_handle(j),'YData',ac_B);
            % 'cross-cov AB'
            cc_AB = b_ry(S_pairs(j,1),S_pairs(j,2),:);
            set(p2_handle(j),'ZData',cc_AB);
            
            % coordinates of last sample
            ac_A_l = ac_A(end);
            ac_B_l = ac_B(end);
            cc_AB_l = cc_AB(end);

            if params.show_details
                % indices on last sample
                set(tx2_handle(j),'Position',[ac_A_l ac_B_l cc_AB_l]);
                set(tx2_handle(j),'String',num2str(S_pairs(j,:)));
            else
                % spot on last sample
                set(sc2_handle(j),'XData',ac_A_l);
                set(sc2_handle(j),'YData',ac_B_l);
                set(sc2_handle(j),'ZData',cc_AB_l);
            end

            if params.show_local
                % calculate local limits
                if ac_A_l < X_min;  X_min = ac_A_l;  end
                if ac_A_l > X_max;  X_max = ac_A_l;  end
                if ac_B_l < Y_min;  Y_min = ac_B_l;  end
                if ac_B_l > Y_max;  Y_max = ac_B_l;  end
                if cc_AB_l < Z_min; Z_min = cc_AB_l; end
                if cc_AB_l > Z_max; Z_max = cc_AB_l; end
            end
        end
     
        if params.fixed_axes
            % reset plot3 axes limits
            reset_plot3_lims(s2_handle,params);
        end
                
        if params.show_local
            % update local limits
            set(s2_handle,'XLim',[X_min-1 X_max+1]);
            set(s2_handle,'YLim',[Y_min-1 Y_max+1]);
            set(s2_handle,'ZLim',[Z_min-1 Z_max+1]);
        end

%         % update slider position
%         place_slider(f2_handle,sl2_handle);
%         % update values from sliders
%         params.timescale = round(get(sl2_handle,'value'));
        
        % wait fake sampling time (account for timescale)
        pause(((n_rate/params.Fs) - toc) / params.timescale);
        
        drawnow

        % reset fake sampling time
        tic
    end
    
    if params.do_loop && (i == (Nsamp-1))
        i = 1;
    end

end

close all


function place_textbox(f_handle,b_handle,value)
value = num2str(value);
% font ratio to smallest figure dimension
font_ratio = 0.10;
fig_p = get(f_handle,'Position');
fig_w = fig_p(3);
fig_h = fig_p(4);
% add 1 to avoid crush on degenarate resizing
font_size = min((fig_w*font_ratio)+1, (fig_h*font_ratio)+1);
% calculate box size and position
box_w = length(value)*font_size;
box_h = 1.5*font_size;
box_x = fig_w/3 - box_w;
box_y = fig_h - box_h;
% set properties
set(b_handle,'Position',[box_x box_y box_w box_h]);
% update string
set(b_handle,'String',value,'FontSize',font_size)

function reset_plot3_lims(handle,params)
% reset axes limits
set(handle,'XLim',[0 params.ac_lim]);
set(handle,'YLim',[0 params.ac_lim]);
set(handle,'ZLim',[-params.ac_lim params.ac_lim]);

% function place_slider(figure_handle,slider_handle)
% fig_p = get(figure_handle,'Position');
% fig_w = fig_p(3);
% fig_h = fig_p(4);
% box_w = fig_w*(0.78);
% box_h = fig_h*(0.02);
% box_x = fig_w*(0.13);
% box_y = fig_h*(0.03);
% set(slider_handle,'Position',[box_x box_y box_w box_h]);
