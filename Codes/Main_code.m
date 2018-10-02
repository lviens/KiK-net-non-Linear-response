%% L. Viens 10/01/2018
% Nonlinear response of a KiK-net station borehole (Japan) versus its linear response

clear all
close all
clc

directory = '/Users/loic/Desktop/Jupyter_notebook/KiKnet_nonlinearity/Stations' ; % Needs be changed, but the subfolders have to be as described earlier
station = {'IBUH01'};

% Filter parameters
cut   = .05; % in Hz
cut2  = 30;  % in Hz
order = 4;
disp(['Parameters for Butterworth filter: Order: ' num2str(order) ' Filter limits: ' num2str(cut) '-' num2str(cut2) ' Hz'])

%% if more that one station, replace "i = 1" by the following line:
for i =1 :length(station) %(also need to uncomment the line after saving the files)
    
    % Search for the earthquake in the Mainshock folder.
    dir_fold_main = dir([directory '/' station{i} '/Mainshock/' station{i} '*']);
    % Computes the response of the borehole for this station.
    [dat_m] = KiKnet_borehole_response(dir_fold_main, station{i}, [cut cut2], order, 1, 1);
    
    dir_fold_small = dir([directory '/' station{i} '/Small_events/' station{i} '*']);
    dat_s = KiKnet_borehole_response(dir_fold_small, station{i}, [cut cut2], order, 1, 1);
    
    dat = dat_s.ratio;
    [ave, std1, std2] = Geo_average_KiK(dat);
    
    fig = figure;
    semilogx(dat_s.f, ave, 'b', 'linewidth', 2)
    hold on
    
    semilogx(dat_m.f, dat_m.ratio, 'r', 'linewidth', 2)
    
    dat_s.f(1) = 1.0000e-03;
    fill([dat_s.f fliplr(dat_s.f)], [std1(2,:) fliplr(std1(1,:))], [0.7 0.7 .7], 'FaceAlpha', 0.4)
    
    set(gca, 'XScale', 'log')
    xlim([.5 30])
    
    plot(dat_s.f, std1(1,:), 'k')
    semilogx(dat_s.f, std1(2,:), 'k')
    title([station{i} ' station'])
    xlabel('Frequency (Hz)')
    ylabel('Borehole response')
    grid on
    
    leg = legend('Linear response', 'Mainshock')  ;
    set(leg, 'location', 'northeast')
    set(gca, 'XTickLabel', {'1';'10'});
    
    
    data.dat_m = dat_m;
    data.dat_s = dat_s;
    data.average = ave;
    data.std1 = std1;
    data.std2 = std2;
    save([ directory '/' station{i} '/'  station{i} '_Borehole_response.mat'] , 'data')
    
end
