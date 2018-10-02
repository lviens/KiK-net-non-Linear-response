function  [dat_s] = KiKnet_borehole_response(dir_fold, station, bpcut , order, trust, plt)

cp = 1;
for ii = 1: length(dir_fold)
    listHZ =  [dir([dir_fold(ii).folder '/'  dir_fold(ii).name '/' station '*EW*']); ...
        dir([dir_fold(ii).folder  '/'  dir_fold(ii).name '/' station '*NS*'])];
    listUD =  dir([dir_fold(ii).folder  '/'  dir_fold(ii).name '/' station '*UD*']);
    
    cut_time = 10;
    cut_time2 = cut_time + 20 ;
    
    for iii =1:length(listUD)
        filename = [listUD(iii).folder '/' listUD(iii).name ];
        [h, t, delta] = read_KiK_net(filename);
        h = h-mean(h);
        h = BP_butter_Filter(h,order,delta,bpcut(1),bpcut(2));
        dataUD(:,iii) = h;
        [loc(iii), snr_db(iii)] = Simple_PphasePicker(dataUD(cut_time*delta :cut_time2*delta, iii), 1/delta);
    end
        
    if  (abs(loc(1) - loc(2)))<1  && loc(2)>1 &&loc(1)>1 && delta==100 && length(dataUD(:,2))>10000
        mean_loc = loc(2)+cut_time- 1;
        
        if trust ==1
            for ik = 1:length(listHZ)
                filename = [listHZ(ik).folder '/' listHZ(ik).name ];
                [h_f,~, ~] = read_KiK_net(filename);
                h_f = BP_butter_Filter(h_f, order, delta, bpcut(1), bpcut(2));
                dat_s.dataHZ{cp}(:,ik) = h_f;
            end
            
            dat_s.Keep_event(cp) = ii;
            dat_s.mean_loc_f(cp) = mean_loc*delta;
            dat_s.len_sign(cp) = length(dataUD(:,2));
            cp = cp+1;
        else
            figure
            subplot(3,1,1)
            plot( t(1:length(dataUD(cut_time*delta:cut_time2*delta,1)))+cut_time, dataUD(cut_time*delta:cut_time2*delta,1) )
            hold on
            plot(loc(iii)+cut_time,0, 'xr')
            xlim([cut_time cut_time2])
            title( [num2str(ii) ' / ' num2str(length(dir_fold))] )
            
            subplot(3,1,2)
            plot( t(1:length(dataUD(cut_time*delta:cut_time2*delta,2))) + cut_time, dataUD(cut_time*delta:cut_time2*delta,2) )
            hold on
            plot(loc(1)+cut_time,0, 'xr')
            xlim([cut_time cut_time2])
            subplot(3,1,3)
            plot(t, dataUD(:,2))
            hold on
            plot(loc(2) + cut_time,0, 'xr')
            while(1)
                m=input('Press [n] if you want to skip the event, enter or any other key to accept: ','s');
                if m=='n'
                    break
                else
                    for ik = 1:length(listHZ)
                        filename = [listHZ(ik).folder '/' listHZ(ik).name ];
                        [h_f,~, ~] = read_KiK_net( filename );
                        h_f = BP_butter_Filter( h_f, order, delta, bpcut(1), bpcut(2) );
                        dat_s.dataHZ{cp}(:,ik) = h_f;
                    end
                    
                    dat_s.Keep_event(cp) = ii;
                    dat_s.mean_loc_f(cp) = mean_loc*delta;
                    dat_s.len_sign(cp) = length(dataUD(:,2));
                    cp = cp+1;
                    break
                end
            end
            close all
        end
    end
    dataUD = [];
    h = [];  
end

dif_len = min(dat_s.len_sign - dat_s.mean_loc_f);


for ij = 1:length(dat_s.Keep_event)
    
    lim_ini = round(dat_s.mean_loc_f(ij));
    
    hann = hanning( .5*delta );
    dataf{ij}(:,iii) = zeros(500*delta, 1);
    for iii=1:length(listHZ)
        dataf{ij}(1:dif_len+1,iii) = dat_s.dataHZ{ij}( lim_ini:lim_ini+dif_len, iii ) ;
        dataf{ij}(1:length(hann)/2,iii) = dataf{ij}(1:length(hann)/2,iii).*hann(1:end/2);
    end
    datacompl{ij}(:,1) = complex( dataf{ij}(:,1), dataf{ij}(:,3) );
    datacompl{ij}(:,2) = complex( dataf{ij}(:,2), dataf{ij}(:,4) );
    [ dataftt{ij}(:,1), ~] = FFT_amplitude( datacompl{ij}(:,1), delta ) ;
    [ dataftt{ij}(:,2), f] = FFT_amplitude( datacompl{ij}(:,2), delta ) ;
    for iii =1 :2
        [dataftt{ij}(:,iii)] = liskonno(f, dataftt{ij}(:,iii), 100);
    end
    dat_s.ratio(ij,:) = abs(dataftt{ij}(:,2))./abs(dataftt{ij}(:,1));
    dat_s.f = f;
    
    if plt ==1
        semilogx(f, dat_s.ratio(ij,:))
        hold on
        xlim([.5 30])
        title([num2str(ij) ' / ' num2str(length(dat_s.Keep_event)) ])
        xlabel('Frequency (Hz)')
        ylabel('Amplification')
        grid on
%         pause(.1)
    end
    
end