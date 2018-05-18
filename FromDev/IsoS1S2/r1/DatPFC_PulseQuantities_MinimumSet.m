function ee = DatPFC_PulseQuantities_MinimumSet(ee,XML_Settings,Debug)
% Excerp from PulseQuantities_MinimumSet.m
if ~exist('Debug','var')
    Debug = 0;
end

%% Initialize variables
    pmt_chs = 1:122;

    pulse_event_size = [XML_Settings.max_num_pulses 1];
    per_channel_event_size = [XML_Settings.max_num_pulses length(pmt_chs) 1];

    % should make these sparse!

    % Per sum
    ee.rqs.pulse_area_phe                       = zeros(pulse_event_size);
    ee.rqs.skinny_pulse_area_phe                = zeros(pulse_event_size);
    ee.rqs.pulse_area_positive_phe              = zeros(pulse_event_size);
    ee.rqs.pulse_area_negative_phe              = zeros(pulse_event_size);
    ee.rqs.pulse_height_phe_per_sample          = zeros(pulse_event_size);
    ee.rqs.prompt_fraction                      = zeros(pulse_event_size);
    ee.rqs.prompt_fraction_tlx                  = zeros(pulse_event_size);
    ee.rqs.top_bottom_ratio                     = zeros(pulse_event_size);
    ee.rqs.top_bottom_asymmetry                 = zeros(pulse_event_size);

    ee.rqs.exp_fit_amplitude_phe_per_sample     = zeros(pulse_event_size);
    ee.rqs.exp_fit_tau_fall_samples             = zeros(pulse_event_size);
    ee.rqs.exp_fit_time_offset_samples          = zeros(pulse_event_size);
    ee.rqs.exp_fit_tau_rise_samples             = zeros(pulse_event_size);
    ee.rqs.exp_fit_chisq                        = zeros(pulse_event_size);
    ee.rqs.exp_fit_dof                          = zeros(pulse_event_size);

    ee.rqs.gaus_fit_amplitude_phe_per_sample    = zeros(pulse_event_size);
    ee.rqs.gaus_fit_mu_samples                  = zeros(pulse_event_size);
    ee.rqs.gaus_fit_sigma_samples               = zeros(pulse_event_size);
    ee.rqs.gaus_fit_chisq                       = zeros(pulse_event_size);
    ee.rqs.gaus_fit_dof                         = zeros(pulse_event_size);

    ee.rqs.s2filter_max_area_diff               = zeros(pulse_event_size);
    ee.rqs.s2filter_max_s2_area                 = zeros(pulse_event_size);
    ee.rqs.s2filter_max_s1_area                 = zeros(pulse_event_size);

    ee.rqs.pulse_mean_samples                       = zeros(pulse_event_size);
    ee.rqs.rms_width_samples                        = zeros(pulse_event_size);
    ee.rqs.pulse_skewness                           = zeros(pulse_event_size);
    ee.rqs.pulse_kurtosis                           = zeros(pulse_event_size);

    ee.rqs.rms_width_samples                    = zeros(pulse_event_size);

    ee.rqs.pulse_length_samples                 = zeros(pulse_event_size);

    preBins     = XML_Settings.prompt_fraction.preBins_samples;
    windowBins  = XML_Settings.prompt_fraction.windowBins_samples;

    s1window_samples = XML_Settings.s2filter.s1window_samples;
    s2window_samples = XML_Settings.s2filter.s2window_samples;

    % Per channel
    ee.rqs.peak_area_phe                        = zeros(per_channel_event_size);
    ee.rqs.peak_area_negative_phe               = zeros(per_channel_event_size);
    ee.rqs.peak_height_phe_per_sample           = zeros(per_channel_event_size);
    ee.rqs.mean_first_last_pts_phe_per_sample   = zeros(per_channel_event_size);
    ee.rqs.peak_height_mV                       = zeros(per_channel_event_size);
    ee.rqs.daq_saturation_flag                  = zeros(per_channel_event_size);
    ee.rqs.pmt_2pct_saturation_flag             = zeros(per_channel_event_size);
    ee.rqs.peak_detect_pretrigger_mean_phe      = zeros(per_channel_event_size);
    ee.rqs.peak_detect_posttrigger_mean_phe     = zeros(per_channel_event_size);
    ee.rqs.peak_detect_pretrigger_std_phe       = zeros(per_channel_event_size);
    ee.rqs.peak_detect_posttrigger_std_phe      = zeros(per_channel_event_size);
    ee.rqs.baseline_daq_mV                      = zeros(per_channel_event_size);

    % new RQs --pfs
    ee.rqs.peak_area_thr_phe                    = zeros(per_channel_event_size);
    ee.rqs.skinny_peak_height_mV				= zeros(per_channel_event_size);
    ee.rqs.skinny_peak_area_phe					= zeros(per_channel_event_size);

    % RQs originally in TRC  - SS
    ee.rqs.pre_pulse_area_positive_phe          = zeros(pulse_event_size);
    ee.rqs.pre_pulse_area_negative_phe          = zeros(pulse_event_size);
    ee.rqs.post_pulse_area_positive_phe         = zeros(pulse_event_size);
    ee.rqs.post_pulse_area_negative_phe         = zeros(pulse_event_size);
    ee.rqs.amis1_fraction                       = zeros(pulse_event_size);

    ee.rqs.pulse_length_sumpods                 = zeros(pulse_event_size);

%try

    %% Compute RQs
    % Check if fields exist and has data
    if isfield(ee,'sumpod_data_phe_per_sample') && ee.empty == 0
        loop_max = min([sum(isfinite(ee.rqs.index_kept_sumpods(:,1))) XML_Settings.max_num_pulses]);
        if loop_max == 0
            if Debug
                disp('Fail to run pulse basic quantities.')
            end
            ee.info.PulseQuantitySuccess = 0;
            ee.info.PulseQuantityError = 'No pulse found';
            return
        end
        % for every pulse (NOT pod)
        for pp = 1:loop_max
            %dis('pp=%ee.rqs,1=%ee.rqs ... index_kept=%ee.rqs\n',pp,1,ee.rqs.index_kept_sumpods(pp,1))
            % -1 and +1 signs are needed to get the right time cut
            %pulse_cut = inrange(ee.sumpod_time_samples,ee.rqs.pulse_start_samples(pp,1)-1,ee.rqs.pulse_end_samples(pp,1)+1);

            if ee.rqs.index_kept_sumpods(pp,1) == 0
                continue
            end

            % use explicit range cuts, not "inrange" function
            pulse_cut = (ee.sumpod_time_samples >= ee.rqs.pulse_start_samples(pp,1)) ...
                      & (ee.sumpod_time_samples <= ee.rqs.pulse_end_samples(pp,1));
            pulse_data_phe = ee.sumpod_data_phe_per_sample(pulse_cut);
            pulse_data_phe_full = ee.sumpod_data_phe_per_sample(ee.sumpod_time_samples < ee.info.EventWindow(2)+1e5);
            pulse_time_samples = ee.sumpod_time_samples(pulse_cut);
            pulse_time_samples_full = ee.sumpod_time_samples(ee.sumpod_time_samples < ee.info.EventWindow(2)+1e5);
            pulse_length = ee.rqs.pulse_end_samples(pp,1) - ee.rqs.pulse_start_samples(pp,1)+1;
             % ***Calculate PULSE quantities per pulse here***

            % Pulse area -- see below, it is calculated from peak_area_phe 
    %            ee.rqs.pulse_area_phe(pp,1) = sum(pulse_data_phe);
    %            ee.rqs.pulse_area_thr_phe(pp,1) = sum(pulse_data_phe(pulse_data_phe>XML_Settings.thr_phe));
    %            ee.rqs.pulse_area_positive_phe(pp,1) = sum(pulse_data_phe(pulse_data_phe >= 0));
    %            ee.rqs.pulse_area_negative_phe(pp,1) = sum(pulse_data_phe(pulse_data_phe < 0));
    %Skinny Boxcar
            fullBoxArea = 0;
            for ii = 1:pulse_length
            if ii <= numel(pulse_data_phe)
              fullBoxArea = fullBoxArea + pulse_data_phe(ii);
            else 
              break
            end
            end
            % Initiliaze values
            maxSkinnyBoxArea = 0;
            tt_maxSkinnyBoxArea = 0;

            % Find skinny box with maximum area
            for tt = 1:(pulse_length-XML_Settings.skinnyBoxSamples+1);
            skinnyBoxArea = 0;
            for ii=tt:tt+XML_Settings.skinnyBoxSamples-1
              if ii <= numel(pulse_data_phe)
                skinnyBoxArea = skinnyBoxArea + pulse_data_phe(ii);
              else 
                break
              end
            end
              if skinnyBoxArea > maxSkinnyBoxArea
                maxSkinnyBoxArea = skinnyBoxArea;
                tt_maxSkinnyBoxArea = tt;
              end 
            end

            % When pulse is shorter than the skinny box:
            if XML_Settings.skinnyBoxSamples > pulse_length
             tt_maxSkinnyBoxArea = 1;
             maxSkinnyBoxArea = fullBoxArea;     
            end

            % Record values
            ee.rqs.skinny_pulse_start_samples(pp,1) = ee.rqs.pulse_start_samples(pp,1) + tt_maxSkinnyBoxArea-1;
            ee.rqs.skinny_pulse_end_samples(pp,1)   = ee.rqs.pulse_start_samples(pp,1) + tt_maxSkinnyBoxArea + XML_Settings.skinnyBoxSamples-1;

            % amis1 - measure of S1-ness 0=>S2, 1=>S1 (ish)
            ee.rqs.amis1_fraction(pp,1) = maxSkinnyBoxArea/fullBoxArea;
            if isnan(ee.rqs.amis1_fraction(pp,1))
            ee.rqs.amis1_fraction(pp,1) = 0; % reasonable asymptote
            end
            if isinf(ee.rqs.amis1_fraction(pp,1))
            ee.rqs.amis1_fraction(pp,1) = -1; % reasonable asymptote 
            end
            if ee.rqs.amis1_fraction(pp,1) > 10 || ee.rqs.amis1_fraction(pp,1) < -10
            ee.rqs.amis1_fraction(pp,1) = -2; % error code
            end


            % Pulse height
            [ee.rqs.pulse_height_phe_per_sample(pp,1) loc] = max(pulse_data_phe);

            ee.rqs.pulse_height_timing_samples(pp,1) = loc;

            % Prompt fraction
            tleft = find(pulse_time_samples == (ee.rqs.aft_t0_samples(pp,1))); % I added the -1, based on the event profiles. we may need to change this as algorithms evolve ..pfs
            tlx = find(pulse_time_samples == (ee.rqs.aft_tlx_samples(pp,1)));
            ee.rqs.prompt_fraction(pp,1) = LUXPromptFraction_framework(pulse_data_phe,preBins,windowBins,tleft);
            ee.rqs.prompt_fraction_tlx(pp,1) = LUXPromptFraction_framework(pulse_data_phe,preBins,windowBins,tlx);

            try
            % Exponential fits
            if isfield(XML_Settings,'exponential_fit')
                efit_settings = XML_Settings.exponential_fit;
            else
                efit_settings = [];
            end
            [efit_params, ee.rqs.exp_fit_chisq(pp,1) ee.rqs.exp_fit_dof(pp,1)] = LUXExponentialFit_framework(1:length(pulse_data_phe),pulse_data_phe,efit_settings);
            catch;efit_params=[0 0 0 0];end; % kludge since I don't understand the error (yet)

            ee.rqs.exp_fit_amplitude_phe_per_sample(pp,1) = efit_params(1);
            ee.rqs.exp_fit_tau_fall_samples(pp,1) = efit_params(2);
            ee.rqs.exp_fit_time_offset_samples(pp,1) = efit_params(3); % this is relative to the start of the PULSE!
            ee.rqs.exp_fit_tau_rise_samples(pp,1) = efit_params(4);
            clear efit_params;

            % Gaussian fits
            if isfield(XML_Settings,'gaussian_fit')
                gfit_settings = XML_Settings.gaussian_fit;
            else
                gfit_settings = [];
            end
            [gfit_params, ee.rqs.gaus_fit_chisq(pp,1) ee.rqs.gaus_fit_dof(pp,1)] = LUXGaussianFit_framework(1:length(pulse_data_phe),pulse_data_phe,gfit_settings);

            ee.rqs.gaus_fit_amplitude_phe_per_sample(pp,1) = gfit_params(1);
            ee.rqs.gaus_fit_mu_samples(pp,1) = gfit_params(2);
            ee.rqs.gaus_fit_sigma_samples(pp,1) = gfit_params(3);
            clear gfit_params;

            % S2Filter
            ext_pulse = zeros(1,length(pulse_data_phe) + 2*s2window_samples); % extend pulse so that filters are not longer than data
            ext_pulse((s2window_samples+1) : (s2window_samples+length(pulse_data_phe))) = pulse_data_phe;

            [max_area_diff max_s2_area max_s1_area] = LUXS2FilterMatlab_framework(ext_pulse,s2window_samples,s1window_samples);

            ee.rqs.s2filter_max_area_diff(pp,1) = max_area_diff;
            ee.rqs.s2filter_max_s2_area(pp,1) = max_s2_area;
            ee.rqs.s2filter_max_s1_area(pp,1) = max_s1_area;

            % rms_width_samples
            rms_width_range = inrange(pulse_time_samples,[ee.rqs.aft_t0_samples(pp,1), ee.rqs.aft_t2_samples(pp,1)]');
            ee.rqs.rms_width_samples(pp,1) = sum(pulse_data_phe(rms_width_range).*(pulse_time_samples(rms_width_range) - ee.rqs.aft_t1_samples(pp,1)).^2) ./ sum(pulse_data_phe(rms_width_range));
            ee.rqs.pulse_mean_samples(pp,1) = sum(pulse_data_phe(rms_width_range).*pulse_time_samples(rms_width_range))./ sum(pulse_data_phe(rms_width_range));

            if(ee.rqs.rms_width_samples(pp,1)>0)
                ee.rqs.rms_width_samples(pp,1)=sqrt(ee.rqs.rms_width_samples(pp,1));
            else
                ee.rqs.rms_width_samples(pp,1)=-1.0*sqrt(-1.0*ee.rqs.rms_width_samples(pp,1));
            end
            ee.rqs.pulse_skewness(pp,1) = sum(pulse_data_phe(rms_width_range).*(pulse_time_samples(rms_width_range) - ee.rqs.pulse_mean_samples(pp,1)).^3) ./ sum(pulse_data_phe(rms_width_range));
            ee.rqs.pulse_skewness(pp,1)=ee.rqs.pulse_skewness(pp,1)/(ee.rqs.rms_width_samples(pp,1)^3);
            ee.rqs.pulse_kurtosis(pp,1) = sum(pulse_data_phe(rms_width_range).*(pulse_time_samples(rms_width_range) - ee.rqs.pulse_mean_samples(pp,1)).^4) ./ sum(pulse_data_phe(rms_width_range));
            ee.rqs.pulse_kurtosis(pp,1)=ee.rqs.pulse_kurtosis(pp,1)/(ee.rqs.rms_width_samples(pp,1)^4)-3.0;



            %Summing of pre and post pulse PMT signals

            tmp_pre_t0 = max(pulse_time_samples_full(1), ee.rqs.pulse_start_samples(pp,1)-XML_Settings.preBoxSamples);  
            tmp_pre_range = find((pulse_time_samples_full >= tmp_pre_t0) & (pulse_time_samples_full <= (ee.rqs.pulse_start_samples(pp,1)-1)));

            tmp_post_t2 = min(pulse_time_samples_full(end),ee.rqs.pulse_end_samples(pp,1) + XML_Settings.preBoxSamples);
            tmp_post_range = find((pulse_time_samples_full >= (ee.rqs.pulse_end_samples(pp,1)+1)) & (pulse_time_samples_full<=tmp_post_t2));    

            tmp_pre_pulse_data_phe = pulse_data_phe_full(tmp_pre_range);
            tmp_post_pulse_data_phe = pulse_data_phe_full(tmp_post_range);

            ee.rqs.pre_pulse_area_positive_phe(pp,1) = sum(tmp_pre_pulse_data_phe(tmp_pre_pulse_data_phe>=0));
            ee.rqs.pre_pulse_area_negative_phe(pp,1) = sum(tmp_pre_pulse_data_phe(tmp_pre_pulse_data_phe<0));
            ee.rqs.post_pulse_area_positive_phe(pp,1) = sum(tmp_post_pulse_data_phe(tmp_post_pulse_data_phe>=0));
            ee.rqs.post_pulse_area_negative_phe(pp,1) = sum(tmp_post_pulse_data_phe(tmp_post_pulse_data_phe<0));

            % If the following rq is shorter than pulse_length_samples there are missing sumpods in pulse_data_phe
            ee.rqs.pulse_length_sumpods(pp,1) = length(pulse_data_phe);  

            total_sum_peak_area_phe = 0;

            % PER CH QUANTITIES
            ch_map = find(~[ee.ch(:).empty]);

            if ~isempty(ch_map)

                for ch = ch_map

                        % Make slice of data
                        peak_cut = inrange(ee.ch(ch).pod_time_samples,ee.rqs.pulse_start_samples(pp,1)-1,ee.rqs.pulse_end_samples(pp,1)+1);

                    if ~isempty(peak_cut) && sum(peak_cut) > 0

                        peak_data_phe = ee.ch(ch).pod_data_phe_per_sample(peak_cut);
                        peak_time_samples = ee.ch(ch).pod_time_samples(peak_cut);
                        peak_data_mV = ee.ch(ch).pod_data_mV(permute(peak_cut,[2 1]));

                        % This is to find which baseline_mV to use. Since
                        % we are cutting pods in time, it's not
                        % straightforward to know which pod baseline 'twas                            
                        % But this statement takes care of that
                        peak_inds = unique(cumsum(ismember(ee.ch(ch).pod_time_samples,ee.ch(ch).pod_start_samples)));

                        ee.rqs.baseline_daq_mV(pp,ch,1) = mean(ee.ch(ch).pod_baseline_mV(peak_inds));

                        % Peak area
                        ee.rqs.peak_area_phe(pp,ch,1) = sum(peak_data_phe);
                        ee.rqs.peak_area_negative_phe(pp,ch,1) = sum(peak_data_phe(peak_data_phe<0));
                        ee.rqs.peak_area_thr_phe(pp,ch,1) = sum(peak_data_phe(peak_data_phe>XML_Settings.thr_phe));

                        % new RQ: skinnyBox peak area
                        skinny_peak_cut = inrange(ee.ch(ch).pod_time_samples,ee.rqs.skinny_pulse_start_samples(pp,1)-1,ee.rqs.skinny_pulse_end_samples(pp,1)+1);				
                        ee.rqs.skinny_peak_area_phe(pp,ch,1) = sum(ee.ch(ch).pod_data_phe_per_sample(skinny_peak_cut));
                        tmp = max(ee.ch(ch).pod_data_mV(skinny_peak_cut));
                        if isempty(tmp); tmp = NaN; end;
                        ee.rqs.skinny_peak_height_mV(pp,ch,1) = tmp;

                        % For debugging
                        total_sum_peak_area_phe = total_sum_peak_area_phe + sum(peak_data_phe);

                        % Peak height
                        ee.rqs.peak_height_phe_per_sample(pp,ch,1) = max(peak_data_phe);
                        ee.rqs.peak_height_mV(pp,ch,1) = max(peak_data_mV);

                        ee.rqs.mean_first_last_pts_phe_per_sample(pp,ch,1) = mean(peak_data_phe([1 end]));

            % Pre and post trigger baseline rqs. Corrected to use POD array instead of peak array - AL - 2013/08/02
                        if length(ee.ch(ch).pod_data_phe_per_sample) >= 24 % added to fix occational crash - JRV - 2013-03-04
                            ee.rqs.peak_detect_pretrigger_mean_phe(pp,ch,1) = mean(ee.ch(ch).pod_data_phe_per_sample(1:24));
                            ee.rqs.peak_detect_pretrigger_std_phe(pp,ch,1)  = std(ee.ch(ch).pod_data_phe_per_sample(1:24));
                        end

                        if length(ee.ch(ch).pod_data_phe_per_sample) >= 31 % added to fix occational crash - JRV - 2013-03-04
                            ee.rqs.peak_detect_posttrigger_mean_phe(pp,ch,1) = mean(ee.ch(ch).pod_data_phe_per_sample((end-31+1):end));
                            ee.rqs.peak_detect_posttrigger_std_phe(pp,ch,1) = std(ee.ch(ch).pod_data_phe_per_sample((end-31+1):end));
                        end

                        % Saturation flags
                        ee.rqs.daq_saturation_flag(pp,ch,1) = ee.rqs.peak_height_mV(pp,ch,1) > XML_Settings.saturation_flags.daq_saturation_threshold_mV;
                        ee.rqs.pmt_2pct_saturation_flag(pp,ch,1) = ee.rqs.peak_height_mV(pp,ch,1) > XML_Settings.saturation_flags.pmt_saturation_2pct_threshold_mV;
                    end % fi


                end % for ch

                ee.rqs.pulse_area_phe(pp,1) = squeeze(sum(ee.rqs.peak_area_phe(pp,:,1),2));
                ee.rqs.skinny_pulse_area_phe(pp,1) = squeeze(sum(ee.rqs.skinny_peak_area_phe(pp,:,1),2));

                % For debugging sum(peak_area_phe) ~= pulse_area_phe
                if 0
                    if ~inrange(total_sum_peak_area_phe./ee.rqs.pulse_area_phe(pp,1),0.99,1.001)

                        figure(43423); clf;
                        h1 = subplot(2,1,1); hold on
                        plot(pulse_time_samples,pulse_data_phe,'k.-')
                        plot(ee.sumpod_time_samples,ee.sumpod_data_phe_per_sample,'ro-')
                        title(sum(pulse_data_phe));

                        h2 = subplot(2,1,2); hold on

                        for ch = 1:122
                            peak_cut = inrange(ee.ch(ch).pod_time_samples,ee.rqs.pulse_start_samples(pp,1)-1,ee.rqs.pulse_end_samples(pp,1)+1);
                            %if sum(peak_cut); %{fprintf('%ee.rqs ',ch);}% end
                            peak_data_phe = ee.ch(ch).pod_data_phe_per_sample(peak_cut);
                            peak_time_samples = ee.ch(ch).pod_time_samples(peak_cut);
                            plot(peak_time_samples,peak_data_phe,'k.-')
                            plot(ee.ch(ch).pod_time_samples,ee.ch(ch).pod_data_phe_per_sample,'ro-')
                        end
                        title(total_sum_peak_area_phe);
                        linkaxes([h1 h2],'xy');
                        xlim([ee.rqs.pulse_start_samples(pp,1)-10 ee.rqs.pulse_end_samples(pp,1)+10])


                    end
                end


            end % fi            

            top = sum(squeeze(ee.rqs.peak_area_phe(pp,[1:60 121],1)),2);
            bot = sum(squeeze(ee.rqs.peak_area_phe(pp,[61:120 122],1)),2);

            ee.rqs.top_bottom_ratio(pp,1) =  (top)./(bot);
            ee.rqs.top_bottom_asymmetry(pp,1) =  (top-bot)./(top+bot);

        end % for pulse pp
        ee.info.PulseQuantitySuccess = 1;
        ee.info.PulseQuantityError = '';
    else
        if Debug
            disp('Fail to run pulse basic quantities.')
        end
        ee.info.PulseQuantitySuccess = 0;
        ee.info.PulseQuantityError = ee.info.SumPODError;
    end % fi sumpod
%{    
catch exception
    if Debug
        disp('Fail to run pulse basic quantities.')
    end
    ee.info.PulseQuantitySuccess = 0;
    ee.info.PulseQuantityError = exception.identifier;
    
end
%}
