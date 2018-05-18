function ee = DatPFC_PulseFinder_TRC(ee,XML_Settings,SR_Version,Debug)

if ~exist('SR_Version','var')
    disp('SR_Version is not specified. Use default SR2.0.')
    SR_Version = 'SR2.0';
end  
if ~exist('Debug','var')
    Debug = 0;
end
%try
    if strcmp(SR_Version,'SR2.1')
        %% Excerp from PulseFinder_TransparentRubiksCube.m
        %% Initialize variables
        pulse_event_size = [XML_Settings.max_num_pulses 1];
        sum_thresh = XML_Settings.noiseThre;

        %optional_extra_iterations = 3; % 3 is the magic number (yes it is)

        % 5 lines copied from PulseFinder_Naive.m
        ee.rqs.pulse_start_samples          = nan(pulse_event_size);
        ee.rqs.pulse_end_samples            = nan(pulse_event_size);
        ee.rqs.index_kept_sumpods           = nan(pulse_event_size);	
        ee.rqs.num_pulses_found             = 0;

        %% Compute timing RQs
        clear range pulse_time_samples* pulse_data_phe*;

        % Check if fields exist and has data
        if isfield(ee,'sumpod_data_phe_per_sample') && ee.empty == 0
            pulse_time_samples = ee.sumpod_time_samples(ee.sumpod_time_samples<(XML_Settings.EventWindow(end) + 1e5));
            pulse_time_samples_full = ee.sumpod_time_samples(ee.sumpod_time_samples<(XML_Settings.EventWindow(end) + 1e5));
            pulse_data_phe = ee.sumpod_data_phe_per_sample(ee.sumpod_time_samples<(XML_Settings.EventWindow(end) + 1e5));

            ee.rqs.n_samples_in_evt     = length(pulse_time_samples);
            ee.rqs.full_evt_area_phe    = sum(pulse_data_phe);

            done_searching_region_a = false;
                    pp = 1;
            fullBoxSamples = XML_Settings.fullBoxSamples;

            while pp <= XML_Settings.max_num_pulses % loop over PULSES
                %TextProgress(pp,XML_Settings.max_num_pulses,'Count')
        %			if pp>1 & flags(1)<mymodule_settings.skinnyBoxSamples % the pulse just found was determined to have a width smaller than the skinny box
        %				fullBoxSamples = mymodule_settings.skinnyBoxSamples % so search again, using the skinnyBox as the fullBox
        %			end
                if pp==3 && ~exist('range','var') % after finding the first two pulses, define search regions a and b
                    %disp('** beginning search of region_a')
                    range = pulse_time_samples<max(ee.rqs.pulse_start_samples(1,1),ee.rqs.pulse_start_samples(2,1)); % min or max, a minor quandry
                    % divide up the pulse so you only look before the 2nd (in time) of the two largest pulses found
                    pulse_time_samples_a = pulse_time_samples(range);
                    pulse_data_phe_a = pulse_data_phe(range);
                    pulse_time_samples_b = pulse_time_samples(~range);
                    pulse_data_phe_b = pulse_data_phe(~range);

                    pulse_time_samples = pulse_time_samples_a;
                    pulse_data_phe = pulse_data_phe_a;
                end

                % DO NOT REMOVE FOLLOWING LINE -pfs
                if length(pulse_time_samples)<10; break; end; % need a cutoff, this one is a bit arbitrary, should check it...
                % perusePeeks looks for an S2 or S2-like pulse	
                % perusePeeks returns the actual sample times, not the indices to the times!

                [afTiming, pp_areas] = ...
                        perusePeeks_sr21( pulse_time_samples ,pulse_data_phe ,[...
                                    fullBoxSamples ...
                                    XML_Settings.maximumGap ...
                                    XML_Settings.nLookAhead ...
                                    XML_Settings.nLookBehind ...
                                    XML_Settings.noiseThre ...
                                    ] );

                % next check if the pulse should be thrown out
                % any pulse with only 3 samples are the result of fluctuations following S2's
                if (afTiming(2) - afTiming(1)) > 3 ...
                   || (afTiming(2) ~= afTiming(1)+1) && (pp_areas(1)/(afTiming(2) - afTiming(1))) > XML_Settings.noiseThre
                    keep_pulse = 1;
                else
                    keep_pulse = 0;
                end
                if keep_pulse
                    % next line needed for PulseQuantities (note it is NOT related to pods)
                    ee.rqs.index_kept_sumpods(pp,1) = 1; 

                    % next, we need the "true" pulse edges in time. this is the region that will be 
                    % integrated to give the pulse area
                    ee.rqs.pulse_start_samples(pp,1) 			= double((afTiming(1)));
                    ee.rqs.pulse_end_samples(pp,1) 			= double((afTiming(2)));


                else
                    %disp('going to toss pulse...');%keyboard;
                end

                % zero part of waveform where the pulse was found
                % we HAVE to use the afTiming array for this info, as dp.pulse_start_samples and dp.pulse_end_samples
                % might not be filled
                zeroMe = (pulse_time_samples>=afTiming(1)) ...
                        &(pulse_time_samples<=afTiming(2));
                pulse_data_phe(zeroMe) = 0;		
                if ~done_searching_region_a
                    pulse_data_phe_a(zeroMe) = 0;
                end

                % if we have not exhausted XML_Settings.max_num_pulses, go and search the region *after* the 2nd (in time) largest pulse
                if (pp_areas(1)<sum_thresh) && ~done_searching_region_a
                %if (~keep_pulse) & ~done_searching_region_a
                    done_searching_region_a = true;
                    if exist('pulse_time_samples_b','var')
                        pulse_time_samples = pulse_time_samples_b;
                        pulse_data_phe = pulse_data_phe_b;
                        %disp('** beginning search of region_b')
                    else
                        break; % nothing left to find
                    end
                elseif ~(keep_pulse) && done_searching_region_a
                    break; % nothing of consequence left to find
                end

                %if (pp==XML_Settings.max_num_pulses) 
                %	break; % call it a day
                %end
                        % make sure to increment the pulse number if appropriate
                if keep_pulse
                    pp = pp + 1;
                end

            if pp == XML_Settings.max_num_pulses+1 && pp_areas(1) > 0 && sum(pulse_data_phe_a) < 80 && ~done_searching_region_a  %if last pulse is >0phe and remaining event area is less than bad area cut
                    priority = 0;                                                                                        % also if last pulse is in region a
                    mightbeSE = 0;
                    mightbeSPE = 0;
                    if afTiming(2)-afTiming(1)+1 > 50 && pp_areas(1) >= 5 % looks like SE 
                        mightbeSE = 1;
                    end

                   ch_map = find(~[ee.ch(:).empty]);
                   ch_map = intersect(ch_map,find(~[ee.ch(:).empty])); 
                   peak_height_phe_per_sample = zeros(1,122); % initialise each time to prevent data from previous pulses staying in PMT channels

                   if ~isempty(ch_map)
                        for ch = ch_map
                            peak_cut = inrange(ee.ch(ch).pod_time_samples,afTiming(1)-1,afTiming(2)+1);  

                            if ~isempty(peak_cut) && sum(peak_cut) > 0

                                peak_data_phe = ee.ch(ch).pod_data_phe_per_sample(peak_cut); 
                                peak_height_phe_per_sample(1,ch) = max(peak_data_phe);
                            end
                        end 
                        if sum(peak_height_phe_per_sample(1,:)>0.09) == 1 %only one PMT reaches threshold
                          mightbeSPE = 1;
                        end 
                        if mightbeSE ==1 || mightbeSPE ==1 % if last pulse looks like an SPE or SE look again
                            while priority == 0 
                            %set pulse time samples and data back to region a

                                [afTiming, pp_areas] = ...
                                perusePeeks_sr21( pulse_time_samples_a, pulse_data_phe_a ,[...
                                    fullBoxSamples ...
                                    XML_Settings.maximumGap ...
                                    XML_Settings.nLookAhead ...
                                    XML_Settings.nLookBehind ...
                                    XML_Settings.noiseThre ...
                                    ] );
                                if (pp_areas(1)<sum_thresh)
                                    break %stop looking if nothing left to find
                                end 

                                % now check for 2PMTs
                                peak_height_phe_per_sample = zeros(1,122);
                                for ch = ch_map
                                    peak_cut = inrange(ee.ch(ch).pod_time_samples,afTiming(1),afTiming(2));    
                                    if ~isempty(peak_cut) && sum(peak_cut) > 0
                                        peak_data_phe = ee.ch(ch).pod_data_phe_per_sample(peak_cut); 
                                        peak_height_phe_per_sample(1,ch) = max(peak_data_phe);
                                    end
                                end
                                if sum(peak_height_phe_per_sample(1,:)>0.09) > 1 && (afTiming(2)-afTiming(1)+1) < 100 %length cut will remove most SEs
                                    priority = 1;  %prioritise this pulse as it could be an S1
                                else 
                                    zeroMe = (pulse_time_samples_a>=afTiming(1))&(pulse_time_samples_a<=afTiming(2));
                                    pulse_data_phe_a(zeroMe) = 0;	 % pulse not an S1, zero it
                                end

                            end
                        end


                        if  priority == 1  %now do usual quality check on pulse
                            if (afTiming(2) - afTiming(1)) > 3 ...
                            || (afTiming(2) ~= afTiming(1)+1) && (pp_areas(1)/(afTiming(2) - afTiming(1))) > XML_Settings.noiseThre
                                ee.rqs.pulse_start_samples(pp-1,1) 			= double((afTiming(1))); %reset pulse start and end times of last pulse to priotised one
                                ee.rqs.pulse_end_samples(pp-1,1)              = double((afTiming(2)));			
                            end
                        end
                    end
                end

            end % for pulse pp

            ee.rqs.num_pulses_found = sum(ee.rqs.index_kept_sumpods==1); % convenience RQ
            % Sort all pulse-level quantities by time (previously variables to be 
            % sorted were listed explicitly, this caused a crash if event-level 
            % varibales were added but not included in the list. Instead we now
            % only sort dp varibales with identical dimensions to pulse_event_size)                
            fields_all = fieldnames(ee.rqs);
            [Y,I]=sort(ee.rqs.pulse_start_samples(:,1),1,'ascend');
            for f=1:length(fieldnames(ee.rqs))
                if isequal(size(ee.rqs.(fields_all{f})), pulse_event_size)
                    ee.rqs.(fields_all{f})(:,1) = ee.rqs.(fields_all{f})(I,1);
                end 
            end % for fieldnames


            % Buffer for SPE tails (set in xml settings)
            pMax = ee.rqs.num_pulses_found(1);
            if pMax > 0 

                for p = 1:pMax-1

                    if (ee.rqs.pulse_end_samples(p,1)+ XML_Settings.extendPulse < ee.rqs.pulse_start_samples(p+1,1))
                        ee.rqs.pulse_end_samples(p,1) = ee.rqs.pulse_end_samples(p,1) + XML_Settings.extendPulse;
                    end
                    if (ee.rqs.pulse_end_samples(p,1)+ XML_Settings.extendPulse >= ee.rqs.pulse_start_samples(p+1,1))
                        ee.rqs.pulse_end_samples(p,1) = ee.rqs.pulse_start_samples(p+1,1)-1;
                    end

                end

                if (ee.rqs.pulse_end_samples(pMax,1) + XML_Settings.extendPulse <= max(pulse_time_samples_full))
                    ee.rqs.pulse_end_samples(pMax,1) = ee.rqs.pulse_end_samples(pMax,1) + XML_Settings.extendPulse;
                end
                if (ee.rqs.pulse_end_samples(pMax,1) + XML_Settings.extendPulse > max(pulse_time_samples_full))
                    ee.rqs.pulse_end_samples(pMax,1) = max(pulse_time_samples_full);
                end

            end
            
            
            ee.rqs.pulse_start_samples          = ee.rqs.pulse_start_samples(1:ee.rqs.num_pulses_found);
            ee.rqs.pulse_end_samples            = ee.rqs.pulse_end_samples(1:ee.rqs.num_pulses_found);
            ee.rqs.index_kept_sumpods           = ee.rqs.index_kept_sumpods(1:ee.rqs.num_pulses_found);

            ee.info.PulseFinderSuccess = 1;
            ee.info.PulseFinderError = '';
            ee.info.PulseFinderVersion = SR_Version;
        else
            if Debug
                disp('Fail to run pulse finder.')
            end
            ee.rqs.pulse_start_samples          = nan(1,1);
            ee.rqs.pulse_end_samples            = nan(1,1);
            ee.rqs.index_kept_sumpods           = nan(1,1);	
            ee.rqs.num_pulses_found             = 0;

            ee.info.PulseFinderSuccess = 0;
            ee.info.PulseFinderError = ee.info.SumPODError;
        end % fi
      

    elseif strcmp(SR_Version,'SR2.0')
        %% Excerp from PulseFinder_TransparentRubiksCube.m
        %% Initialize variables
        pulse_event_size = [XML_Settings.max_num_pulses 1];
        sum_thresh = XML_Settings.noiseThre;

        %optional_extra_iterations = 3; % 3 is the magic number (yes it is)

        % 5 lines copied from PulseFinder_Naive.m
        ee.rqs.pulse_start_samples          = nan(pulse_event_size);
        ee.rqs.pulse_end_samples            = nan(pulse_event_size);
        ee.rqs.index_kept_sumpods           = nan(pulse_event_size);	
        ee.rqs.num_pulses_found             = 0;
        
        %% Compute timing RQs
        clear range pulse_time_samples* pulse_data_phe*;

        % Check if fields exist and has data
        if isfield(ee,'sumpod_data_phe_per_sample') && ee.empty == 0
            pulse_time_samples = ee.sumpod_time_samples(ee.sumpod_time_samples<(XML_Settings.EventWindow(end) + 1e5));
            pulse_time_samples_full = ee.sumpod_time_samples(ee.sumpod_time_samples<(XML_Settings.EventWindow(end) + 1e5));
            pulse_data_phe = ee.sumpod_data_phe_per_sample(ee.sumpod_time_samples<(XML_Settings.EventWindow(end) + 1e5));

            ee.rqs.n_samples_in_evt     = length(pulse_time_samples);
            ee.rqs.full_evt_area_phe    = sum(pulse_data_phe);

            done_searching_region_a = false;
                    pp = 1;
            fullBoxSamples = XML_Settings.fullBoxSamples;

            while pp <= XML_Settings.max_num_pulses % loop over PULSES
        %			if pp>1 & flags(1)<mymodule_settings.skinnyBoxSamples % the pulse just found was determined to have a width smaller than the skinny box
        %				fullBoxSamples = mymodule_settings.skinnyBoxSamples % so search again, using the skinnyBox as the fullBox
        %			end
                if pp==3 && ~exist('range','var') % after finding the first two pulses, define search regions a and b
                    %disp('** beginning search of region_a')
                    range = pulse_time_samples<max(ee.rqs.pulse_start_samples(1,1),ee.rqs.pulse_start_samples(2,1)); % min or max, a minor quandry
                    % divide up the pulse so you only look before the 2nd (in time) of the two largest pulses found
                    pulse_time_samples_a = pulse_time_samples(range);
                    pulse_data_phe_a = pulse_data_phe(range);
                    pulse_time_samples_b = pulse_time_samples(~range);
                    pulse_data_phe_b = pulse_data_phe(~range);

                    pulse_time_samples = pulse_time_samples_a;
                    pulse_data_phe = pulse_data_phe_a;
                end

                % DO NOT REMOVE FOLLOWING LINE -pfs
                if length(pulse_time_samples)<10; break; end; % need a cutoff, this one is a bit arbitrary, should check it...
                % perusePeeks looks for an S2 or S2-like pulse	
                % perusePeeks returns the actual sample times, not the indices to the times!

                [afTiming, pp_areas] = ...
                        perusePeeks_sr20( pulse_time_samples ,pulse_data_phe ,[...
                                    XML_Settings.preBoxSamples ...
                                    fullBoxSamples ...
                                    XML_Settings.postBoxSamples ...
                                    XML_Settings.edgeFraction ...
                                    XML_Settings.txFraction ...
                                    XML_Settings.skinnyBoxSamples ...
                                    XML_Settings.maximumGap ...
                                    XML_Settings.nLookAhead ...
                                    XML_Settings.nLookBehind ...
                                    XML_Settings.noiseThre ...
                                    ] );

                % next check if the pulse should be thrown out
                % any pulse with only 3 samples are the result of fluctuations following S2's
                if (afTiming(2) - afTiming(1)) > 3 ...
                   || (afTiming(2) ~= afTiming(1)+1) && (pp_areas(2)/(afTiming(2) - afTiming(1))) > XML_Settings.noiseThre
                    keep_pulse = 1;
                else
                    keep_pulse = 0;
                end
                if keep_pulse
                    % next line needed for PulseQuantities (note it is NOT related to pods)
                    ee.rqs.index_kept_sumpods(pp,1) = 1; 

                    % next, we need the "true" pulse edges in time. this is the region that will be 
                    % integrated to give the pulse area
                    ee.rqs.pulse_start_samples(pp,1) 			= double((afTiming(1)));
                    ee.rqs.pulse_end_samples(pp,1) 			= double((afTiming(2)));


                else
                    %disp('going to toss pulse...');%keyboard;
                end

                % zero part of waveform where the pulse was found
                % we HAVE to use the afTiming array for this info, as dp.pulse_start_samples and dp.pulse_end_samples
                % might not be filled
                zeroMe = (pulse_time_samples>=afTiming(1)) ...
                        &(pulse_time_samples<=afTiming(2));
                pulse_data_phe(zeroMe) = 0;		

                % if we have not exhausted XML_Settings.max_num_pulses, go and search the region *after* the 2nd (in time) largest pulse
                if (pp_areas(2)<sum_thresh) && ~done_searching_region_a
                %if (~keep_pulse) & ~done_searching_region_a
                    done_searching_region_a = true;
                    if exist('pulse_time_samples_b','var')
                        pulse_time_samples = pulse_time_samples_b;
                        pulse_data_phe = pulse_data_phe_b;
                        %disp('** beginning search of region_b')
                    else
                        break; % nothing left to find
                    end
                elseif ~(keep_pulse) && done_searching_region_a
                    break; % nothing of consequence left to find
                end

                %if (pp==XML_Settings.max_num_pulses) 
                %	break; % call it a day
                %end
                        % make sure to increment the pulse number if appropriate
                if keep_pulse
                    pp = pp + 1;
                end
            end % for pulse pp

            ee.rqs.num_pulses_found = sum(ee.rqs.index_kept_sumpods==1); % convenience RQ
            % Sort all pulse-level quantities by time (previously variables to be 
            % sorted were listed explicitly, this caused a crash if event-level 
            % varibales were added but not included in the list. Instead we now
            % only sort dp varibales with identical dimensions to pulse_event_size)                
            fields_all = fieldnames(ee.rqs);
            [Y,I]=sort(ee.rqs.pulse_start_samples(:,1),1,'ascend');
            for f=1:length(fieldnames(ee.rqs))
                if isequal(size(ee.rqs.(fields_all{f})), pulse_event_size)
                    ee.rqs.(fields_all{f})(:,1) = ee.rqs.(fields_all{f})(I,1);
                end 
            end % for fieldnames


            % Buffer for SPE tails (set in xml settings)
            pMax = ee.rqs.num_pulses_found(1);
            if pMax > 0 

                for p = 1:pMax-1

                    if (ee.rqs.pulse_end_samples(p,1)+ XML_Settings.extendPulse < ee.rqs.pulse_start_samples(p+1,1))
                        ee.rqs.pulse_end_samples(p,1) = ee.rqs.pulse_end_samples(p,1) + XML_Settings.extendPulse;
                    end
                    if (ee.rqs.pulse_end_samples(p,1)+ XML_Settings.extendPulse >= ee.rqs.pulse_start_samples(p+1,1))
                        ee.rqs.pulse_end_samples(p,1) = ee.rqs.pulse_start_samples(p+1,1)-1;
                    end

                end

                if (ee.rqs.pulse_end_samples(pMax,1) + XML_Settings.extendPulse <= max(pulse_time_samples_full))
                    ee.rqs.pulse_end_samples(pMax,1) = ee.rqs.pulse_end_samples(pMax,1) + XML_Settings.extendPulse;
                end
                if (ee.rqs.pulse_end_samples(pMax,1) + XML_Settings.extendPulse > max(pulse_time_samples_full))
                    ee.rqs.pulse_end_samples(pMax,1) = max(pulse_time_samples_full);
                end

            end
            ee.rqs.pulse_start_samples          = ee.rqs.pulse_start_samples(1:ee.rqs.num_pulses_found);
            ee.rqs.pulse_end_samples            = ee.rqs.pulse_end_samples(1:ee.rqs.num_pulses_found);
            ee.rqs.index_kept_sumpods           = ee.rqs.index_kept_sumpods(1:ee.rqs.num_pulses_found);

            ee.info.PulseFinderSuccess = 1;
            ee.info.PulseFinderError = '';
            ee.info.PulseFinderVersion = SR_Version;
        else
            if Debug
                disp('Fail to run pulse finder.')
            end
            ee.rqs.pulse_start_samples          = nan(1,1);
            ee.rqs.pulse_end_samples            = nan(1,1);
            ee.rqs.index_kept_sumpods           = nan(1,1);	
            ee.rqs.num_pulses_found             = 0;

            ee.info.PulseFinderSuccess = 0;
            ee.info.PulseFinderError = ee.info.SumPODError;
        end % fi
        
    else
        error('Do not recognize SR version.');
    end
%{
catch exception
    if Debug
        disp('Fail to run pulse finder.')
    end
    ee.rqs.pulse_start_samples          = nan(1,1);
    ee.rqs.pulse_end_samples            = nan(1,1);
    ee.rqs.index_kept_sumpods           = nan(1,1);	
    ee.rqs.num_pulses_found             = 0;
    
    ee.info.PulseFinderSuccess = 0;
    ee.info.PulseFinderVersion = SR_Version;
    if isfield(ee.info,'SumPODSuccess') && isfield(ee.info,'SumPODError')
        if ~ee.info.SumPODSuccess
            ee.info.PulseFinderError = ee.info.SumPODError;
        else
            ee.info.PulseFinderError = exception.identifier;
        end
    else
        ee.info.PulseFinderError = exception.identifier;
    end
end
%}

    
    
