function timing = DatPFC_PulseTiming_HeightTiming_Converted_AreaTiming(pulse_data_phe)
    timing = nan(10,1);
    threshold = 0.15;
    last = length(pulse_data_phe)-1;
    
    [~, max_idx_temp] = max(pulse_data_phe);
    max_idx = max_idx_temp - 1;
    %This is a special cumulative area
    
    pre_above = vertcat(1,pulse_data_phe(1:end-1) > threshold);
    above = pulse_data_phe > threshold;
    pos_above = vertcat(pulse_data_phe(2:end) > threshold,1);

    cumulative_area = cumsum(pulse_data_phe.*( (pre_above | above | pos_above) & (pulse_data_phe >= 0) ));    
    area = sum(pulse_data_phe.*( (pre_above | above | pos_above) & (pulse_data_phe >= 0) ));    
    
    if area <= 0
        timing(5) = max_idx;
        if max_idx > 0
            timing(1) = max_idx-1;
            timing(2) = max_idx-1;
            timing(3) = max_idx-1;
            timing(4) = max_idx-1;
        else
            timing(1) = max_idx;
            timing(2) = max_idx;
            timing(3) = max_idx;
            timing(4) = max_idx;
        end
        if max_idx < last
            timing(6) = max_idx+1;
            timing(7) = max_idx+1;
            timing(8) = max_idx+1;
            timing(9) = max_idx+1;
        else
            timing(6) = max_idx;
            timing(7) = max_idx;
            timing(8) = max_idx;
            timing(9) = max_idx;
        end
        timing(10) = cumulative_area(timing(9)+1) - cumulative_area(timing(1)+1);
    end
    
    cumulative_fraction_area = cumulative_area/area;
    t01_idx_temp = find(cumulative_fraction_area <= 0.01 ,1,'last');
    if ~isempty(t01_idx_temp)
        t01_idx = t01_idx_temp - 1;
    else
        t01_idx = 0;
    end
    t05_idx_temp = find(cumulative_fraction_area <= 0.05 ,1,'last');
    if ~isempty(t05_idx_temp)
        t05_idx = t05_idx_temp - 1;
    else
        t05_idx = 0;
    end
    t10_idx_temp = find(cumulative_fraction_area <= 0.10 ,1,'last');
    if ~isempty(t10_idx_temp)
        t10_idx = t10_idx_temp - 1;
    else
        t10_idx = 0;
    end
    t25_idx_temp = find(cumulative_fraction_area <= 0.25 ,1,'last');
    if ~isempty(t25_idx_temp)
        t25_idx = t25_idx_temp - 1;
    else
        t25_idx = 0;
    end
    
    t50_idx_temp = find(cumulative_fraction_area <= 0.50 ,1,'last');
    if ~isempty(t50_idx_temp)
        t50_idx = t50_idx_temp - 1;
        if t50_idx == last
            t50_idx = last - 1;
        end
    else
        t50_idx = 0;
    end
    
    t75_idx_temp = find(cumulative_fraction_area >= 0.75 ,1,'first');
    if ~isempty(t75_idx_temp)
        t75_idx = t75_idx_temp - 1;
    else
        t75_idx = last;
    end
    t90_idx_temp = find(cumulative_fraction_area >= 0.90 ,1,'first');
    if ~isempty(t90_idx_temp)
        t90_idx = t90_idx_temp - 1;
    else
        t90_idx = last;
    end
    t95_idx_temp = find(cumulative_fraction_area >= 0.95 ,1,'first');
    if ~isempty(t95_idx_temp)
        t95_idx = t95_idx_temp - 1;
    else
        t95_idx = last;
    end
    t99_idx_temp = find(cumulative_fraction_area >= 0.99 ,1,'first');
    if ~isempty(t99_idx_temp)
        t99_idx = t99_idx_temp - 1;
    else
        t99_idx = last;
    end
    
    timing(1) = t01_idx;
    timing(2) = t05_idx;
    timing(3) = t10_idx;
    timing(4) = t25_idx;
    timing(5) = t50_idx;
    timing(6) = t75_idx;
    timing(7) = t90_idx;
    timing(8) = t95_idx;
    timing(9) = t99_idx;
    timing(10) = cumulative_area(timing(9)+1) - cumulative_area(timing(1)+1);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end