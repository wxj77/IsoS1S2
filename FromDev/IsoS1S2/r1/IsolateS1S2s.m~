function [isolated_s1_rate,isolated_s2_rate,s1_spect,s2_spect,livetime] = IsolateS1S2s(ee_merge,MaxDriftLength)
isolated_s1s=[];
isolated_s2s=[];
isolated_s1s_pulse_area=[];
isolated_s2s_pulse_area=[];


	for jj = 1:length(ee_merge.rqs.pulse_start_samples)
		delta_t_p1 =MaxDriftLength + 1;
		delta_t_m1 =MaxDriftLength + 1;
		if jj ==1
			delta_t_p1 =  ee_merge.rqs.pulse_start_samples(jj+1) -  ee_merge.rqs.pulse_start_samples(jj);
	
		elseif jj ==length( ee_merge.rqs.pulse_start_samples)
			delta_t_m1 =  ee_merge.rqs.pulse_start_samples(jj) -  ee_merge.rqs.pulse_start_samples(jj-1);
		else
			delta_t_p1 =  ee_merge.rqs.pulse_start_samples(jj+1) -  ee_merge.rqs.pulse_start_samples(jj);
			delta_t_m1 =  ee_merge.rqs.pulse_start_samples(jj) -  ee_merge.rqs.pulse_start_samples(jj-1);
			if ((delta_t_p1 >MaxDriftLength) && (delta_t_m1 >MaxDriftLength))
				if  (ee_merge.rqs.pulse_classification(jj) == 1)
					isolated_s1s = [isolated_s1s, jj];
					isolated_s1s_pulse_area=[isolated_s1s_pulse_area, ee_merge.rqs.pulse_area_phe(jj)];
				end
				if  (ee_merge.rqs.pulse_classification(jj) == 2)
					isolated_s2s = [isolated_s2s, jj];
					isolated_s2s_pulse_area=[isolated_s2s_pulse_area, ee_merge.rqs.pulse_area_phe(jj)];
				end
			end
		end
	end
	livetime = (ee_merge.livetime.end - ee_merge.livetime.latch)*10e-9;
	isolated_s1_rate = length(isolated_s1s);
	isolated_s2_rate = length(isolated_s2s);
	dummy = histogram(isolated_s1s_pulse_area,(0:1:1000));
	s1_spect = dummy.Values;
	dummy = histogram(isolated_s2s_pulse_area,(0:10:10000));
	s2_spect = dummy.Values;
end
