function [isolated_s1_rate,isolated_s2_rate,isolated_s1_spect,isolated_s2_spect,s1_rate,s2_rate,s1_spect,s2_spect,livetime,isolated_livetime,allquiet_livetime] = IsolateS1S2s_IncludeNew(ee_merge,MaxDriftLength)

set(gcf,'Visible','off')
isolated_s1s=[];
isolated_s2s=[];

isolated_s1s_pulse_area=[];
isolated_s2s_pulse_area=[];
listOfS1S2s = find(ee_merge.rqs.pulse_classification<3);
listOfS1s = find(ee_merge.rqs.pulse_classification==1);
listOfS2s = find(ee_merge.rqs.pulse_classification==2);
listOfGS2s = find(ee_merge.rqs.pulse_classification==2 & ee_merge.rqs.pulse_area_phe>160);
listOfGS1s = find(ee_merge.rqs.pulse_classification==1 & ee_merge.rqs.pulse_area_phe<50);
listOfGS1S2s = sort(vertcat(listOfGS1s,listOfGS2s));
listOfS3s = find(ee_merge.rqs.pulse_classification==3);
listOfS4s = find(ee_merge.rqs.pulse_classification==4);
s1_rate = length(listOfS1s);
s2_rate = length(listOfS2s);
s1_pulse_area = ee_merge.rqs.pulse_area_phe(listOfS1s);
s2_pulse_area = ee_merge.rqs.pulse_area_phe(listOfS2s);
s1_spect= histogram(s1_pulse_area,(0:1:1000));
s1_spect = s1_spect.Values;
s2_spect= histogram(s2_pulse_area,(0:10:10000));
s2_spect = s2_spect.Values;
Times = diff(ee_merge.rqs.pulse_start_samples);
TimesGTD = Times(Times>2*MaxDriftLength);
allquiet_livetime = sum(TimesGTD)*10e-9;

TimesS1S2 = diff(ee_merge.rqs.pulse_start_samples(listOfS1S2s));
TimesGTDS1S2 = TimesS1S2(TimesS1S2>2*MaxDriftLength);
isolated_livetime = sum(TimesGTDS1S2) * 10e-9;

	for jj = 1:(length(listOfS1S2s)-1)
		delta_t_p1 =MaxDriftLength + 1;
		delta_t_m1 =MaxDriftLength - 1;
		if jj ==1
			delta_t_p1 =  ee_merge.rqs.pulse_start_samples(listOfS1S2s(jj+1)) -  ee_merge.rqs.pulse_start_samples(listOfS1S2s(jj));
	
		elseif jj ==length( ee_merge.rqs.pulse_start_samples)
			delta_t_m1 =  ee_merge.rqs.pulse_start_samples(listOfS1S2s(jj)) -  ee_merge.rqs.pulse_start_samples(listOfS1S2s(jj-1));
		else
			delta_t_p1 =  ee_merge.rqs.pulse_start_samples(listOfS1S2s(jj+1)) -  ee_merge.rqs.pulse_start_samples(listOfS1S2s(jj));
			delta_t_m1 =  ee_merge.rqs.pulse_start_samples(listOfS1S2s(jj)) -  ee_merge.rqs.pulse_start_samples(listOfS1S2s(jj-1));
			if ((delta_t_p1 >MaxDriftLength) && (delta_t_m1 >MaxDriftLength))
				if  (ee_merge.rqs.pulse_classification(listOfS1S2s(jj)) == 1)
					isolated_s1s = [isolated_s1s, jj];
					isolated_s1s_pulse_area=[isolated_s1s_pulse_area, ee_merge.rqs.pulse_area_phe(jj)];
				end
				if  (ee_merge.rqs.pulse_classification(listOfS1S2s(jj)) == 2)

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
	isolated_s1_spect = dummy.Values;
	dummy = histogram(isolated_s2s_pulse_area,(0:10:10000));
	isolated_s2_spect = dummy.Values;
end
