function ee = DatPFC_PulseClassifier_MultiDimensional(ee,PCM_Settings,SR_Version,Debug)
if ~exist('SR_Version','var')
    disp('SR_Version is not specified. Use default SR2.0.')
    SR_Version = 'SR2.0';
end 
if ~exist('Debug','var')
    Debug = 0;
end
try
    if strcmp(SR_Version,'SR2.0')
        %% Initialize variables

        pmt_chs = 1:122;

        [A B] = size(ee.rqs.pulse_area_phe);

        
        pulse_event_size = [PCM_Settings.max_num_pulses 1];
        if PCM_Settings.max_num_pulses == 0
            ee.rqs.pulse_classification           = NaN(pulse_event_size); % initialize
            if Debug
                disp('Fail to run pulse classifier.')
            end
            ee.info.PulseClassifierSuccess = 0;
            ee.info.PulseClassifierError = 'No pulse found';
            ee.info.PulseClassifierVersion = SR_Version;
            return
        end
        
        %% Define Categories %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %  This is where you define the cuts for S1 (cut cS1), S2 (cS2),
        %  sphe (cut_sphe), SE (cut_SE) and other (cut_nota). 
        %
        %


        % no pulse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cut_no_pulse = ee.rqs.pulse_area_phe == 0;



        % S1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % band 1: double-boxcar

        bc_min_S1 = -0.01+(-0.5*exp(-1.2*(ee.rqs.pulse_area_phe)));

        % modify to allow for 83mKr S1's
        %bc_max_S1 = 0.055 + 9e-2*(log(ee.rqs.pulse_area_phe).^-1);  bc_max_S1(ee.rqs.pulse_area_phe<1) = 1;  bc_max_S1(bc_max_S1>1) = 1;
        %bc_max_S1 = 0.075 + 6.5e-2*(log(ee.rqs.pulse_area_phe).^-1.4);  bc_max_S1(ee.rqs.pulse_area_phe<1) = 1;  bc_max_S1(bc_max_S1>1) = 1;
        bc_max_S1 = 0.07 + 9e-2*(log(ee.rqs.pulse_area_phe).^-0.9);  bc_max_S1(ee.rqs.pulse_area_phe<1) = 1;  bc_max_S1(bc_max_S1>1) = 1;
        bc_max_S1(inrange(ee.rqs.pulse_area_phe,[100,500])) = .4;


        cS1bc =  (ee.rqs.s2filter_max_area_diff./ee.rqs.pulse_area_phe > bc_min_S1) & ...
                 (ee.rqs.s2filter_max_area_diff./ee.rqs.pulse_area_phe < bc_max_S1) ;



        % band 2: prompt fraction at 10% (tlx)

        %pf_min_S1_1 = 0.7+(-1*exp(-0.27*(ee.rqs.pulse_area_phe-1)));
        %pf_min_S1_1 = 0.62+(-1*exp(-0.33*(ee.rqs.pulse_area_phe-1)));
        pf_min_S1_1 = 0.56+(-1.2*exp(-0.26*(ee.rqs.pulse_area_phe+0.2)));
        %pf_min_S1_2 = 0.75*(1-0.8*sigmf(log10(ee.rqs.pulse_area_phe), [3 2]));
        %pf_min_S1_2 = 0.68*(1-0.8*sigmf(log10(ee.rqs.pulse_area_phe), [3 2]));
        pf_min_S1_2 = 0.68*(1-0.8*sigmf(log10(ee.rqs.pulse_area_phe), [2.6 2]));
        %pf_energy = 18.3;
        %pf_energy = 20.5;
        pf_energy = 32.8;

        cS1pf_1 = ((ee.rqs.prompt_fraction_tlx >  pf_min_S1_1) & (ee.rqs.pulse_area_phe<pf_energy));
        cS1pf_2 = ((ee.rqs.prompt_fraction_tlx >  pf_min_S1_2) & (ee.rqs.pulse_area_phe>=pf_energy));


        cS1pf = cS1pf_1 | cS1pf_2;

        % band 3: top-bottom asymmetry 
        z_min_S1 = -0.55 - (0.5*log(ee.rqs.pulse_area_phe)).^-0.7;  z_min_S1(z_min_S1 < -1.1) = -1.1;  z_min_S1(ee.rqs.pulse_area_phe < 1) = -1.1;
        z_max_S1 = -0.35 + (0.3*log(ee.rqs.pulse_area_phe)).^-2.2;  z_max_S1(z_max_S1 > 1.1) = 1.1; z_max_S1(ee.rqs.pulse_area_phe < 1) = 1.1;


        cS1z  = (ee.rqs.top_bottom_asymmetry > z_min_S1) & ...
                (ee.rqs.top_bottom_asymmetry < z_max_S1) ;


        %band 4: ratios of width LR

        %w_max_S1 = exp(-0.035*ee.rqs.pulse_area_phe)+0.28-8e-1*(1./(exp(-1*(log10(ee.rqs.pulse_area_phe)-6))));
        w_max_S1 = exp(-0.025*ee.rqs.pulse_area_phe)+0.29-8e-1*(1./(exp(-1*(log10(ee.rqs.pulse_area_phe)-6))));
        w_max_S1(inrange(ee.rqs.pulse_area_phe, [150, 500])) = .34;

        cS1w = double(ee.rqs.aft_t1_samples-ee.rqs.aft_t0_samples)./double(ee.rqs.aft_t2_samples-ee.rqs.aft_t0_samples)<w_max_S1;


        %additional cuts LR

        singlephe_peak_threshold=0.09;
        cut_which_tubes_hit_height = (ee.rqs.peak_height_phe_per_sample > singlephe_peak_threshold);
        numtubes_height = squeeze(sum(cut_which_tubes_hit_height,2));

        cut_two_tubes_height = (numtubes_height>=2);

        %if there's only one event in the evt some arrays me be transposed, this
        %messes up the combining the cuts
        %if size(cut_two_tubes_height,2) == 1
        if size(cut_two_tubes_height,2) ~= 1 %Changed for DatPFC
            cut_two_tubes_height=cut_two_tubes_height';
        end

        skinny_area_threshold=0.3;
        cut_which_tubes_hit_sarea = (ee.rqs.skinny_peak_area_phe > skinny_area_threshold);
        numtubes_sarea = squeeze(sum(cut_which_tubes_hit_sarea,2));

        cut_two_tubes_sarea = (numtubes_sarea>=2);

        %if there's only one event in the evt some arrays me be transposed, this
        %messes up the combining the cuts
        %if size(cut_two_tubes_sarea,2) == 1
        if size(cut_two_tubes_sarea,2) ~= 1 %Changed for DatPFC
        cut_two_tubes_sarea=cut_two_tubes_sarea';
        end

        cut_two_tubes = cut_two_tubes_height & cut_two_tubes_sarea;


        % combine all four bands
        cS1 = cS1bc & cS1pf & cS1w & cS1z;
        cS1a = cS1bc & cS1pf & cS1w & cS1z & cut_two_tubes;



        % single phe
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        cut_one_tube = (numtubes_height == 1);


        %if there's only one event in the evt some arrays me be transposed, this
        %messes up the combining the cuts
        %if size(cut_one_tube,2) == 1
        if size(cut_one_tube,2) ~= 1 %Change for DatPFC
            cut_one_tube=cut_one_tube';
        end

        cut_sphe = cut_one_tube & cS1;




        % S2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %band 1: double boxcar

        bc_min_S2 = 0.12 - (1.35*log(ee.rqs.pulse_area_phe).^-1.5);  bc_min_S2(bc_min_S2 < 2e-4) = 2e-4;  bc_min_S2(ee.rqs.pulse_area_phe < 1) = 2e-4;
        bc_max_S2 = ones(size(ee.rqs.pulse_area_phe));

        cS2bc =  (ee.rqs.s2filter_max_area_diff./ee.rqs.pulse_area_phe > bc_min_S2) & ...
                 (ee.rqs.s2filter_max_area_diff./ee.rqs.pulse_area_phe < bc_max_S2) ;


        % band 2: prompt fraction


        pf_min_S2 = -0.015+(-0.6*exp(-0.3*(ee.rqs.pulse_area_phe)));
        pf_max_S2 = 0.7*(1-0.9*sigmf(log10(ee.rqs.pulse_area_phe), [3 2]));


        cS2pf_first = (ee.rqs.prompt_fraction >  pf_min_S2) & ...
                (ee.rqs.prompt_fraction <  pf_max_S2) ;


        % alterante cut from Jeremy, modified LR - includes oultier population of
        % S2s, i.e. most of the here selected S2s are multiples scatters clustered together
        %- should be re-tuned post new pulse finding algorithm
        pf_max_S2_alt = 0.2;
        pf_max_S2_alt_energy = 3;

        cS2pf_alt = (ee.rqs.prompt_fraction > pf_min_S2) & (log10(ee.rqs.pulse_area_phe) >= pf_max_S2_alt_energy) & (ee.rqs.prompt_fraction < pf_max_S2_alt);


        cS2pf = cS2pf_first | cS2pf_alt;    

        % band 3: top-bottom asymmetry


        %S2mean = 0.1-0.50*sigmf(log10(ee.rqs.pulse_area_phe), [2 4.8]); z_min_S2 = S2mean - (0.4*log(0.5*ee.rqs.pulse_area_phe)).^-1.5; z_min_S2(z_min_S2 < -1.1) = -1.1; z_min_S2(ee.rqs.pulse_area_phe < 4) = -1.1 ;
        S2mean = 0.05-0.50*sigmf(log10(ee.rqs.pulse_area_phe), [2 4.8]); z_min_S2 = S2mean - (0.4*log(0.5*ee.rqs.pulse_area_phe)).^-1.5; z_min_S2(z_min_S2 < -1.1) = -1.1; z_min_S2(ee.rqs.pulse_area_phe < 4) = -1.1;
        S2mean = 0.15*(1-sigmf(log10(ee.rqs.pulse_area_phe), [2 5]));  z_max_S2 = S2mean + (0.5*log(ee.rqs.pulse_area_phe)).^-0.5;  z_max_S2(z_max_S2 > 1.1) = 1.1; z_max_S2(ee.rqs.pulse_area_phe < 1) = 1.1;


        cS2z  = (ee.rqs.top_bottom_asymmetry > z_min_S2) & ...
                (ee.rqs.top_bottom_asymmetry < z_max_S2) ;


        %additional cuts LR

        %possible width cut?
        %w_min = 10;

        %cs2w = (ee.rqs.gaus_fit_sigma_samples > w_min);



        % in addition to SE, LR
        h_min = 1;
        cS2h = (ee.rqs.pulse_height_phe_per_sample > h_min);

        %S2_min_e = 50;
        %S2_min_e = 100;
        S2_min_e = 33;

        cS2e = (ee.rqs.pulse_area_phe >= S2_min_e);

        % combine all three bands
        cS2 = cS2bc & cS2pf & cS2z;

        cS2a = cS2bc & cS2pf & cS2z & cS2h & cS2e;




        % single electron %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        Single_e_min = 5;
        %Single_e_max = 100;
        Single_e_max = 33;

        cut_se =  cS2 & (ee.rqs.pulse_area_phe >= Single_e_min) & ...
                        (ee.rqs.pulse_area_phe < Single_e_max) ;   




        % neither S1 nor S2 (nor no pulse) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cut_nota = ~cS1a & ~cut_sphe & ~cS2a & ~cut_se & ~cut_no_pulse;





        %% Create output object %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % assign each pulse a numerical category, using this system:

        % 0 means no pulse (empty entry in matrix)
        % 1 means S1
        % 2 means S2
        % 3 means single phe
        % 4 means single electron S2
        % 5 means none of the above (nota)

        ee.rqs.pulse_classification           = nan(pulse_event_size); % initialize

        ee.rqs.pulse_classification(cS2a     ) = 2;   % first, assign both categories of S2
        ee.rqs.pulse_classification(cut_se  ) = 4;   % 

        ee.rqs.pulse_classification(cS1a     ) = 1;   % then, assign S1 (since low-energy events are in both cuts,
        ee.rqs.pulse_classification(cut_sphe) = 3;   %                  are indistinguishible, and are more properly
                                                 %                  called S1 than S2)

        ee.rqs.pulse_classification(cut_nota) = 5;   % (order of this one shouldn't matter)

        ee.rqs.pulse_classification(cut_no_pulse) = 0; %last, just in case the S1 and S2 definitions include zero-area pulses.
        
        ee.info.PulseClassifierSuccess = 1;
        ee.info.PulseClassifierError = '';
        ee.info.PulseClassifierVersion = SR_Version;
    elseif strcmp(SR_Version,'SR2.1')
       
        %% Initialize variables

        pmt_chs = 1:122;

        [A B] = size(ee.rqs.pulse_area_phe);

        pulse_event_size = [PCM_Settings.max_num_pulses 1];
        if PCM_Settings.max_num_pulses == 0
            ee.rqs.pulse_classification           = NaN(pulse_event_size); % initialize
            if Debug
                disp('Fail to run pulse classifier.')
            end
            ee.info.PulseClassifierSuccess = 0;
            ee.info.PulseClassifierError = 'No pulse found';
            ee.info.PulseClassifierVersion = SR_Version;
            return
        end
        %% Define Categories %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %  This is where you define the cuts for S1 (cut cS1), S2 (cS2),
        %  sphe (cut_sphe), SE (cut_SE) and other (cut_nota). 
        %
        %


        % no pulse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cut_no_pulse = ee.rqs.pulse_area_phe == 0;



        % S1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % band 1: double-boxcar
        bc_min_S1 = -0.01+(-0.5*exp(-1.2*(ee.rqs.pulse_area_phe)));
        % modify to allow for 83mKr S1's

        bc_max_S1 = 0.07 + 9e-2*(log(ee.rqs.pulse_area_phe).^-0.9); 
        bc_max_S1(ee.rqs.pulse_area_phe<1) = 1;  bc_max_S1(bc_max_S1>1) = 1;
        bc_max_S1(inrange(ee.rqs.pulse_area_phe,[100,500])) = .4;


        cS1bc =  (ee.rqs.s2filter_max_area_diff./ee.rqs.pulse_area_phe > bc_min_S1) & ...
                 (ee.rqs.s2filter_max_area_diff./ee.rqs.pulse_area_phe < bc_max_S1) ;



        % band 2: prompt fraction at 5% (tlx)
        pf_min_S1_1 = 0.56+(-1.2*exp(-0.26*(ee.rqs.pulse_area_phe+0.2)));
        pf_min_S1_2 = 0.68*(1-0.8*sigmf(log10(ee.rqs.pulse_area_phe), [2.6 2]));
        pf_energy = 32.8;

        cS1pf_1 = ((ee.rqs.prompt_fraction_tlx >  pf_min_S1_1) & (ee.rqs.pulse_area_phe<pf_energy));
        cS1pf_2 = ((ee.rqs.prompt_fraction_tlx >  pf_min_S1_2) & (ee.rqs.pulse_area_phe>=pf_energy));


        cS1pf = cS1pf_1 | cS1pf_2;

        % band 3: top-bottom asymmetry 
        z_min_S1 = -0.55 - (0.5*log(ee.rqs.pulse_area_phe)).^-0.7;  z_min_S1(z_min_S1 < -1.1) = -1.1;  z_min_S1(ee.rqs.pulse_area_phe < 1) = -1.1;
        z_max_S1 = -0.35 + (0.3*log(ee.rqs.pulse_area_phe)).^-2.2;  z_max_S1(z_max_S1 > 1.1) = 1.1; z_max_S1(ee.rqs.pulse_area_phe < 1) = 1.1;


        cS1z  = (ee.rqs.top_bottom_asymmetry > z_min_S1) & ...
                (ee.rqs.top_bottom_asymmetry < z_max_S1) ;


        %band 4: ratios of width LR

        w_max_S1 = exp(-0.025*ee.rqs.pulse_area_phe)+0.29-8e-1*(1./(exp(-1*(log10(ee.rqs.pulse_area_phe)-6))));
        w_max_S1(inrange(ee.rqs.pulse_area_phe, [100, 500])) = .34;

        cS1w = (ee.rqs.aft_t1_samples-ee.rqs.aft_t0_samples)./(ee.rqs.aft_t2_samples-ee.rqs.aft_t0_samples)<w_max_S1;


        %additional cuts LR

        %modified LR
        singlephe_peak_threshold=0.09;
        cut_which_tubes_hit_height = (ee.rqs.peak_height_phe_per_sample > singlephe_peak_threshold);
        numtubes_height = squeeze(sum(cut_which_tubes_hit_height,2));

        cut_two_tubes_height = (numtubes_height>=2);


        %if there's only one event in the evt some arrays me be transposed, this
        %messes up the combining the cuts
        if size(cut_two_tubes_height,2) == 1
        %    cut_two_tubes_height=cut_two_tubes_height';
        end

        skinny_area_threshold=0.3;
        cut_which_tubes_hit_sarea = (ee.rqs.skinny_peak_area_phe > skinny_area_threshold);
        numtubes_sarea = squeeze(sum(cut_which_tubes_hit_sarea,2));

        cut_two_tubes_sarea = (numtubes_sarea>=2);

        %if there's only one event in the evt some arrays me be transposed, this
        %messes up the combining the cuts
        if size(cut_two_tubes_sarea,2) == 1
        %    cut_two_tubes_sarea=cut_two_tubes_sarea';
        end

        cut_two_tubes = cut_two_tubes_height & cut_two_tubes_sarea;

        %LR, replace top-bottom_asymmetry with width ratio cut

        % combine all three bands
        cS1 = cS1bc & cS1pf & cS1w & cS1z;

        cS1a = cS1bc & cS1pf & cS1w & cS1z & cut_two_tubes;

        % single phe
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        cut_one_tube = (numtubes_height == 1);
        %cut_singlephe_energy = (ee.rqs.pulse_area_phe < singlephe_pulse_max); 

        %if there's only one event in the evt some arrays me be transposed, this
        %messes up the combining the cuts
        if size(cut_one_tube,2) == 1
        %    cut_one_tube=cut_one_tube';
        end

        cut_sphe = cut_one_tube & cS1;

        % S2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %band 1: double boxcar

        bc_min_S2 = 0.12 - (1.35*log(ee.rqs.pulse_area_phe).^-1.5);  bc_min_S2(bc_min_S2 < 2e-4) = 2e-4;  bc_min_S2(ee.rqs.pulse_area_phe < 1) = 2e-4;
        bc_max_S2 = ones(size(ee.rqs.pulse_area_phe));

        cS2bc =  (ee.rqs.s2filter_max_area_diff./ee.rqs.pulse_area_phe > bc_min_S2) & ...
                 (ee.rqs.s2filter_max_area_diff./ee.rqs.pulse_area_phe < bc_max_S2) ;


        % band 2: prompt fraction

        %pf_min_S2 = 1e-2*(1-0.8*sigmf(log10(ee.rqs.pulse_area_phe), [3 2]));
        %modified LR
        pf_min_S2 = -0.015+(-0.6*exp(-0.3*(ee.rqs.pulse_area_phe)));
        pf_max_S2 = 0.7*(1-0.9*sigmf(log10(ee.rqs.pulse_area_phe), [3 2]));
        %modified - adapted from Jeremy LR
        %pf_max_S2 = 0.7*(1-0.83*sigmf(log10(ee.rqs.pulse_area_phe), [3.5 2]));


        cS2pf_first = (ee.rqs.prompt_fraction >  pf_min_S2) & ...
                (ee.rqs.prompt_fraction <  pf_max_S2) ;


        % alterante cut from Jeremy, modified LR - includes oultier population of
        % S2s, i.e. most of the here sleected S2s are multiples scatters clustered together
        %- should be re-tuned post new pulse finding algorithm
        pf_max_S2_alt = 0.2;
        pf_max_S2_alt_energy = 3;

        cS2pf_alt = (ee.rqs.prompt_fraction > pf_min_S2) & (log10(ee.rqs.pulse_area_phe) >= pf_max_S2_alt_energy) & (ee.rqs.prompt_fraction < pf_max_S2_alt);


        cS2pf = cS2pf_first | cS2pf_alt;    

        % band 3: top-bottom asymmetry

        S2mean = 0.05-0.50*sigmf(log10(ee.rqs.pulse_area_phe), [2 4.8]); z_min_S2 = S2mean - (0.4*log(0.5*ee.rqs.pulse_area_phe)).^-1.5; z_min_S2(z_min_S2 < -1.1) = -1.1; z_min_S2(ee.rqs.pulse_area_phe < 4) = -1.1;
        S2mean = 0.15*(1-sigmf(log10(ee.rqs.pulse_area_phe), [2 5]));  z_max_S2 = S2mean + (0.5*log(ee.rqs.pulse_area_phe)).^-0.5;  z_max_S2(z_max_S2 > 1.1) = 1.1; z_max_S2(ee.rqs.pulse_area_phe < 1) = 1.1;


        cS2z  = (ee.rqs.top_bottom_asymmetry > z_min_S2) & ...
                (ee.rqs.top_bottom_asymmetry < z_max_S2) ;


        %additional cuts LR

        %possible width cut?
        %w_min = 10;

        %cs2w = (ee.rqs.gaus_fit_sigma_samples > w_min);



        % in addition to SE, LR
        h_min = 1;
        cS2h = (ee.rqs.pulse_height_phe_per_sample > h_min);

        %S2_min_e = 50;
        S2_min_e = 55;
        cS2e = (ee.rqs.pulse_area_phe > S2_min_e);

        % cut for gas events SS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t25 = double(ee.rqs.aft_t25_samples-ee.rqs.aft_t0_samples);
        t75 = double(ee.rqs.aft_t75_samples-ee.rqs.aft_t0_samples);
        t99 = double(ee.rqs.aft_t2_samples-ee.rqs.aft_t0_samples);

        cGas1 = t75./(0.75*t99) < 0.8;
        cGas2 = (t75./(0.75*t99)) > 0.8 & (t25./(0.25*t99)) > 1.97*(t75./(0.75*t99)) - 1.6;

        cut_gas = cGas1 | cGas2; 


        % combine all three bands
        cS2 = cS2bc & cS2pf & cS2z;

        cS2a = cS2bc & cS2pf & cS2z & cS2h & cS2e & cut_gas;




        % single electron %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % that doesn't seem right to require a minimum threshold - at least then we
        % to define a minimum threshold for S2s as well!!!!

        Single_e_min = 5;
        %Single_e_max = 70;
        %modified LR
        %Single_e_max = 50; % fit from Jeremy 27+/-7 - mu+- 3 sigma? aparerntly set to 50 currently
        Single_e_max = 55;

        cut_se =  cS2 & (ee.rqs.pulse_area_phe > Single_e_min) & ...
                        (ee.rqs.pulse_area_phe < Single_e_max) ;   




        % neither S1 nor S2 (nor no pulse) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cut_nota = ~cS1a & ~cut_sphe & ~cS2a & ~cut_se & ~cut_no_pulse;
        %think about cut_nota again!! does is need cS1 a and cs2a?

        %potential gas events / merged multiple scatters
        cut_area = ee.rqs.pulse_area_phe > 150;
        cut_merged = ~cut_gas & cut_area & ~cS1a;




        %% Create output object %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % assign each pulse a numerical category, using this system:

        % 0 means no pulse (empty entry in matrix)
        % 1 means S1
        % 2 means S2
        % 3 means single phe
        % 4 means single electron S2
        % 5 means none of the above (nota)

        ee.rqs.pulse_classification           = NaN(pulse_event_size); % initialize

        ee.rqs.pulse_classification(cS2a     ) = 2;   % first, assign both categories of S2
        ee.rqs.pulse_classification(cut_se  ) = 4;   % 

        ee.rqs.pulse_classification(cut_merged) = 9; % indicates possible merged pulses 

        ee.rqs.pulse_classification(cS1a     ) = 1;   % then, assign S1 (since low-energy events are in both cuts,
        ee.rqs.pulse_classification(cut_sphe) = 3;   %                  are indistinguishible, and are more properly
                                                 %                  called S1 than S2)

        ee.rqs.pulse_classification(cut_nota) = 5;   % (order of this one shouldn't matter)

        ee.rqs.pulse_classification(cut_no_pulse) = 0; %last, just in case the S1 and S2 definitions include zero-area pulses.

        ee.info.PulseClassifierSuccess = 1;
        ee.info.PulseClassifierError = '';
        ee.info.PulseClassifierVersion = SR_Version;
    else
        error('Do not recognize SR version.');
    end
catch exception
    ee.rqs.pulse_classification           = NaN(pulse_event_size); % initialize
    if Debug
        disp('Fail to run pulse classifier.')
    end
    
    ee.info.PulseClassifierSuccess = 0;
    ee.info.PulseClassifierError = exception.identifier;
    ee.info.PulseClassifierVersion = SR_Version;
end