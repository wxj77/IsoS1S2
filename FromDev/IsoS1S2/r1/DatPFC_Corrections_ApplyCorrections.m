function ee = DatPFC_Corrections_ApplyCorrections(ee,XML_Settings,Iq,SR_Version,Debug)

if ~exist('SR_Version','var')
    disp('SR_Version is not specified. Use default SR2.0.')
    SR_Version = 'SR2.0';
end 
if ~exist('Debug','var')
    Debug = 0;
end

%% Initialize variables

    pmt_chs = 1:122; % must take from daq settings

    pulse_event_size = [XML_Settings.max_num_pulses 1];

    %livetime = LUXGetLivetime_framework(filename_evt,data_path_evt);
    %livetime = ee.rqs.admin.livetime;

    %  Correcting possible issue with top_bottom_ratio

    %{ 
    %Another dirty edit DatPFC probably doesn't want this

    [a1 b1] = size(ee.rqs.pulse_area_phe);  
    [a2 b2] = size(ee.rqs.top_bottom_ratio);  
    adif = a1-a2;

    blanks = [];

    for i=1:b1
       blanks = [blanks 0]; 
    end

    if adif>0
       for i=1:adif
         ee.rqs.top_bottom_ratio = vertcat(ee.rqs.top_bottom_ratio,blanks);    
       end
    end
    %}
    % Initializing corrected values of area

    ee.rqs.z_corrected_pulse_area_all_phe   = ee.rqs.pulse_area_phe;
    ee.rqs.xyz_corrected_pulse_area_all_phe = ee.rqs.pulse_area_phe;

    ee.rqs.z_corrected_pulse_area_bot_phe   = ee.rqs.pulse_area_phe./(1+ee.rqs.top_bottom_ratio);
    ee.rqs.xyz_corrected_pulse_area_bot_phe = ee.rqs.pulse_area_phe./(1+ee.rqs.top_bottom_ratio);

    ee.rqs.correction_electron_lifetime = nan(pulse_event_size);
    ee.rqs.correction_s1_z_dependence   = nan(pulse_event_size);
    ee.rqs.correction_s1_xyz_dependence = nan(pulse_event_size);
    ee.rqs.correction_s2_xy_dependence  = nan(pulse_event_size);
    % Initialize variable with negative value so that we're sure it goes through the code 
    %{ 
    Really ?
    ee.rqs.s1_correction_is_3d= zeros(1,1,'int8');
    ee.rqs.s1_correction_is_3d(:) = -99;
    %}

disp(XML_Settings)
disp(XML_Settings.Setting.wrapper)
    s1_ref_z_ns = XML_Settings.detector_centre_z_ns;
    allowed_gs  = XML_Settings.allowed_gs;

    
    if XML_Settings.max_num_pulses == 0
        if Debug
            disp('Fail to run pulse area correction.')
        end
        ee.info.PulseCorrSuccess = 0;
        ee.info.PulseCorrError = 'No pulse found.';
        ee.info.PulseCorrVersion = SR_Version;
        return
    end
try

if strcmp(SR_Version,'SR2.0')
    %-----------------------------------------------------------------------------------------------------------------
    %-----------------------------------------------------------------------------------------------------------------
    %% Reading iq values from the IQs xml ----------------------------------------------------------------------------
    %-----------------------------------------------------------------------------------------------------------------
    %----------------------------------------------------------------------------------------------------------------- 
    % electron_lifetime - Interpolate - DONE
    % z_dep_par_all - Nearest - DONE
    % z_dep_par_bot - Nearest - DONE
    % x_bins,y_bins,z_bins,s1_map_all - Nearest
    % x_bins,y_bins,z_bins,s1_map_bot - Nearest
    % x_bins,y_bins,s2_map_all - Nearest
    % x_bins,y_bins,s2_map_bot - Nearest

    %[a num_iqs] = size(lug_iqs_xml.iq);

    %-----------------------------------
    % Finding the electron lifetime
    %-----------------------------------

    lifetime_values = [];
    dataset_times = [];
    filename_prefixs = {};
    cp_number = {};
    %Why do this ???
    %{
     for i=1:num_iqs
      if isfield(PulseCorrIq.electron_lifetime.iq.correction,'fit')
        if isfield(PulseCorrIq.electron_lifetime.iq.correction.fit,'electron_attenuation_us')...
                && strncmp(PulseCorrIq.electron_lifetime.iq.global.algorithm_name,'LUXkrypCal',10)       
               lifetime_values = [lifetime_values PulseCorrIq.electron_lifetime.iq.correction.fit.electron_attenuation_us];
               dataset_times = [dataset_times filename2epoch_framework(PulseCorrIq.electron_lifetime.iq.global.filename_prefix) ];
               filename_prefixs = [filename_prefixs PulseCorrIq.electron_lifetime.iq.global.filename_prefix];
               cp_number = [cp_numberPulseCorrIq.electron_lifetime.iq.global.cp_number];
        end
      end  
     end
        current_data_set_time=filename2epoch_framework(filename_prefix );

     if inrange(current_data_set_time,[min(dataset_times), max(dataset_times)] )
         [electron_lifetime] = InterpIQ(filename_prefix,dataset_times,lifetime_values);
     else    
         [index electron_lifetime_time] = NearestIQ(filename_prefix,dataset_times); 
         electron_lifetime = lifetime_values(index);
     end
    %}

    if strcmp(Iq.electron_lifetime_1.iq.global.filename_prefix,Iq.electron_lifetime_2.iq.global.filename_prefix)
        electron_lifetime = Iq.electron_lifetime_1.iq.correction.fit.electron_attenuation_us;
    else
        electron_lifetime = interp1(...
            [filename2epoch_framework(Iq.electron_lifetime_1.iq.global.filename_prefix) filename2epoch_framework(Iq.electron_lifetime_2.iq.global.filename_prefix)],...
            [Iq.electron_lifetime_1.iq.correction.fit.electron_attenuation_us Iq.electron_lifetime_2.iq.correction.fit.electron_attenuation_us],...
            filename2epoch_framework(XML_Settings.data_set),'linear');
    end
    %{  
    /// Commented Out, maybe use Z bins for S2 correction in the future -AD

    %------------------------------------------
    % Finding the S2 Z correction map values in each Z bin
    %------------------------------------------

    s2_all_z0_vals = [];
    s2_bot_z0_vals = [];
    dataset_times = [];
    filename_prefixs = {};
    cp_number = {};

     for i=1:num_iqs
      if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

        if (isfield(lug_iqs_xml.iq(i).correction.fit,'s2_bottom_z0_phe').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
    %       if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
    %           fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);
               s2_all_z0_vals = [s2_both_z0 lug_iqs_xml.iq(i).correction.fit.s2_both_z0_phe];
               s2_bot_z0_vals = [s2_both_z0 lug_iqs_xml.iq(i).correction.fit.s2_bottom_z0_phe];
               dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix ) ];
               filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
               cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
     %      else
     %          fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
     %      end

        end
      end  
     end

    current_data_set_time=filename2epoch_framework(filename_prefix );

     if inrange(current_data_set_time,[min(dataset_times), max(dataset_times)] )
         [s2_all_z0] = InterpIQ(filename_prefix,dataset_times,s2_all_z0_vals);
         [s2_bot_z0] = InterpIQ(filename_prefix,dataset_times,s2_bot_z0_vals);
     else    
         [index z0_times] = NearestIQ(filename_prefix,dataset_times); 
         s2_all_z0 = s2_all_z0_vals(index);
         s2_bot_z0 = s2_bot_z0_vals(index);
     end


    s2_z_index = [];
    dataset_times = [];
    filename_prefixs = {};
    cp_number = {};

     for i=1:num_iqs
      if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

        if (isfield(lug_iqs_xml.iq(i).correction.fit,'bin_means_s2_bottom').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
    %       if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
    %           fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);
               s2_z_index = [s2_z_index i];
               dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix ) ];
               filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
               cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
     %      else
     %          fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
     %      end

        end
      end  
     end

       [index iq_time_s2zDep] = NearestIQ(filename_prefix,dataset_times);

       s2_z_bins = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.bin_means_s2_bottom;
       s2_z_all = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.means_s2_bottom;
       s2_z_bottom = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.means_s2_both;

     %}


    %-----------------------------------
    % Finding the S1 Z-dep values
    %-----------------------------------
    %{
    z_dep_both_values = zeros(0,3);
    z_dep_bottom_values = zeros(0,3);
    dataset_times = [];
    filename_prefixs = {};
    cp_number = {};

     for i=1:num_iqs
      if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

        if (isfield(lug_iqs_xml.iq(i).correction.fit,'s1_both_zdep_quad_fit').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
    %       if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
    %           fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);
               z_dep_both_values = vertcat(z_dep_both_values,lug_iqs_xml.iq(i).correction.fit.s1_both_zdep_quad_fit);
               z_dep_bottom_values = vertcat(z_dep_bottom_values,lug_iqs_xml.iq(i).correction.fit.s1_bottom_zdep_quad_fit);

               dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix )];
               filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
               cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
    %       else
    %           fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
    %       end

        end
      end  
     end

       [index iq_time_zDep] = NearestIQ(filename_prefix,dataset_times);
       %}

    z_dep_par_all = Iq.z_dep_s1_correction.iq.correction.fit.s1_both_zdep_quad_fit;
    z_dep_par_bot = Iq.z_dep_s1_correction.iq.correction.fit.s1_bottom_zdep_quad_fit;


    %------------------------------------------
    % Finding the S2 xy correction map values
    %------------------------------------------
    %{
    s2_xy_index = [];
    dataset_times = [];
    filename_prefixs = {};
    cp_number = {};

     for i=1:num_iqs
      if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

        if (isfield(lug_iqs_xml.iq(i).correction.fit,'norm_s2_both').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
    %       if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
    %           fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);
               s2_xy_index = [s2_xy_index i];
               dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix ) ];
               filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
               cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
     %      else
     %          fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
     %      end

        end
      end  
     end

       [index iq_time_s2xyDep] = NearestIQ(filename_prefix,dataset_times);

       s2_x_bins = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.x_bin_center;
       s2_y_bins = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.y_bin_center;
       s2_map_all = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.norm_s2_both;
       s2_map_bottom = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.norm_s2_bottom;
     %}  

    s2_x_bins = Iq.s2_xy_correction.iq.correction.fit.x_bin_center;
    s2_y_bins = Iq.s2_xy_correction.iq.correction.fit.y_bin_center;
    s2_map_all = Iq.s2_xy_correction.iq.correction.fit.norm_s2_both;
    s2_map_bottom = Iq.s2_xy_correction.iq.correction.fit.norm_s2_bottom;



    %------------------------------------------
    % Finding the S1 xy correction map values
    %------------------------------------------
    %{   
    s1_xy_index = [];
    dataset_times = [];
    filename_prefixs = {};
    cp_number = {};

     for i=1:num_iqs
      if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

        if (isfield(lug_iqs_xml.iq(i).correction.fit,'norm_s1_both').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1       
    %       if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
     %          fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);
               s1_xy_index = [s1_xy_index i];                      
               dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix ) ];
               filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
               cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
     %      else
     %          fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
     %      end

        end
      end  
     end

       [index iq_time] = NearestIQ(filename_prefix,dataset_times);

       s1_x_bins = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.x_bin_center;
       s1_y_bins = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.y_bin_center;
       s1_map_all = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.norm_s1_both;
       s1_map_bottom = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.norm_s1_bottom;
    %}

    s1_x_bins = Iq.s1_xy_correction.iq.correction.fit.x_bin_center;
    s1_y_bins = Iq.s1_xy_correction.iq.correction.fit.y_bin_center;
    s1_map_all = Iq.s1_xy_correction.iq.correction.fit.norm_s1_both;
    s1_map_bottom = Iq.s1_xy_correction.iq.correction.fit.norm_s1_bottom;



    % %------------------------------------------
    % % Finding the S1 xyz correction map values
    % %------------------------------------------
    %{   
    s1_xyz_index = [];
    dataset_times = [];
    filename_prefixs = {};
    cp_number = {};

     for i=1:num_iqs
      if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

        if (isfield(lug_iqs_xml.iq(i).correction.fit,'norm_s1_both_xyz').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1       
     %      if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
     %          fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);

               s1_xyz_index = [s1_xyz_index i];                      
               dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix )];
               filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
               cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
    %       else
    %           fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
    %       end

        end
      end  
     end

       [index iq_time_s1xyzDep] = NearestIQ(filename_prefix,dataset_times);

       s1_xyz_x_bins = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.x_bin_center;
       s1_xyz_y_bins = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.y_bin_center;
       s1_xyz_z_bins = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.z_bin_center;

       xx = size(s1_xyz_x_bins); yy = size(s1_xyz_y_bins); zz = size(s1_xyz_z_bins);

       s1_xyz_map_all = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_both_xyz; s1_xyz_map_all = reshape(s1_xyz_map_all,xx(2),yy(2),zz(2));
       s1_xyz_map_bottom = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_bottom_xyz;   s1_xyz_map_bottom = reshape(s1_xyz_map_bottom,xx(2),yy(2),zz(2));

    %}
    if isfield(Iq,'s1_xyz_correction')
        s1_xyz_x_bins = Iq.s1_xyz_correction.iq.correction.fit.x_bin_center;
        s1_xyz_y_bins = Iq.s1_xyz_correction.iq.correction.fit.y_bin_center;
        s1_xyz_z_bins = Iq.s1_xyz_correction.iq.correction.fit.z_bin_center;

        xx = size(s1_xyz_x_bins); 
        yy = size(s1_xyz_y_bins); 
        zz = size(s1_xyz_z_bins);

        s1_xyz_map_all = Iq.s1_xyz_correction.iq.correction.fit.norm_s1_both_xyz;
        s1_xyz_map_all = reshape(s1_xyz_map_all,xx(2),yy(2),zz(2));
        s1_xyz_map_bottom = Iq.s1_xyz_correction.iq.correction.fit.norm_s1_bottom_xyz;   
        s1_xyz_map_bottom = reshape(s1_xyz_map_bottom,xx(2),yy(2),zz(2));
    end


    %-----------------------------------------------------------------------------------------------------------------
    %-----------------------------------------------------------------------------------------------------------------
    %% Calculating corrections ----------------------------------------------------------------------------
    %-----------------------------------------------------------------------------------------------------------------
    %----------------------------------------------------------------------------------------------------------------- 

    %% Finding the drift time, x and y associated with S1

    % Finding the drift time and xy position from the largest S2 pulse in the
    % pairing

    drift = ee.rqs.z_drift_samples;
    drift(isnan(drift)) = 0.0;
    s1_drift_ns = nan(size(ee.rqs.pulse_classification));
    s1_x_cm = nan(size(ee.rqs.pulse_classification));
    s1_y_cm = nan(size(ee.rqs.pulse_classification));
    s2_phe =  ee.rqs.pulse_area_phe.*( ismember(ee.rqs.pulse_classification,XML_Settings.s2_class_to_correct) );
    s2_phe((isnan(s2_phe))) = 0.0;

    TotalPairedS1 = ee.rqs.s1s2_pairing(find(ee.rqs.pulse_classification(:,1)==1,1,'last'));
    for ii_s1 = 1:TotalPairedS1
         [~, r] = max(s2_phe .* (ee.rqs.s1s2_pairing == ii_s1));% value and index (r for row) of Max S2
         %[c1, ~, ~] = find(ee.rqs.pulse_classification == 1,1,'first');
         %[c2 r2 v2] = find(ee.rqs.pulse_classification(:,i)==2,1,'first');
         s1_drift_ns((ee.rqs.pulse_classification == 1) & (ee.rqs.s1s2_pairing == ii_s1),1) = 10.*drift(r,1); %row c1,r and column i
         s1_x_cm((ee.rqs.pulse_classification == 1) & (ee.rqs.s1s2_pairing == ii_s1),1) = ee.rqs.x_cm(r,1);
         s1_y_cm((ee.rqs.pulse_classification == 1) & (ee.rqs.s1s2_pairing == ii_s1),1) = ee.rqs.y_cm(r,1);  
    end

    %fprintf('Done finding S1 depths :) \n');

    s1_drift_ns(s1_drift_ns==0) = nan;
    s1_x_cm(s1_x_cm==0) = nan;
    s1_y_cm(s1_y_cm==0) = nan;

    %--------------------------------------------------------------------------
    % Calculate Z corrections
    %--------------------------------------------------------------------------

    % Calculate S1 Z-correction

    s1_z_correction = polyval(z_dep_par_all,s1_ref_z_ns./1000)./polyval(z_dep_par_all,s1_drift_ns./1000);
    s1_z_correction(isnan(s1_z_correction)) = 1.0;

    s1_z_correction_bot = polyval(z_dep_par_bot,s1_ref_z_ns./1000)./polyval(z_dep_par_bot,s1_drift_ns./1000);
    s1_z_correction_bot(isnan(s1_z_correction_bot)) = 1.0;

    % Calculate electron lifetime correction (S2 Z-correction)
    % Reading the values of electron lifetime from the LUG and interpolating
    % between the two either side

    electron_lifetime_correction = exp(double(ee.rqs.z_drift_samples)./(100.*electron_lifetime));
    electron_lifetime_correction(isnan(electron_lifetime_correction) | ~ismember(ee.rqs.pulse_classification,XML_Settings.s2_class_to_correct)) = 1.0;

    %--------------------------------------------------------------------------
    % Calculate XYZ corrections
    %--------------------------------------------------------------------------

    % Calculate S1 xyz corrections from map stored in IQ

    %s1_x_cm(isnan(s1_x_cm)) = 0; s1_y_cm(isnan(s1_y_cm)) = 0; s1_drift_ns(isnan(s1_drift_ns)) = 0;

    s1xyz_correction = interp3(s1_xyz_x_bins,s1_xyz_y_bins,s1_xyz_z_bins,s1_xyz_map_all,s1_x_cm,s1_y_cm,s1_drift_ns./1000,'spline');
    s1xyz_correction((isnan(s1xyz_correction))) = 1.0;

    s1xyz_correction_bot = interp3(s1_xyz_x_bins,s1_xyz_y_bins,s1_xyz_z_bins,s1_xyz_map_bottom,s1_x_cm,s1_y_cm,s1_drift_ns./1000,'spline');
    s1xyz_correction_bot((isnan(s1xyz_correction_bot))) = 1.0;


    % Calculate S1 XY-only corrections

    % s2x = dp.x_cm.*(+(dp.pulse_classification==2)); s2x(find(s2x==0)) = nan;
    % s2y = dp.y_cm.*(+(dp.pulse_classification==2)); s2y(find(s2y==0)) = nan;
    % 
    % s2xy_correction = interp2(s1_x_bins,s1_y_bins,s2_map_all,dp.x_cm,dp.y_cm,'spline');
    % s2xy_correction = s2xy_correction.*(+dp.pulse_classification==2);
    % s2xy_correction(find(s2xy_correction==0))=1.0;  s2xy_correction(find(isnan(s2xy_correction)))=1.0;
    % 
    % s2xy_correction_bot = interp2(s1_x_bins,s1_y_bins,s2_map_bottom,dp.x_cm,dp.y_cm,'spline');
    % s2xy_correction_bot = s2xy_correction_bot.*(+dp.pulse_classification==2);
    % s2xy_correction_bot(find(s2xy_correction_bot==0))=1.0;  s2xy_correction_bot(find(isnan(s2xy_correction_bot)))=1.0;




    % Calculate S2 XY corrections

    s2x = ee.rqs.x_cm.*(+(ee.rqs.pulse_classification==2)); s2x((s2x==0)) = nan;
    s2y = ee.rqs.y_cm.*(+(ee.rqs.pulse_classification==2)); s2y((s2y==0)) = nan;

    s2xy_correction = interp2(s2_x_bins,s2_y_bins,s2_map_all,ee.rqs.x_cm,ee.rqs.y_cm,'spline');
    s2xy_correction = s2xy_correction.*(+ee.rqs.pulse_classification==2);
    s2xy_correction((s2xy_correction==0))=1.0;  s2xy_correction((isnan(s2xy_correction)))=1.0;

    s2xy_correction_bot = interp2(s2_x_bins,s2_y_bins,s2_map_bottom,ee.rqs.x_cm,ee.rqs.y_cm,'spline');
    s2xy_correction_bot = s2xy_correction_bot.*(+ee.rqs.pulse_classification==2);
    s2xy_correction_bot((s2xy_correction_bot==0))=1.0;  s2xy_correction_bot((isnan(s2xy_correction_bot)))=1.0;


    %% add RQs of correction factors

    ee.rqs.correction_electron_lifetime     = electron_lifetime_correction;
    ee.rqs.correction_s1_z_dependence       = s1_z_correction;
    ee.rqs.correction_s1_z_dependence_bot   = s1_z_correction_bot;
    ee.rqs.correction_s1_xyz_dependence     = s1xyz_correction; %this is either 3D or 1D+2D
    ee.rqs.correction_s1_xyz_dependence_bot = s1xyz_correction_bot; %this is either 3D or 1D+2D
    ee.rqs.correction_s2_xy_dependence      = s2xy_correction;
    ee.rqs.correction_s2_xy_dependence_bot  = s2xy_correction_bot;

    ee.rqs.admin.corrections.electron_lifetime        = electron_lifetime;
    %ee.rqs.admin.corrections.s1_z_dependence.iq_time  = iq_time_zDep;
    ee.rqs.admin.corrections.s1_z_dependence.all      = z_dep_par_all;
    ee.rqs.admin.corrections.s1_z_dependence.bottom   = z_dep_par_bot;
    %ee.rqs.admin.corrections.s1_xy_dependence.iq_time = iq_time;
    %ee.rqs.admin.corrections.s2_xy_dependence.iq_time = iq_time_s2xyDep;


    %-----------------------------------------------------------------------------------------------------------------
    %-----------------------------------------------------------------------------------------------------------------
    %% Applying corrections ----------------------------------------------------------------------------
    %-----------------------------------------------------------------------------------------------------------------
    %----------------------------------------------------------------------------------------------------------------- 

    %--------------------------------------------------------------------------
    % Apply Z corrections
    %--------------------------------------------------------------------------

    ee.rqs.z_corrected_pulse_area_all_phe = ee.rqs.z_corrected_pulse_area_all_phe.*electron_lifetime_correction.*s1_z_correction; 
    ee.rqs.z_corrected_pulse_area_bot_phe = ee.rqs.z_corrected_pulse_area_bot_phe.*electron_lifetime_correction.*s1_z_correction_bot; 

    %--------------------------------------------------------------------------
    % Apply XYZ corrections
    %--------------------------------------------------------------------------

    ee.rqs.xyz_corrected_pulse_area_all_phe = ee.rqs.xyz_corrected_pulse_area_all_phe.*electron_lifetime_correction.*s2xy_correction.*s1xyz_correction;
    ee.rqs.xyz_corrected_pulse_area_bot_phe = ee.rqs.xyz_corrected_pulse_area_bot_phe.*electron_lifetime_correction.*s2xy_correction_bot.*s1xyz_correction_bot;

    ee.info.PulseCorrSuccess = 1;
 	ee.info.PulseCorrError = '';
 	ee.info.PulseCorrVersion = SR_Version;
elseif strcmp(SR_Version,'SR2.1')
    %-----------------------------------------------------------------------------------------------------------------
    %-----------------------------------------------------------------------------------------------------------------
    %% Reading iq values from the IQs xml ----------------------------------------------------------------------------
    %-----------------------------------------------------------------------------------------------------------------
    %----------------------------------------------------------------------------------------------------------------- 
    % electron_lifetime - Interpolate - DONE
    % z_dep_par_all - Nearest - DONE
    % z_dep_par_bot - Nearest - DONE
    % x_bins,y_bins,z_bins,s1_map_all - Nearest
    % x_bins,y_bins,z_bins,s1_map_bot - Nearest
    % x_bins,y_bins,s2_map_all - Nearest
    % x_bins,y_bins,s2_map_bot - Nearest

    %[a num_iqs] = size(lug_iqs_xml.iq);

    %-----------------------------------
    % Finding the electron lifetime
    %-----------------------------------

    lifetime_values = [];
    dataset_times = [];
    filename_prefixs = {};
    cp_number = {};
    %Why do this ???
    %{
     for i=1:num_iqs
      if isfield(PulseCorrIq.electron_lifetime.iq.correction,'fit')
        if isfield(PulseCorrIq.electron_lifetime.iq.correction.fit,'electron_attenuation_us')...
                && strncmp(PulseCorrIq.electron_lifetime.iq.global.algorithm_name,'LUXkrypCal',10)       
               lifetime_values = [lifetime_values PulseCorrIq.electron_lifetime.iq.correction.fit.electron_attenuation_us];
               dataset_times = [dataset_times filename2epoch_framework(PulseCorrIq.electron_lifetime.iq.global.filename_prefix) ];
               filename_prefixs = [filename_prefixs PulseCorrIq.electron_lifetime.iq.global.filename_prefix];
               cp_number = [cp_numberPulseCorrIq.electron_lifetime.iq.global.cp_number];
        end
      end  
     end
        current_data_set_time=filename2epoch_framework(filename_prefix );

     if inrange(current_data_set_time,[min(dataset_times), max(dataset_times)] )
         [electron_lifetime] = InterpIQ(filename_prefix,dataset_times,lifetime_values);
     else    
         [index electron_lifetime_time] = NearestIQ(filename_prefix,dataset_times); 
         electron_lifetime = lifetime_values(index);
     end
    %}
    if strcmp(Iq.electron_lifetime_1.iq.global.filename_prefix,Iq.electron_lifetime_2.iq.global.filename_prefix)
        electron_lifetime = Iq.electron_lifetime_1.iq.correction.fit.electron_attenuation_us;
    else
        electron_lifetime = interp1(...
            [filename2epoch_framework(Iq.electron_lifetime_1.iq.global.filename_prefix) filename2epoch_framework(Iq.electron_lifetime_2.iq.global.filename_prefix)],...
            [Iq.electron_lifetime_1.iq.correction.fit.electron_attenuation_us Iq.electron_lifetime_2.iq.correction.fit.electron_attenuation_us],...
            filename2epoch_framework(XML_Settings.data_set),'linear');
    end
    %{  
    /// Commented Out, maybe use Z bins for S2 correction in the future -AD

    %------------------------------------------
    % Finding the S2 Z correction map values in each Z bin
    %------------------------------------------

    s2_all_z0_vals = [];
    s2_bot_z0_vals = [];
    dataset_times = [];
    filename_prefixs = {};
    cp_number = {};

     for i=1:num_iqs
      if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

        if (isfield(lug_iqs_xml.iq(i).correction.fit,'s2_bottom_z0_phe').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
    %       if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
    %           fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);
               s2_all_z0_vals = [s2_both_z0 lug_iqs_xml.iq(i).correction.fit.s2_both_z0_phe];
               s2_bot_z0_vals = [s2_both_z0 lug_iqs_xml.iq(i).correction.fit.s2_bottom_z0_phe];
               dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix ) ];
               filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
               cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
     %      else
     %          fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
     %      end

        end
      end  
     end

    current_data_set_time=filename2epoch_framework(filename_prefix );

     if inrange(current_data_set_time,[min(dataset_times), max(dataset_times)] )
         [s2_all_z0] = InterpIQ(filename_prefix,dataset_times,s2_all_z0_vals);
         [s2_bot_z0] = InterpIQ(filename_prefix,dataset_times,s2_bot_z0_vals);
     else    
         [index z0_times] = NearestIQ(filename_prefix,dataset_times); 
         s2_all_z0 = s2_all_z0_vals(index);
         s2_bot_z0 = s2_bot_z0_vals(index);
     end


    s2_z_index = [];
    dataset_times = [];
    filename_prefixs = {};
    cp_number = {};

     for i=1:num_iqs
      if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

        if (isfield(lug_iqs_xml.iq(i).correction.fit,'bin_means_s2_bottom').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
    %       if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
    %           fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);
               s2_z_index = [s2_z_index i];
               dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix ) ];
               filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
               cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
     %      else
     %          fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
     %      end

        end
      end  
     end

       [index iq_time_s2zDep] = NearestIQ(filename_prefix,dataset_times);

       s2_z_bins = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.bin_means_s2_bottom;
       s2_z_all = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.means_s2_bottom;
       s2_z_bottom = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.means_s2_both;

     %}


    %-----------------------------------
    % Finding the S1 Z-dep values
    %-----------------------------------
    %{
    z_dep_both_values = zeros(0,3);
    z_dep_bottom_values = zeros(0,3);
    dataset_times = [];
    filename_prefixs = {};
    cp_number = {};

     for i=1:num_iqs
      if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

        if (isfield(lug_iqs_xml.iq(i).correction.fit,'s1_both_zdep_quad_fit').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
    %       if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
    %           fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);
               z_dep_both_values = vertcat(z_dep_both_values,lug_iqs_xml.iq(i).correction.fit.s1_both_zdep_quad_fit);
               z_dep_bottom_values = vertcat(z_dep_bottom_values,lug_iqs_xml.iq(i).correction.fit.s1_bottom_zdep_quad_fit);

               dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix )];
               filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
               cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
    %       else
    %           fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
    %       end

        end
      end  
     end

       [index iq_time_zDep] = NearestIQ(filename_prefix,dataset_times);
       %}

    z_dep_par_all = Iq.z_dep_s1_correction.iq.correction.fit.s1_both_zdep_quad_fit;
    z_dep_par_bot = Iq.z_dep_s1_correction.iq.correction.fit.s1_bottom_zdep_quad_fit;


    %------------------------------------------
    % Finding the S2 xy correction map values
    %------------------------------------------
    %{
    s2_xy_index = [];
    dataset_times = [];
    filename_prefixs = {};
    cp_number = {};

     for i=1:num_iqs
      if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

        if (isfield(lug_iqs_xml.iq(i).correction.fit,'norm_s2_both').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
    %       if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
    %           fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);
               s2_xy_index = [s2_xy_index i];
               dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix ) ];
               filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
               cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
     %      else
     %          fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
     %      end

        end
      end  
     end

       [index iq_time_s2xyDep] = NearestIQ(filename_prefix,dataset_times);

       s2_x_bins = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.x_bin_center;
       s2_y_bins = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.y_bin_center;
       s2_map_all = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.norm_s2_both;
       s2_map_bottom = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.norm_s2_bottom;
     %}  

    s2_x_bins = Iq.s2_xy_correction.iq.correction.fit.x_bin_center;
    s2_y_bins = Iq.s2_xy_correction.iq.correction.fit.y_bin_center;
    s2_map_all = Iq.s2_xy_correction.iq.correction.fit.norm_s2_both;
    s2_map_bottom = Iq.s2_xy_correction.iq.correction.fit.norm_s2_bottom;



    %------------------------------------------
    % Finding the S1 xy correction map values
    %------------------------------------------
    %{   
    s1_xy_index = [];
    dataset_times = [];
    filename_prefixs = {};
    cp_number = {};

     for i=1:num_iqs
      if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

        if (isfield(lug_iqs_xml.iq(i).correction.fit,'norm_s1_both').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1       
    %       if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
     %          fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);
               s1_xy_index = [s1_xy_index i];                      
               dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix ) ];
               filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
               cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
     %      else
     %          fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
     %      end

        end
      end  
     end

       [index iq_time] = NearestIQ(filename_prefix,dataset_times);

       s1_x_bins = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.x_bin_center;
       s1_y_bins = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.y_bin_center;
       s1_map_all = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.norm_s1_both;
       s1_map_bottom = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.norm_s1_bottom;
    %}

    s1_x_bins = Iq.s1_xy_correction.iq.correction.fit.x_bin_center;
    s1_y_bins = Iq.s1_xy_correction.iq.correction.fit.y_bin_center;
    s1_map_all = Iq.s1_xy_correction.iq.correction.fit.norm_s1_both;
    s1_map_bottom = Iq.s1_xy_correction.iq.correction.fit.norm_s1_bottom;



    % %------------------------------------------
    % % Finding the S1 xyz correction map values
    % %------------------------------------------
    %{   
    s1_xyz_index = [];
    dataset_times = [];
    filename_prefixs = {};
    cp_number = {};

     for i=1:num_iqs
      if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

        if (isfield(lug_iqs_xml.iq(i).correction.fit,'norm_s1_both_xyz').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1       
     %      if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
     %          fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);

               s1_xyz_index = [s1_xyz_index i];                      
               dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix )];
               filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
               cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
    %       else
    %           fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
    %       end

        end
      end  
     end

       [index iq_time_s1xyzDep] = NearestIQ(filename_prefix,dataset_times);

       s1_xyz_x_bins = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.x_bin_center;
       s1_xyz_y_bins = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.y_bin_center;
       s1_xyz_z_bins = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.z_bin_center;

       xx = size(s1_xyz_x_bins); yy = size(s1_xyz_y_bins); zz = size(s1_xyz_z_bins);

       s1_xyz_map_all = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_both_xyz; s1_xyz_map_all = reshape(s1_xyz_map_all,xx(2),yy(2),zz(2));
       s1_xyz_map_bottom = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_bottom_xyz;   s1_xyz_map_bottom = reshape(s1_xyz_map_bottom,xx(2),yy(2),zz(2));

    %}
    if isfield(Iq,'s1_xyz_correction')
        s1_xyz_x_bins = Iq.s1_xyz_correction.iq.correction.fit.x_bin_center;
        s1_xyz_y_bins = Iq.s1_xyz_correction.iq.correction.fit.y_bin_center;
        s1_xyz_z_bins = Iq.s1_xyz_correction.iq.correction.fit.z_bin_center;

        xx = size(s1_xyz_x_bins); 
        yy = size(s1_xyz_y_bins); 
        zz = size(s1_xyz_z_bins);

        s1_xyz_map_all = Iq.s1_xyz_correction.iq.correction.fit.norm_s1_both_xyz;
        s1_xyz_map_all = reshape(s1_xyz_map_all,xx(2),yy(2),zz(2));
        s1_xyz_map_bottom = Iq.s1_xyz_correction.iq.correction.fit.norm_s1_bottom_xyz;   
        s1_xyz_map_bottom = reshape(s1_xyz_map_bottom,xx(2),yy(2),zz(2));
    end


    %-----------------------------------------------------------------------------------------------------------------
    %-----------------------------------------------------------------------------------------------------------------
    %% Calculating corrections ----------------------------------------------------------------------------
    %-----------------------------------------------------------------------------------------------------------------
    %----------------------------------------------------------------------------------------------------------------- 

    %% Finding the drift time, x and y associated with S1

    % Finding the drift time and xy position from the largest S2 pulse in the
    % pairing

    drift = ee.rqs.z_drift_samples;
    drift(isnan(drift)) = 0.0;
    s1_drift_ns = nan(size(ee.rqs.pulse_classification));
    s1_x_cm = nan(size(ee.rqs.pulse_classification));
    s1_y_cm = nan(size(ee.rqs.pulse_classification));
    s2_phe =  ee.rqs.pulse_area_phe.*( ismember(ee.rqs.pulse_classification,XML_Settings.s2_class_to_correct) );
    s2_phe((isnan(s2_phe))) = 0.0;

    TotalPairedS1 = ee.rqs.s1s2_pairing(find(ee.rqs.pulse_classification(:,1)==1,1,'last'));
    for ii_s1 = 1:TotalPairedS1
         [~, r] = max(s2_phe .* (ee.rqs.s1s2_pairing == ii_s1));% value and index (r for row) of Max S2
         %[c1, ~, ~] = find(ee.rqs.pulse_classification == 1,1,'first');
         %[c2 r2 v2] = find(ee.rqs.pulse_classification(:,i)==2,1,'first');
         s1_drift_ns((ee.rqs.pulse_classification == 1) & (ee.rqs.s1s2_pairing == ii_s1),1) = 10.*drift(r,1); %row c1,r and column i
         s1_x_cm((ee.rqs.pulse_classification == 1) & (ee.rqs.s1s2_pairing == ii_s1),1) = ee.rqs.x_cm(r,1);
         s1_y_cm((ee.rqs.pulse_classification == 1) & (ee.rqs.s1s2_pairing == ii_s1),1) = ee.rqs.y_cm(r,1);  
    end

    %fprintf('Done finding S1 depths :) \n');

    s1_drift_ns(s1_drift_ns==0) = nan;
    s1_x_cm(s1_x_cm==0) = nan;
    s1_y_cm(s1_y_cm==0) = nan;

    %--------------------------------------------------------------------------
    % Calculate Z corrections
    %--------------------------------------------------------------------------

    % Calculate S1 Z-correction

    s1_z_correction = polyval(z_dep_par_all,s1_ref_z_ns./1000)./polyval(z_dep_par_all,s1_drift_ns./1000);
    s1_z_correction(isnan(s1_z_correction)) = 1.0;

    s1_z_correction_bot = polyval(z_dep_par_bot,s1_ref_z_ns./1000)./polyval(z_dep_par_bot,s1_drift_ns./1000);
    s1_z_correction_bot(isnan(s1_z_correction_bot)) = 1.0;

    % Calculate electron lifetime correction (S2 Z-correction)
    % Reading the values of electron lifetime from the LUG and interpolating
    % between the two either side

    electron_lifetime_correction = exp(double(ee.rqs.z_drift_samples)./(100.*electron_lifetime));
    electron_lifetime_correction(isnan(electron_lifetime_correction) | ~ismember(ee.rqs.pulse_classification,XML_Settings.s2_class_to_correct)) = 1.0;

    %--------------------------------------------------------------------------
    % Calculate XY corrections
    %--------------------------------------------------------------------------


    % Calculate S1 XY-only corrections

    %s2x = ee.rqs.x_cm.*((ee.rqs.pulse_classification==2) & (ee.rqs.pulse_classification==4)); 
    %s2x(s2x==0) = nan;
    %s2y = ee.rqs.y_cm.*((ee.rqs.pulse_classification==2) & (ee.rqs.pulse_classification==4)); 
    %s2y(s2y==0) = nan;

    s1xy_correction = interp2(s1_x_bins,s1_y_bins,s1_map_all,s1_x_cm,s1_y_cm,'spline',1);
    s1xy_correction = s1xy_correction.*(ee.rqs.pulse_classification == 1);
    s1xy_correction(s1xy_correction == 0)   = 1.0;  
    s1xy_correction(isnan(s1xy_correction)) = 1.0;

    s1xy_correction_bot = interp2(s1_x_bins,s1_y_bins,s1_map_bottom,s1_x_cm,s1_y_cm,'spline',1);
    s1xy_correction_bot = s1xy_correction_bot.*(ee.rqs.pulse_classification == 1);
    s1xy_correction_bot(s1xy_correction_bot == 0)   = 1.0;  
    s1xy_correction_bot(isnan(s1xy_correction_bot)) = 1.0;

    % Calculate S2 XY corrections

    s2xy_correction = interp2(s2_x_bins,s2_y_bins,s2_map_all,ee.rqs.x_cm,ee.rqs.y_cm,'spline',1);
    s2xy_correction = s2xy_correction.*( ismember(ee.rqs.pulse_classification,XML_Settings.s2_class_to_correct) );
    s2xy_correction(s2xy_correction == 0)   = 1.0;  
    s2xy_correction(isnan(s2xy_correction)) = 1.0;

    s2xy_correction_bot = interp2(s2_x_bins,s2_y_bins,s2_map_bottom,ee.rqs.x_cm,ee.rqs.y_cm,'spline',1);
    s2xy_correction_bot = s2xy_correction_bot.*( ismember(ee.rqs.pulse_classification,XML_Settings.s2_class_to_correct) ) ;
    s2xy_correction_bot(s2xy_correction_bot == 0)   = 1.0;  
    s2xy_correction_bot(isnan(s2xy_correction_bot)) = 1.0;

    %--------------------------------------------------------------------------
    % Calculate XYZ corrections
    %--------------------------------------------------------------------------


    % Use S1 XYZ correction if available, otherwise use S1 XY + Z correction
    if isfield(Iq,'s1_xyz_correction')    
        s1xyz_correction = interp3(s1_xyz_x_bins,s1_xyz_y_bins,s1_xyz_z_bins,s1_xyz_map_all,s1_x_cm,s1_y_cm,s1_drift_ns./1000,'spline');
        s1xyz_correction(isnan(s1xyz_correction)) = 1.0;

        s1xyz_correction_bot = interp3(s1_xyz_x_bins,s1_xyz_y_bins,s1_xyz_z_bins,s1_xyz_map_bottom,s1_x_cm,s1_y_cm,s1_drift_ns./1000,'spline');
        s1xyz_correction_bot(isnan(s1xyz_correction_bot)) = 1.0;
        ee.rqs.s1_correction_is_3d(:) = 1;
    else
        s1xyz_correction=s1xy_correction.*s1_z_correction;
        s1xyz_correction_bot=s1xy_correction_bot.*s1_z_correction_bot;
        ee.rqs.s1_correction_is_3d(:) = 0;
    end

    %% add RQs of correction factors

    ee.rqs.correction_electron_lifetime     = electron_lifetime_correction;
    ee.rqs.correction_s1_z_dependence       = s1_z_correction;
    ee.rqs.correction_s1_z_dependence_bot   = s1_z_correction_bot;
    ee.rqs.correction_s1_xyz_dependence     = s1xyz_correction; %this is either 3D or 1D+2D
    ee.rqs.correction_s1_xyz_dependence_bot = s1xyz_correction_bot; %this is either 3D or 1D+2D
    ee.rqs.correction_s1_xy_dependence      = s1xy_correction;
    ee.rqs.correction_s1_xy_dependence_bot  = s1xy_correction_bot;
    ee.rqs.correction_s2_xy_dependence      = s2xy_correction;
    ee.rqs.correction_s2_xy_dependence_bot  = s2xy_correction_bot;

    ee.rqs.admin.corrections.electron_lifetime        = electron_lifetime;
    %ee.rqs.admin.corrections.s1_z_dependence.iq_time  = iq_time_zDep;
    ee.rqs.admin.corrections.s1_z_dependence.all      = z_dep_par_all;
    ee.rqs.admin.corrections.s1_z_dependence.bottom   = z_dep_par_bot;
    %ee.rqs.admin.corrections.s1_xy_dependence.iq_time = iq_time;
    %ee.rqs.admin.corrections.s2_xy_dependence.iq_time = iq_time_s2xyDep;


    %-----------------------------------------------------------------------------------------------------------------
    %-----------------------------------------------------------------------------------------------------------------
    %% Applying corrections ----------------------------------------------------------------------------
    %-----------------------------------------------------------------------------------------------------------------
    %----------------------------------------------------------------------------------------------------------------- 

    %--------------------------------------------------------------------------
    % Apply Z corrections
    %--------------------------------------------------------------------------

    ee.rqs.z_corrected_pulse_area_all_phe = ee.rqs.z_corrected_pulse_area_all_phe.*electron_lifetime_correction.*s1_z_correction; 
    ee.rqs.z_corrected_pulse_area_bot_phe = ee.rqs.z_corrected_pulse_area_bot_phe.*electron_lifetime_correction.*s1_z_correction_bot; 

    %--------------------------------------------------------------------------
    % Apply XYZ corrections
    %--------------------------------------------------------------------------

    ee.rqs.xyz_corrected_pulse_area_all_phe = ee.rqs.xyz_corrected_pulse_area_all_phe.*electron_lifetime_correction.*s2xy_correction.*s1xyz_correction;
    ee.rqs.xyz_corrected_pulse_area_bot_phe = ee.rqs.xyz_corrected_pulse_area_bot_phe.*electron_lifetime_correction.*s2xy_correction_bot.*s1xyz_correction_bot;

    ee.info.PulseCorrSuccess = 1;
    ee.info.PulseCorrError = '';
 	ee.info.PulseCorrVersion = SR_Version;
else
	error('Do not recognize SR version.');
end   
catch exception
    if Debug
        disp('Fail to run pulse area correction.')
    end
    
    if strcmp(Iq.electron_lifetime_1.iq.global.filename_prefix,Iq.electron_lifetime_2.iq.global.filename_prefix)
        electron_lifetime = Iq.electron_lifetime_1.iq.correction.fit.electron_attenuation_us;
    else
        electron_lifetime = interp1(...
            [filename2epoch_framework(Iq.electron_lifetime_1.iq.global.filename_prefix) filename2epoch_framework(Iq.electron_lifetime_2.iq.global.filename_prefix)],...
            [Iq.electron_lifetime_1.iq.correction.fit.electron_attenuation_us Iq.electron_lifetime_2.iq.correction.fit.electron_attenuation_us],...
            filename2epoch_framework(XML_Settings.data_set),'linear');
    end
    z_dep_par_all = Iq.z_dep_s1_correction.iq.correction.fit.s1_both_zdep_quad_fit;
    z_dep_par_bot = Iq.z_dep_s1_correction.iq.correction.fit.s1_bottom_zdep_quad_fit;
    
    ee.rqs.admin.corrections.electron_lifetime        = electron_lifetime;
    ee.rqs.admin.corrections.s1_z_dependence.all      = z_dep_par_all;
    ee.rqs.admin.corrections.s1_z_dependence.bottom   = z_dep_par_bot;
    
    ee.info.PulseCorrSuccess = 0;
    ee.info.PulseCorrError = exception.identifier;
    ee.info.PulseCorrVersion = SR_Version;
end
