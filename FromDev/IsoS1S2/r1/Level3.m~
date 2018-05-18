function [ee] = DatPulseFinderClassifierExtended_BH(raw, TrgTime, EventWindow, Setting, DataPath, DataSet)
%For debug
%
%DataPath = 'G:\LUXData\';
%DataSet = 'lux10_20130728T2255';
%DataSet = 'lux10_20150108T0855';
%DataSet = 'lux10_20150930T0132';
%{
Debug.load_evt = 0;
Debug.load_rqs = 0;
Setting.debug.verbose = 0;

if exist([DataPath '\dat\' DataSet ],'dir')
    DatPath = [DataPath '\dat\' DataSet filesep];
    DatFileList = dir([DatPath '*.dat']);
else
    disp('No dat files found')
    return
end
%%
if Debug.load_evt
    disp('Loading evt file')
    if exist([DataPath '\evt\' DataSet ],'dir')
        EvtPath = [DataPath '\evt\' DataSet filesep];
        EvtFileList = dir([EvtPath '*.evt']);
        [event_struct setting] = LUXEventLoader_framework(EvtPath,EvtFileList(1).name);
        Setting.wrapper.xml_file_path = 'D:\Matlab\Code\Scratch\DatPFC_Extended\SettingXMLFiles\';
        PMTGain = XMLReader_framework([Setting.wrapper.xml_file_path 'PMTVUVGainV1_3.xml']);
        amp_gain = 5*1.5;
        pmt_gains_mVns_per_phe = [PMTGain.iq(1).fit.channel.mVns_per_phe];
        event_struct = DatPFC_LUXBaselineZen_framework(event_struct);
        event_struct(1).thr = 0;
        event_struct = DatPFC_LUXCalibratePulses_framework(event_struct,pmt_gains_mVns_per_phe,amp_gain);
        event_struct = DatPFC_LUXSumPOD_framework(event_struct);
        for evt = 1:length(event_struct)
            event_struct(evt).sumpod_time_samples = permute(event_struct(evt).sumpod_time_samples,[2 1]);
        end
    else
        disp('No evt files found')
        return
    end
end
%%
if Debug.load_rqs
    disp('Loading rqs file')
    RqsFolderList = dir([DataPath '\rqs\']);
    Temp = struct2cell(RqsFolderList);
    ID = find(cellfun(@isempty,strfind(Temp(1,:),DataSet)) == 0);
    if ~isempty(ID)
        RqsPath = [DataPath '\rqs\' RqsFolderList(ID).name filesep];
        [d] = LUXLoadMultipleRQMs_framework(RqsPath,'all',10);
    else
        disp('No rqs files found')
        return
    end
end
%% First load all setting files
%ff = 3;
%Setting.wrapper.max_num_pulses = 20000;%Need to find a better place to put this line
%Setting.wrapper.xml_file_path = '/home/mmoongwe/Code/DatPFC/';
if Setting.debug.verbose
    disp('Loading dat file.')
end
[raw, livetime, ~] = LUXLoadRawDataFile(DatFileList(ff).name,DatPath);

clear ee
TrgTime = livetime.latch;
%TrgTime = raw.ch(128).pod(1).timestamp;
%EventWindow = [0 1e7];
EventWindow = [0 (livetime.end-livetime.latch+1)];
%EventWindow = [0 (livetime.end-livetime.latch+1)];
%EventWindow = [-1e5 +1e5];

Setting.wrapper.xml_file_path = 'D:\Matlab\Code\Scratch\DatPFC_Extended\SettingXMLFiles\';
%}


if Setting.debug.verbose
    disp('Constructing waveform from dat file.')
end

ee = DatPFC_WaveformMaker(raw,TrgTime,EventWindow,Setting.dp_module_setting.sr_version);

if all([ee.ch.empty] == 1)
    ee.empty = 1;
    ee.info.EventConstructionSuccess = 0;
    ee.info.EventConstructionError = 'No pods found in the file.';
end

gs_xml = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.gs_file_name]);
for ii_gs = 1:length(gs_xml.data_processing_settings.module)
   ModuleSetting.(gs_xml.data_processing_settings.module(ii_gs).module_name).parameters = ...
       gs_xml.data_processing_settings.module(ii_gs).parameters;
end

%% BaselineZen
if Setting.dp_module_setting.LUXBaselineZen.run_module
    if Setting.debug.verbose
        disp('Adjusting baseline.')
    end
    ee = DatPFC_LUXBaselineZen_framework(ee,Setting.debug.verbose);
end
%% Calibrate pulses
if Setting.dp_module_setting.pulse_calibration.run_module
    if Setting.debug.verbose
        disp('Calibrating pulses.')
    end
    PMTGain = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.pulse_calibration.pmt_gain_xml]);
    amp_gain = 5*1.5;
    pmt_gains_mVns_per_phe = [PMTGain.iq(1).fit.channel.mVns_per_phe];
    ee.thr = 0;
    ee = DatPFC_LUXCalibratePulses_framework(ee,pmt_gains_mVns_per_phe,amp_gain,Setting.debug.verbose);
end
%% LUXSumPOD
if Setting.dp_module_setting.LUXSumPOD.run_module
    if Setting.debug.verbose
        disp('Sum waveform.')
    end
    ee = DatPFC_LUXSumPOD_framework(ee,1:122,Setting.dp_module_setting.sr_version,Setting.debug.verbose);
    if ee.empty == 1
        ee.sumpod_time_samples = [];
    end
    ee.sumpod_time_samples = ee.sumpod_time_samples';
end
%% PulseFinder_TRC.
if Setting.dp_module_setting.pulse_finder.run_module
    if Setting.debug.verbose
        disp('Running pulse finder.')
    end

    %PulseFinder_XML = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.pulse_finder.setting_xml]);
    PulseFinder_Settings = ModuleSetting.PulseFinder_TransparentRubiksCube.parameters;
    PulseFinder_Settings.EventWindow = ee.info.EventWindow;
    PulseFinder_Settings.max_num_pulses = Setting.wrapper.max_num_pulses;

    ee = DatPFC_PulseFinder_TRC(ee,PulseFinder_Settings,Setting.dp_module_setting.sr_version,Setting.debug.verbose);%Make sure to double check the xml and perusepeek.c version

    if Setting.debug.verbose
        disp([num2str(ee.rqs.num_pulses_found) ' pulses found.'])
    end

    %%So far up to this part seems to work almost okay.
    %%Two standing problems
    %%1.The start/end time is not 100% match. Still a few samples of from time
    %%to time.
    %%2.Some small pulses are missed. -> Need benchmark.
    if ee.rqs.num_pulses_found == 0
        ee.info.PulseFinderSuccess = 0;
        ee.info.PulseFinderError = 'No pulse found in this dat file.';
    end
end


%% PulseTiming
if Setting.dp_module_setting.pulse_timing.run_module
    if Setting.debug.verbose
        disp('Calculate pulse timings.')
    end

    %PulseTiming_XML = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.pulse_timing.setting_xml]);
    PulseTiming_Settings = ModuleSetting.PulseTiming_HeightTiming.parameters;
    PulseTiming_Settings.max_num_pulses = ee.rqs.num_pulses_found;

    ee = DatPFC_PulseTiming_HeightTiming_Converted(ee,PulseTiming_Settings,Setting.debug.verbose);
end
%% PulseQuantities
if Setting.dp_module_setting.pulse_quantities.run_module
    if Setting.debug.verbose
        disp('Calculate pulse basic quantities.')
    end
    %PulseQuantities_XML = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.pulse_quantity.setting_xml]);
    PulseQuantities_Settings = ModuleSetting.PulseQuantities_MinimumSet.parameters;
    PulseQuantities_Settings.max_num_pulses = ee.rqs.num_pulses_found;

    ee = DatPFC_PulseQuantities_MinimumSet(ee,PulseQuantities_Settings,Setting.debug.verbose);
end
%% Next real pulse classifier. Excerp from
if Setting.dp_module_setting.pulse_classifier.run_module
    if Setting.debug.verbose
        disp('Running pulse classifier.')
    end
    PulseClassifier_Settings.max_num_pulses = ee.rqs.num_pulses_found;
    ee = DatPFC_PulseClassifier_MultiDimensional(ee,PulseClassifier_Settings,Setting.dp_module_setting.sr_version,Setting.debug.verbose);
end

%% Photon counting
if Setting.dp_module_setting.photon_counting.run_module
    if Setting.debug.verbose
        disp('Running photon counting.')
    end

    %PhotonCounting_XML = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.photon_counting.setting_xml]);
    PhotonCounting_Settings = ModuleSetting.PulseQuantities_PhotonCounting.parameters;
    PhotonCounting_Settings.max_num_pulses = ee.rqs.num_pulses_found;

    ee = DatPFC_PhotonCounting(ee,PhotonCounting_Settings,Setting.debug.verbose);
end
%% S1S2 pairing
if Setting.dp_module_setting.s1s2_pairing.run_module
    if Setting.debug.verbose
        disp('Pairing S1/S2.')
    end

    S1S2Finder_Settings.max_num_pulses = ee.rqs.num_pulses_found;

    ee = DatPFC_S1S2Finder(ee,S1S2Finder_Settings,Setting.debug.verbose);
    %z_drift is simple aft0_s2-aft0s1
end
%% PosRecon
if Setting.dp_module_setting.position_reconstruction.run_module
    if Setting.debug.verbose
        disp('Running position reconstruction.')
    end

    %PosRecon_XML = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.position_reconstruction.setting_xml]);
    PosRecon_Settings = ModuleSetting.PositionReconstruction_MercuryI.parameters;
    PosReconIq.lrfs = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.position_reconstruction.lrfs_xml]);
    PosReconIq.pmt_gains = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.position_reconstruction.pmt_gain_xml]);
    PosRecon_Settings.pmts_off = gs_xml.data_processing_settings.global.pmts_off;
    PosRecon_Settings.file_lrfmat = [DataPath '\rec_fun.mat'];%This will cause problem

    ee = DatPFC_PosRecon(ee,PosRecon_Settings,PosReconIq,Setting.dp_module_setting.sr_version,Setting.debug.verbose);
end

%% Position Correction
if Setting.dp_module_setting.position_correction.run_module
    if Setting.debug.verbose
        disp('Running position correction.')
    end

    %PosCorr_XML = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.position_correction.setting_xml]);
    PosCorr_Settings = ModuleSetting.Corrections_PositionCorrection.parameters;
    PosCorrIq.xy_rec_cor = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.position_correction.xy_correction_xml]);

    ee = DatPFC_Corrections_PositionCorrection(ee,PosCorr_Settings,PosCorrIq,Setting.dp_module_setting.sr_version,Setting.debug.verbose);
end
%% Pulse area correction
if Setting.dp_module_setting.pulse_area_correction.run_module
    if Setting.debug.verbose
        disp('Running pulse area correction.')
    end


    %PulseCorr_XML = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.pulse_area_correction.setting_xml]);
    PulseCorr_Settings = ModuleSetting.Corrections_PositionCorrection.parameters;
    PulseCorr_Settings.Setting.wrapper.max_num_pulses = ee.rqs.num_pulses_found;
    PulseCorr_Settings.s2_class_to_correct = [2 4];
    PulseCorr_Settings.max_num_pulses = ee.rqs.num_pulses_found;
    PulseCorr_Settings.data_set = DataSet;
    
    PulseCorrIq.electron_lifetime_1 = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.pulse_area_correction.electron_lifetime_xml_1]);
    PulseCorrIq.electron_lifetime_2 = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.pulse_area_correction.electron_lifetime_xml_2]);
    PulseCorrIq.z_dep_s1_correction = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.pulse_area_correction.z_dep_s1_xml]);
    PulseCorrIq.s2_xy_correction = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.pulse_area_correction.s2_xy_xml]);
    PulseCorrIq.s1_xy_correction = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.pulse_area_correction.s1_xy_xml]);
    PulseCorrIq.s1_xyz_correction = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.pulse_area_correction.s1_xyz_xml]);

    ee = DatPFC_Corrections_ApplyCorrections(ee,PulseCorr_Settings,PulseCorrIq,Setting.dp_module_setting.sr_version,Setting.debug.verbose);
end
%% EnergyRecon
if Setting.dp_module_setting.energy_reconstruction.run_module
    if Setting.debug.verbose
        disp('Running energy reconstruction.')
    end

    %EnergyRecon_XML = XMLReader_framework([Setting.wrapper.xml_file_path Setting.dp_module_setting.energy_reconstruction.setting_xml]);
    EnergyRecon_Settings = ModuleSetting.EnergyReconstruction_Naive.parameters;
    EnergyRecon_Settings.max_num_pulses = ee.rqs.num_pulses_found;
    ee = DatPFC_EnergyReconstruction_Naive(ee,EnergyRecon_Settings,Setting.debug.verbose);
end





