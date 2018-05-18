
function Level1_r3_sbatch_1(iij)

pc=parcluster('local');
JOB_ID=getenv('SLURM_JOBID');
CPUS=str2num(getenv('SLURM_JOB_CPUS_PER_NODE'));
pc.JobStorageLocation=strcat('/local_scratch/',JOB_ID);
parpool(pc,CPUS)
DataSets{1,1}='lux10_20151225T0621';
DataSets{1,2}='lux10_20151225T1415';
DataSets{1,3}='lux10_20151225T2205';
DataSets{1,4}='lux10_20151226T0559';
DataSets{1,5}='lux10_20151226T1355';
DataSets{1,6}='lux10_20151226T2143';
DataSets{1,7}='lux10_20151227T0051';
DataSets{1,8}='lux10_20151227T0838';
DataSets{1,9}='lux10_20151227T1629';
DataSets{1,10}='lux10_20151228T0020';
DataSets{1,11}='lux10_20151228T0943';
DataSets{1,12}='lux10_20151228T1517';
DataSets{1,13}='lux10_20151228T2308';
DataSets{1,14}='lux10_20151229T0657';
DataSets{1,15}='lux10_20151229T1647';
DataSets{1,16}='lux10_20151230T0033';
set(gcf,'Visible','off')
addpath('/scratch/dkhaitan/IsolatedS1S2/r1');
RunSetting='/scratch/dkhaitan/IsolatedS1S2/r1/Settings_SR21_BH.xml';


cd '/scratch/dkhaitan/IsolatedS1S2/r1/';
DataOutPath = ['/scratch/dkhaitan/IsolatedS1S2/s4_output/'];
MaxDriftLength = 35000;

	DataSet=DataSets{1,iij};
	disp([DataSet])
	DataPath = ['/scratch/dkhaitan/LUXData/IsolatedS1S2/' DataSet '/'];
	FileList = dir([DataPath '*.dat']);
	FileListToRun = FileList;%datasample(FileList,NumberOfFiles,'Replace',false);

	Setting = XMLReader_framework(RunSetting);
	temp_isolated_s1_rate = 0;
	isolated_s2_rate = 0;
	s1_spect = zeros(1,1000);
	s2_spect = zeros(1,1000);
	livetime=0;
	counter=0;


	parfor ii=1:length(FileListToRun)
		try
			disp([ii,FileListToRun(ii).name]);
			[FileOut,ee_merge] = Level2(FileListToRun(ii).name,DataPath,DataOutPath,DataSet,Setting);
			[temp_is1r,temp_is2r,temp_is1s,temp_is2s,temp_s1r,temp_s2r,temp_s1s,temp_s2s,temp_lt,temp_ilt,temp_aqlt] = IsolateS1S2s_OneS4(ee_merge,MaxDriftLength);

			fid_s1_spect=fopen([DataOutPath DataSet 'iso_s1_spect.csv'],'at');
			fprintf(fid_s1_spect,'%s,',num2str(temp_is1r));
			fprintf(fid_s1_spect,'%s,',num2str(temp_ilt));
			fprintf(fid_s1_spect,'%s,',num2str(temp_aqlt));
			dlmwrite([DataOutPath DataSet 'iso_s1_spect.csv'],temp_is1s,'-append','delimiter',',');
			fclose(fid_s1_spect);

			fid_s1_spect=fopen([DataOutPath DataSet 's1_spect.csv'],'at');
			fprintf(fid_s1_spect,'%s,',num2str(temp_s1r));
			fprintf(fid_s1_spect,'%s,',num2str(temp_lt));
			dlmwrite([DataOutPath DataSet 's1_spect.csv'],temp_s1s,'-append','delimiter',',');
			fclose(fid_s1_spect);

			fid_s2_spect=fopen([DataOutPath DataSet 'iso_s2_spect.csv'],'at');
			fprintf(fid_s2_spect,'%s,',num2str(temp_is2r));
			fprintf(fid_s2_spect,'%s,',num2str(temp_ilt));
			fprintf(fid_s2_spect,'%s,',num2str(temp_aqlt));
			dlmwrite([DataOutPath DataSet 'iso_s2_spect.csv'],temp_is2s,'-append','delimiter',',');
			fclose(fid_s2_spect);

			fid_s1_spect=fopen([DataOutPath DataSet 's2_spect.csv'],'at');
			fprintf(fid_s2_spect,'%s,',num2str(temp_s2r));
			fprintf(fid_s2_spect,'%s,',num2str(temp_lt));
			dlmwrite([DataOutPath DataSet 's2_spect.csv'],temp_s1s,'-append','delimiter',',');
			fclose(fid_s1_spect);


		catch exception
			disp(['Fail ' FileList(ii).name])
		end
		%clearvars -global -except DataSets DataSet DataPath DataOutPath MaxDriftLength FileList FileListToRun NumberOfFiles Setting isolated_s1_rate isolated_s2_rate s1_spect s2_spect livetime counter
	end
	iso_s1_rate = 0;
	iso_s2_rate = 0;
	iso_s1_spect = zeros(1,1000);
	iso_s2_spect = zeros(1,1000);
	s1_rate = 0;
	s2_rate = 0;
	s1_spect = zeros(1,1000);
	s2_spect = zeros(1,1000);
	livetime=0;
	ilt = 0;
	aqlt = 0;

	disp([DataOutPath DataSet 'iso_s1_spect.csv']);
	dummy = csvread([DataOutPath DataSet 'iso_s1_spect.csv']);
	dummy = sum(dummy);
	iso_s1_rate = dummy(1);
	ilt = dummy(2);
	aqlt = dummy(3);
	iso_s1_spect = dummy(4:length(dummy));

	dummy = csvread([DataOutPath DataSet 's1_spect.csv']);
	dummy = sum(dummy);
	s1_rate = dummy(1);
	lt = dummy(2);
	s1_spect = dummy(3:length(dummy));

	dummy = csvread([DataOutPath DataSet 'iso_s2_spect.csv']);
	dummy = sum(dummy);
	iso_s2_rate = dummy(1);
	ilt = dummy(2);
	aqlt = dummy(3);
	iso_s2_spect = dummy(4:length(dummy));

	dummy = csvread([DataOutPath DataSet 's2_spect.csv']);
	dummy = sum(dummy);
	s2_rate = dummy(1);
	lt = dummy(2);
	s2_spect = dummy(3:length(dummy));

	fid_s1_spect=fopen([DataOutPath 'total_iso_s1_spect_real1.csv'],'at');
	fprintf(fid_s1_spect,'%s,',DataSet);
	fprintf(fid_s1_spect,'%s,',num2str(iso_s1_rate));
	fprintf(fid_s1_spect,'%s,',num2str(ilt));
	fprintf(fid_s1_spect,'%s,',num2str(aqlt));
	dlmwrite([DataOutPath 'total_iso_s1_spect_real1.csv'],iso_s1_spect,'-append','delimiter',',');
	fclose(fid_s1_spect);

	fid_s1_spect=fopen([DataOutPath 'total_s1_spect_real1.csv'],'at');
	fprintf(fid_s1_spect,'%s,',DataSet);
	fprintf(fid_s1_spect,'%s,',num2str(s1_rate));
	fprintf(fid_s1_spect,'%s,',num2str(lt));
	dlmwrite([DataOutPath 'total_s1_spect_real1.csv'],s1_spect,'-append','delimiter',',');
	fclose(fid_s1_spect);


	fid_s2_spect=fopen([DataOutPath 'total_iso_s2_spect_real1.csv'],'at');
	fprintf(fid_s2_spect,'%s,',DataSet);
	fprintf(fid_s2_spect,'%s,',num2str(iso_s2_rate));
	fprintf(fid_s2_spect,'%s,',num2str(ilt));
	fprintf(fid_s2_spect,'%s,',num2str(aqlt));
	dlmwrite([DataOutPath 'total_s2_spect_real1.csv'],iso_s2_spect,'-append','delimiter',',');
	fclose(fid_s2_spect);

	fid_s2_spect=fopen([DataOutPath 'total_s2_spect_real1.csv'],'at');
	fprintf(fid_s2_spect,'%s,',DataSet);
	fprintf(fid_s2_spect,'%s,',num2str(s2_rate));
	fprintf(fid_s2_spect,'%s,',num2str(lt));
	dlmwrite([DataOutPath 'total_iso_s2_spect_real1.csv'],s2_spect,'-append','delimiter',',');
	fclose(fid_s2_spect);

	delete([DataOutPath DataSet 's2_spect.csv'])
	delete([DataOutPath DataSet 's1_spect.csv'])
	delete([DataOutPath DataSet 'iso_s2_spect.csv'])
	delete([DataOutPath DataSet 'iso_s1_spect.csv'])

	
end
