

DataSets{1,1}='lux10_20140903T1918';
DataSets{1,2}='lux10_20140903T2300';
RunSetting='/scratch/dkhaitan/IsolatedS1S2/r1/Settings_SR21_BH.xml';

slurm_profile=mdce_profile();
nworkers=str2num(getenv('SLURM_NTASKS'))
parpool(slurm_profile, nworkers)

addpath('/scratch/dkhaitan/IsolatedS1S2/r1');
cd '/scratch/dkhaitan/IsolatedS1S2/';
DataOutPath = ['/scratch/dkhaitan/IsolatedS1S2/test_out/'];
MaxDriftLength = 35000;


fid_s2_spect=fopen([DataOutPath 'total_s2_spect.csv'],'at');
fprintf(fid_s2_spect,'Acquisition,Rate,');
for kk =1:1000
fprintf(fid_s2_spect,'%s,',num2str((kk-1)*10));
end
fclose(fid_s2_spect);

fid_s1_spect=fopen([DataOutPath 'total_s1_spect.csv'],'at');
fprintf(fid_s1_spect,'Acquisition,Rate,');
for kk =1:1000
fprintf(fid_s1_spect,'%s,',num2str((kk-1)*1));
end
fclose(fid_s1_spect);


disp(size(DataSets,2));
for dn = 1:size(DataSets,2)
	DataSet=DataSets{1,dn};
	disp([DataSet])
	DataPath = ['/scratch/dkhaitan/LUXData/IsolatedS1S2/' DataSet '/'];
	FileList = dir([DataPath '*.dat']);
	NumberOfFiles = 1000;% ceil(length(FileList)/100);
	FileListToRun = datasample(FileList,NumberOfFiles,'Replace',false);

	if ~exist(strrep(DataPath,'dat','DatPFC_Extended'),'dir')
	    mkdir(strrep(DataPath,'dat','DatPFC_Extended'));
	end

	Setting = XMLReader_framework(RunSetting);
	temp_isolated_s1_rate = 0;
	isolated_s2_rate = 0;
	s1_spect = zeros(1,1000);
	s2_spect = zeros(1,1000);
	livetime=0;
	counter=0;


	%Setting = xmlread(RunSetting);
	parfor ii = 1:NumberOfFiles
		try
			%disp([ii,FileListToRun(ii).name]);
			[FileOut,ee_merge] = Level2(FileListToRun(ii).name,DataPath,DataOutPath,DataSet,Setting);
			[temp_isolated_s1_rate,temp_isolated_s2_rate,temp_s1_spect,temp_s2_spect,temp_livetime] = IsolateS1S2s(ee_merge,MaxDriftLength);
			%{
			isolated_s1_rate = isolated_s1_rate + temp_isolated_s1_rate;
			isolated_s2_rate = isolated_s1_rate + temp_isolated_s2_rate;
			s1_spect = s1_spect + temp_s1_spect
			s2_spect = s2_spect + temp_s2_spect
			livetime = livetime+temp_livetime
			%}
			%{
			isolated_s1_rate = (livetime_total*isolated_s1_rate + livetime*temp_isolated_s1_rate)/(livetime_total + livetime);
			isolated_s2_rate = (livetime_total*isolated_s2_rate + livetime*temp_isolated_s2_rate)/(livetime_total + livetime);
			s1_spect = (s1_spect*livetime_total + livetime*temp_s1_spect)/(livetime_total + livetime);
			s2_spect = (s2_spect*livetime_total + livetime*temp_s2_spect)/(livetime_total + livetime);
			livetime_total=livetime_total + livetime;
			%}
			%disp([ii,temp_isolated_s1_rate,temp_isolated_s2_rate,temp_livetime]);
			fid_s2_spect=fopen([DataOutPath DataSet 's2_spect.csv'],'at');
			fprintf(fid_s2_spect,'%s,',num2str(temp_isolated_s2_rate));
			fprintf(fid_s2_spect,'%s,',num2str(temp_livetime));
			dlmwrite([DataOutPath DataSet 's2_spect.csv'],temp_s2_spect,'-append','delimiter',',');
			fclose(fid_s2_spect);

			fid_s1_spect=fopen([DataOutPath DataSet 's1_spect.csv'],'at');
			fprintf(fid_s1_spect,'%s,',num2str(temp_isolated_s1_rate));
			fprintf(fid_s2_spect,'%s,',num2str(temp_livetime));
			dlmwrite([DataOutPath DataSet 's1_spect.csv'],temp_s1_spect,'-append','delimiter',',');
			fclose(fid_s1_spect);
		catch exception
			disp(['Fail ' FileList(ii).name])
		end
		%clearvars -global -except DataSets DataSet DataPath DataOutPath MaxDriftLength FileList FileListToRun NumberOfFiles Setting isolated_s1_rate isolated_s2_rate s1_spect s2_spect livetime counter
	end
	dummy = csvread([DataOutPath DataSet 's2_spect.csv']);
	dummy = sum(dummy);
	isolated_s2_rate = dummy(1)/dummy(2);
	s2_spect = dummy(3:length(dummy))/dummy(2);
	dummy = csvread([DataOutPath DataSet 's1_spect.csv']);
	dummy = sum(dummy);
	isolated_s1_rate = dummy(1)/dummy(2);
	s1_spect = dummy(3:length(dummy))/dummy(2);

	fid_s2_spect=fopen([DataOutPath 'total_s2_spect.csv'],'at');
	fprintf(fid_s2_spect,'%s,',DataSet);
	fprintf(fid_s2_spect,'%s,',num2str(isolated_s2_rate));
	dlmwrite([DataOutPath 'total_s2_spect.csv'],s2_spect,'-append','delimiter',',');
	fclose(fid_s2_spect);

	fid_s1_spect=fopen([DataOutPath 'total_s1_spect.csv'],'at');
	fprintf(fid_s1_spect,'%s,',DataSet);
	fprintf(fid_s1_spect,'%s,',num2str(isolated_s1_rate));
	dlmwrite([DataOutPath 'total_s1_spect.csv'],s1_spect,'-append','delimiter',',');
	fclose(fid_s1_spect);
	delete([DataOutPath DataSet 's2_spect.csv'])
	delete([DataOutPath DataSet 's1_spect.csv'])
end

