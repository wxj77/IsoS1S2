
function Level1_dsp1_1(iij)

DataSets{1,1}='lux10_20140903T2300';
addpath('/home/dkhaitan/IsolatedS1S2/r1');
RunSetting='/home/dkhaitan/IsolatedS1S2/r1/Settings_SR21_BH.xml';


cd '/home/dkhaitan/IsolatedS1S2/';
DataOutPath = ['/home/dkhaitan/IsolatedS1S2/output/'];
MaxDriftLength = 33000;


	DataSet=DataSets{1,iij};
	disp([DataSet])
	DataPath = ['/home/dkhaitan/LUXData/IsolatedS1S2/' DataSet '/'];
	FileList = dir([DataPath '*.dat']);
	FileListToRun = FileList;%datasample(FileList,NumberOfFiles,'Replace',false);

	Setting = XMLReader_framework(RunSetting);
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
	counter=0;

	for ii=1:3%len(FileListToRun)
		%try
			disp([ii,FileListToRun(ii).name]);
			[FileOut,ee_merge] = Level2(FileListToRun(ii).name,DataPath,DataOutPath,DataSet,Setting);
			[temp_is1r,temp_is2r,temp_is1s,temp_is2s,temp_s1r,temp_s2r,temp_s1s,temp_s2s,temp_lt,temp_ilt,temp_aqlt] = IsolateS1S2s_IncludeNew(ee_merge,MaxDriftLength);

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


		%catch exception
		%	disp(['Fail ' FileList(ii).name])
		%end
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
