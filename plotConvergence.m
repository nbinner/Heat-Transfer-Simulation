function plotConvergence()
% Read all .dat files in directory
FilesConv = dir('*-*conv*.dat'); % find all files that have - and conv and .dat

disp('Detected Files: '); 

for i = 1:length(FilesConv) %display detected files
    disp([ '   ' num2str(i) ': ' FilesConv(i).name]);
end

%%

for i = 1:length(FilesConv)
    tmp  = csvread(FilesConv(i).name); %set tmp to hold contents of file
    disp(FilesConv(i).name);

    %resMax    = log(tmp(:,1)); % Pulls the log of 1st column as a vectors from the C printout
    %resRMS    = log(tmp(:,2)); % Pulls the log of 2nd column as a vectors from the C printout
    resMax    = tmp(:,1); % Pulls the log of 1st column as a vectors from the C printout
    resRMS    = tmp(:,2); % Pulls the log of 2nd column as a vectors from the C printout
    Iteration = tmp(:,3); % Pulls the 3rd column as a vectors from the C printout
   
    %plot(Iteration,resMax,'-b',Iteration,resRMS,'-r'); %plots resRMS and resMAX against Iteration
    semilogy(Iteration,resMax,'-b',Iteration,resRMS,'-r'); %plots resRMS and resMAX against Iteration
    legend('Max Residual','RMS Residual'); % Adds legend for each set
    title(strrep(FilesConv(i).name,'.dat','')); %sets title to file name
    xlabel('Iteration','FontSize',14); % adds x-axis label
    ylabel('Log(Residual)','FontSize',14); % adds y-axis label
    set(gca,'XGrid','on');
    set(gca,'YGrid','on');
 
    pause(0.1);
  
    [~,fname,~] = fileparts(FilesConv(i).name); %pull filename from full path
     disp(['   Saving into:' fname]); %display (debugging)
    % adds _Results.tiff to the end of the file name and saves it
     saveas(gcf,[ fname '.png']);    
end
end