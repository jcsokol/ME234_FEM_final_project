function Folder = CreateOutputFolders()
    Folder.Main  = pwd;

    % Compute new folder name
    rep_tot = 1;
    while(isdir([Folder.Main,'/RES_',num2str(rep_tot,'%03d')]))
        rep_tot = rep_tot + 1;
    end
    
    % Create folders
    Folder.Simulation = [Folder.Main,'/RES_',num2str(rep_tot,'%03d')];
    Folder.Output = [Folder.Simulation,'/Output'];

    mkdir(Folder.Simulation);
    mkdir(Folder.Output);
    
    % Copy all .m files
    % mFiles = ls('*.m');
    % mFiles = textscan(mFiles,'%s'); mFiles = mFiles{1};
    % for id=1:length(mFiles)
    %     copyfile(mFiles{id},[Folder.Simulation,'/',mFiles{id}]);
    % end
end