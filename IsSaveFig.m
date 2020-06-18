function IsSaveFig(figName, id, filename, savedir)

    cnt = questdlg('Save Figure?', ...
	'Save Figure', ...
	'Yes','No','No');
    % Handle response
    switch cnt
        case 'Yes'
            disp(['Saving the figure to:', savedir])
            saveas(figure(figName),[savedir,'/',id,filename,'.fig']);
        case 'Exit'
            disp([' ', savedir])
    end
    
close all
end

