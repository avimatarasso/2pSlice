function ROImask = compareROI(meanBL,ROImask,titleFile,justSave)
    answer = 'No, start over'; 
    while strcmp(answer, 'No, start over')
        if ~justSave  
            fig = figure(201);
            fig.Position = [100 100 1000 800];
            % Put both images into subplot
            subplot(1,2,1)
                imshow(meanBL, [], 'Colormap', gray(256));
            subplot(1,2,2)
            ha = imshow(meanBL, [], 'Colormap', gray(256));
            hold on;
            hb = imshow(ROImask, [], 'Colormap', gray(256));
            sgtitle(strrep(titleFile,'_',' '))
            % Set opacity/transparency to something less than 1 (alpha).  
            % 1 is the default and it means the last image is opaque and the image below can't be seen.
            hb.AlphaData = 0.1;

            quest    = 'Does this ROI look good or would you like to try again?';
            dlgtitle = 'ROI validation';
            btn1     = 'Yes, save now'; 
            btn2     = 'No, start over'; defbtn = btn2; 
            answer  = questdlg(quest,dlgtitle,btn1,btn2,defbtn); 
            clf
            if strcmp(answer,btn2)
                ROImask = createROImask(meanBL);                
            end
        else
            answer = 'Yes, save now';
        end
    end

end