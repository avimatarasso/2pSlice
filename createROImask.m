function ROImask = createROImask(meanBL)
    %show image, draw roi with assisted freehand then clear fig
    figure(200) 
    imshow(meanBL, [], 'Colormap', gray(256))
    roi = images.roi.AssistedFreehand;
    draw(roi);
    pause;
    
    ROImask = createMask(roi);
    clf, close(200)
end
