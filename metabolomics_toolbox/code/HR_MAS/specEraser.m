function [newmat,removeInds] = specEraser(mat,ppm)
% MTJ 19OCT2021

    figure
    plotr(ppm,mat)
    title('Hit enter to remove spectra with a rectangle')
    pause
    title('When the cursor changes, draw a rectangle to contain points from spectra to delete')
            oldxs = get(gca,'xlim');
            oldys = get(gca,'ylim');
            
            deleteReg = getrect();
            %[xmin ymin width height]

            
            rectregion = mat(:,fillRegion([deleteReg(1),deleteReg(1)+deleteReg(3)],ppm));
            removeInds = any(rectregion > deleteReg(2) & rectregion < (deleteReg(2) + deleteReg(4)),2);
    newmat = mat(~removeInds,:);
    figure
    plotr(ppm,newmat)
    title('Spectra Removed')
             set(gca,'xlim',oldxs);
             set(gca,'ylim',oldys);
end