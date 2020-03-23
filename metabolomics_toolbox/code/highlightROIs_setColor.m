function [] = highlightROIs_setColor(ROIs,height)

%% highlightROIs
%{
    This function highlights the ppm regions provided in 'ROIs'.

    Hope this helps,
    MJ 2MAY2017
%}

%if exist('color','var')==0
    color = uisetcolor;
%end
   
    %ROIs = features_manual;
            for i = 1:size(ROIs,2)
                                        %llowerpt = [ROIs(1,i),0];
                                        width = (ROIs(2,i)-ROIs(1,i));
                                        %height = max(max(XR));
                                        llx = ROIs(1,i); % 
                                        lux = llx;
                                        rlx = llx + width; % = rux
                                        rux = rlx;
                                        lly = 0;
                                        luy = height;
                                        rly = 0;
                                        ruy = luy;
                                        xcoords = [llx rlx rux lux];
                                        ycoords = [lly rly ruy luy];
                                        %plot(llowerpt(1),llowerpt(2),'r*')
                                        %plot(llowerpt(1)+width,height,'r*')
                                        %rectangle('Position',[llowerpt,width,height]);
                                        p=patch(xcoords,ycoords,color); %light red % pink [1    0.5  1]
                                        set(p,'FaceAlpha',0.1,'EdgeColor','none'); % transparency
            end
end