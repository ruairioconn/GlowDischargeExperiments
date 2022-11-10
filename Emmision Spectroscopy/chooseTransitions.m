function [NIST_lambda, data_lambda, NIST_ix, data_ix, diff, pks, locs, wdths, A] = chooseTransitions(x,y,varargin)

diff = 0;
if exist("varargin","var")
    [pks,locs] = findpeaks(y,x,'MinPeakProminence',0.01*max(y));
    ypks = ones(size(locs))*(-0.05*max(y));
    NIST_data = call_NIST(x,"Table");
    yy = ones(size(NIST_data.lambda))*(-0.1*max(y));

    figure %Plot
    hold on
    plot(x,y-min(y))
%     gscatter(NIST_data.lambda,yy,NIST_data.sp_num)
    scatter(NIST_data.lambda,yy)
    scatter(locs,ypks,'filled','black')
%     legend('Data','Ar I','Ar II','Ar III','AR IV','Data Peaks','AutoUpdate','off')
    legend('Data','Ar I')

    title('Calibration', '\color{blue}Select NIST Transition (Colored Marks)')
    [xg,yg,b] = ginput(1);
    D = pdist2([NIST_data.lambda yy],[xg yg]);
    [~,ixmin] = min(D);
    scatter(NIST_data.lambda(ixmin),yy(ixmin),'MarkerEdgeColor','green')
    text(NIST_data.lambda(ixmin),yy(ixmin),num2str(NIST_data.lambda(ixmin)))
    NIST_lambda = NIST_data.lambda(ixmin);

    title('Calibration', '\color{black}Select Data Transition (Black Marks)')
    [xg,yg,b] = ginput(1);
    D = pdist2([locs' ypks'],[xg yg]);
    [~,ixmin] = min(D);
    scatter(locs(ixmin),ypks(ixmin),'MarkerEdgeColor','green')
    data_lambda = locs(ixmin);
    diff = (data_lambda-NIST_lambda);

    x = x - diff;
end

close all

%% Main Function

[pks,locs,wdths] = findpeaks(y,x,'MinPeakProminence',0.01*max(y));
ypks = ones(size(locs))*(-0.05*max(y));
NIST_data = call_NIST(x,"Table");
yy = ones(size(NIST_data.lambda))*(-0.1*max(y));

figure %Plot
hold on
plot(x,y-min(y))
% gscatter(NIST_data.lambda,yy,NIST_data.sp_num)
scatter(NIST_data.lambda,yy)
scatter(locs,ypks,'filled','black')
% legend('Data','Ar I','Ar II','Ar III','AR IV','Data Peaks','AutoUpdate','off')
legend('Data','Ar I')
set (gcf, 'WindowButtonMotionFcn', @mouseMove);

% NIST Transitions
title('\color{blue}Select NIST Transitions (Colored Marks)')
ix = [];
while 1>0
    [xg,yg,b] = ginput(1);
    if isempty(b)
        % Exit loop if
        break
    elseif b==45
        ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
        axis([xg-width/2 xg+width/2 yg-height/2 yg+height/2]);
        zoom(1/2);
    elseif b==61
        ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
        axis([xg-width/2 xg+width/2 yg-height/2 yg+height/2]);
        zoom(2);
    elseif b==28
        ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
        axis([ax(1)-width/5 ax(2)-width/5 ax(3) ax(4)]);
    elseif b==29
        ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
        axis([ax(1)+width/5 ax(2)+width/5 ax(3) ax(4)]);
    elseif b==30
        ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
        axis([ax(1) ax(2) ax(3)+height/5 ax(4)+height/5]);
    elseif b==31
        ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
        axis([ax(1) ax(2) ax(3)-height/5 ax(4)-height/5]);
    else
        D = pdist2([NIST_data.lambda yy],[xg yg]);
        [~,ixmin] = min(D);
        if ismember(ixmin,ix)
            scatter(NIST_data.lambda(ix(ix==ixmin)),yy(ix(ix==ixmin)),'MarkerEdgeColor','white')
            ix = ix(ix~=ixmin);
        else
            ix = [ix,ixmin];
            scatter(NIST_data.lambda(ix(ix==ixmin)),yy(ix(ix==ixmin)),'MarkerEdgeColor','green')
            text(NIST_data.lambda(ixmin),yy(ixmin),num2str(round(NIST_data.lambda(ixmin),1)))
        end
    end
end
ix = sort(ix);
NIST_lambda = NIST_data.lambda(ix);
NIST_ix = ix;
A = NIST_data.A_ki(NIST_ix)';

clear ix
% Data Transitions
title('\color{black}Select Data Transitions (Black Marks)')
ix = [];
while 1>0
    [xg,yg,b] = ginput(1);
    if isempty(b)
        % Exit loop if
        break
    elseif b==45
        ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
        axis([xg-width/2 xg+width/2 yg-height/2 yg+height/2]);
        zoom(1/2);
    elseif b==61
        ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
        axis([xg-width/2 xg+width/2 yg-height/2 yg+height/2]);
        zoom(2);
    elseif b==28
        ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
        axis([ax(1)-width/5 ax(2)-width/5 ax(3) ax(4)]);
    elseif b==29
        ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
        axis([ax(1)+width/5 ax(2)+width/5 ax(3) ax(4)]);
    elseif b==30
        ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
        axis([ax(1) ax(2) ax(3)+height/5 ax(4)+height/5]);
    elseif b==31
        ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
        axis([ax(1) ax(2) ax(3)-height/5 ax(4)-height/5]);
    else
        D = pdist2([locs' ypks'],[xg yg]);
        [~,ixmin] = min(D);
        if ismember(ixmin,ix)
            scatter(locs(ix(ix==ixmin)),ypks(ix(ix==ixmin)),'MarkerEdgeColor','white')
            ix = ix(ix~=ixmin);
        else
            ix = [ix,ixmin];
            scatter(locs(ix(ix==ixmin)),ypks(ix(ix==ixmin)),'MarkerEdgeColor','green')
        end
    end
end
ix = sort(ix);
data_lambda = locs(ix)';
data_ix = ix;
pks = pks(data_ix);
locs = locs(data_ix);
wdths = wdths(data_ix);
close all

function mouseMove (object, eventdata)
C = get (gca, 'CurrentPoint');
title(gca, ['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);
end
end