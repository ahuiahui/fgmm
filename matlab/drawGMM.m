function h=drawGMM(GMM,dims,color)
%function drawGMM(GMM,dims,color)
%draw the specified dimensions of the GMM and +-3stddev

%TODO:  Alpha comps with prior
%alpha points with weight
if(isempty(GMM))
    return
end
if(nargin<2)
    dims=1:2;
end

if(nargin<3)
    color = [0 0 0];
end

if(ischar(color))%catch plot style strings
    disp('Do not use strings for color, use a triplet')
    color=[0 0 0];
end
%trap for backwards...remove
if(~isfield(GMM,'Type'))
    GMM.Type='Trained';
end
if strcmp(GMM.Type,'Dual')
    drawGMM(GMM.Pos,dims,color+[0 1 0]);
    drawGMM(GMM.Neg,dims,color+[1 0 0]);
    return
end

Prior=GMM.Prior;
Mu=GMM.Mu;
Sigma=GMM.Sigma;
dataset=GMM.data;

color(color>1)=1;
color(color<0)=0;

if(isempty(Mu))
    return;
end

if(isempty(GMM.vars))
    GMM.vars=ones(max(dims),1);
end
if(isempty(GMM.means))
    GMM.means=zeros(max(dims),1);
end

dim=length(dims);
held=ishold;

% if(dim==3)
%     for nt=1:size(Mu,2);
%         [x y z]=sphere;
%         x = reshape(x,441,1);
%         y = reshape(y,441,1);
%         z = reshape(z,441,1);
%         %marginals
%         stdev= sqrtm(Sigma(dims,dims,nt));
%         t = [x y z]*3*real(stdev);
%         x = Mu(dims(1),nt) + reshape(t(:,1),21,21);
%         y = Mu(dims(2),nt) + reshape(t(:,2),21,21);
%         z = Mu(dims(3),nt) + reshape(t(:,3),21,21);
%         s=surf(x,y,z);
%         set(s,'FaceColor',lightcolor,'EdgeColor','none','FaceAlpha',.5)
%         hold on
%     end
%     plot3(Mu(dims(1),:), Mu(dims(2),:), Mu(dims(3),:), 'x', 'lineWidth', 2, 'color', color);
%     if(~isempty(dataset))
%         plot3(dataset(:,dims(1)),dataset(:,dims(2)),dataset(:,dims(3)),'k.');
%     end
% end
if(dim==2)
    numSegments = 40;
    t = linspace(-pi, pi, numSegments);
    for nt=1:size(Mu,2)
        stdev = 3*sqrtm(Sigma(dims,dims,nt));
        X = real(stdev)*[cos(t);sin(t)]+repmat(Mu(dims,nt),1,numSegments);
        %unsphere
        X=X.*repmat(sqrt(GMM.vars(dims)),1,numSegments);
        X=X+repmat(GMM.means(dims),1,numSegments);
        %these may be combineable
        % h=patch(X(1,:), X(2,:), color, 'lineWidth', 1, 'EdgeColor', ...
        %         [0 0 0]);
        % set(h,'FaceAlpha',Prior(nt));%,'EdgeAlpha',0.5);
	plot(X(1,:),X(2,:),'Color',color,'lineWidth',3);
        hold on
    end
    mm=Mu.*repmat(sqrt(GMM.vars),1,size(Mu,2));
    mm=mm+repmat(GMM.means,1,size(Mu,2));
    %plot(mm(dims(1),:), mm(dims(2),:), 'x', 'lineWidth', 2, 'color', color);
    for nt=1:size(Mu,2)
%         text(mm(dims(1),nt),mm(dims(2),nt),num2str(nt),'FontSize',12,'FontWeight','Bold','Color',1-color,'BackgroundColor',color);
    end
    if(~isempty(dataset))
        plot(dataset(dims(1),:),dataset(dims(2),:),'.','Color',color);
    end
end

if held
    hold on
else
    hold off
end
