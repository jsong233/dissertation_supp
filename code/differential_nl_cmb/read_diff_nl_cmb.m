close all
clear
clc

% parameters
r = 1;
N = 80; rowN = 20; columnN = 4;
Nc = 6;
Ni = 4;
Li = 7;
uni_dt = 3 * 1e-3;

% fixed labels
LB = ones(N,1);

% discretize spatial domain
Lx = rowN * 2 * r * 4;      % length of the periodic domain
Ly = columnN * 2 * r * 6;   % height of the periodic domain

% standard time vector
T_MAX = 70; STEP = 10000; 
T = linspace(20,T_MAX,STEP);

% containers
H = [0,1,2,3,4]; % average loop
sampleIter = length(H);
NL = [0,1,2,3,4,5]; % parameter loop
nlIter = length(NL);
RP_cmb = zeros(sampleIter, STEP, nlIter);
RP_avg = zeros(STEP, nlIter);
AR_cmb = zeros(sampleIter, STEP, nlIter);
AR_avg = zeros(STEP, nlIter);


for nlselect = 1:nlIter
    % the nl parameter
    nl = NL(nlselect);
    
    % average over different initilizations under this nl
    for hselect = 1:sampleIter
        
        h = H(hselect);
        sprintf("nl loop = %d, h loop = %d", nl, h)
        t = []; RP = []; AR = [];

        % open file
        filename = "motion" + h + "_nl" + nl + ".dat";
        fid = fopen(filename,'r');

        %figcell = figure; 
        time_count = 0; % gen_num = 0;
        t_block = 0;

        % initial state
        % read matrix nc
        nc_array_init = fscanf(fid,'%d',N*N);
        nc_init = reshape(nc_array_init,[N,N]);
        % read cell label
        LB_init = fscanf(fid,'%d',N);
        vlen_init = N + pre_count(N,N+1,nc_init);
        DATA_init = fscanf(fid,'%f',[2,vlen_init]);
        DATA_init = DATA_init';


        % advance in time --------------------------------------------------
        while (time_count < T_MAX)
            % generation number
            % gen_num = gen_num + 1;

            % read matrix nc
            nc_array = fscanf(fid,'%d',N*N);
            if isempty(nc_array)
                sprintf("the nc matrix is empty")
            end
            nc = reshape(nc_array,[N,N]);

            % read #of c-sites Tnc
            Tnc = fscanf(fid,'%d',1);

            % read time t0
            t0 = fscanf(fid,'%f',1);

            % read cell label
            LB = fscanf(fid,'%d',N);

            if (time_count >= 20)
                % read all the variables left
                vlen = N + pre_count(N,N+1,nc); % number of all the variables
                DATA = fscanf(fid,'%f',[2,vlen * 100]);
                DATA = DATA';

                % read the last block of data
                ti = 100;
                startR = (ti - 1) * vlen + 1;
                endR = startR + vlen - 1;
                V = DATA(startR:endR,:);

                % relative position of red cells
                x0 = 0; red_count = 0;
                % find the frontmost and backmost position
                x1 = V(1,1); x2 = V(1,1);
                for i = 1:N
                    if (V(i,1) < x1)
                        x1 = V(i,1);
                    end
                    if (V(i,1) > x2)
                        x2 = V(i,1);
                    end
                end
                % find the highest and lowest position
                y1 = V(1,2); y2 = V(1,2);
                for i = 1:N
                    if (V(i,2) < y1)
                        y1 = V(i,2);
                    end
                    if (V(i,2) > y2)
                        y2 = V(i,2);
                    end
                end
                % calculate the center of mass of red cells
                for i = 1:N
                    if (LB(i) == 0)
                        x0 = x0 + V(i,1);
                        red_count = red_count + 1;
                    end
                end
                x0 = x0 / red_count;

                % calculate the relative position of red cells
                rp = (x0 - x1) / (x2 - x1);

                % calculate the aspect ratio of the slug hight/width
                ar = (y2 - y1) / (x2 - x1);

                % save information
                time_count = time_count + t0;
                t = [t; time_count];
                RP = [RP; rp];
                AR = [AR; ar];
            else
                time_count = time_count + t0;
            end
        end
        
        % close file
        fclose(fid);

        % save vectors
        outputname = "output" + h + "_nl" + nl + ".mat";
        save(outputname,'t','RP','AR');

        % align data
        for i = 1:STEP
            % find the entry in t closest to T(i)
            Comp = repmat(T(i),[length(t) 1]); % compare vector
            [~,j] = min(abs(t - Comp)); 
            RP_cmb(hselect,i,nlselect) = RP(j);
            AR_cmb(hselect,i,nlselect) = AR(j);
        end
    end
    
    % average over sample iterations
    for i = 1:STEP
        RP_avg(i,nlselect) = mean(RP_cmb(:,i,nlselect));
        AR_avg(i,nlselect) = mean(AR_cmb(:,i,nlselect));
    end
end


% set figure parameters
set(0,'DefaultLineLineWidth',1);

blue = [0.0000    0.4470    0.7410];
red = [0.8500    0.3250    0.0980];
gold = [0.9290    0.6940    0.1250];
teal = [32 178 170]/255;
green= [134, 179, 0]/255;
purple = [153 102 255]/255;

color = {blue red gold green teal purple};
lineSpec = {'-o','-^','-s','-*','-+','-d'};

% Plotting
% RP vs T 
figure(1)
for k = 1:nlIter
    eb = shadedErrorBar(T,RP_cmb(:,:,k),{@mean,@std},'lineprops',...
        {lineSpec{k},'markersize',8});
    eb.patch.FaceColor = color{k};
    eb.mainLineColor = color{k};
    set(eb.edge(1), 'Color', color{k}+(1-color{k})*0.5);
    set(eb.edge(2), 'Color', color{k}+(1-color{k})*0.5);
    hold on;
end
for k = 1:nlIter
    h(k) = plot(T,RP_avg(:,k),lineSpec{k},'markersize',8,'Color',color{k});
    hold on;
end
set(gca,'FontSize',24);
l = legend(h,'$nl = 1$','$nl = 1.5$','$nl = 2$','$nl = 2.5$','$nl = 3$','$nl = 3.5$');
set(l,'Interpreter','latex')
set(l,'FontSize',28);
set(l,'FontName','Times New Roman');
xlim([20,70]);
ylim([0,1]);
xlabel('time','FontSize',36)
ylabel('relative position','FontSize',36)
grid on;
pbaspect([2 1 1])


% AR vs T 
figure(2)
for k = 1:nlIter
    eb = shadedErrorBar(T,AR_cmb(:,:,k),{@mean,@std},'lineprops',...
        {lineSpec{k},'markersize',8});
    eb.patch.FaceColor = color{k};
    eb.mainLineColor = color{k};
    set(eb.edge(1), 'Color', color{k}+(1-color{k})*0.5);
    set(eb.edge(2), 'Color', color{k}+(1-color{k})*0.5);
    hold on;
end
for k = 1:nlIter
    h(k) = plot(T,AR_avg(:,k),lineSpec{k},'markersize',8,'Color',color{k});
    hold on;
end
set(gca,'FontSize',24);
l = legend(h,'$nl = 1$','$nl = 1.5$','$nl = 2$','$nl = 2.5$','$nl = 3$','$nl = 3.5$');
set(l,'Interpreter','latex')
set(l,'FontSize',28);
set(l,'FontName','Times New Roman');
xlim([20,120]);
ylim([0,1]);
xlabel('time','FontSize',36)
ylabel('Aspect Ratio','FontSize',36)
grid on;
pbaspect([2 1 1])