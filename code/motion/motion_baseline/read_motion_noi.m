close all
clear
clc

% read in uniform timestep

% color specification
blue = [0.0000,0.4470,0.7410]; red = [0.8500,0.3250,0.0980];
gold = [0.9290,0.6940,0.1250]; teal = [32 178 170]/255;
green= [134, 179, 0]/255; purple = [153 102 255]/255;

% parameters
r = 1;
N = 80; rowN = 20; columnN = 4;
Nc = 4;
Ni = 4;
Li = 7;
T_MAX = 70;
uni_dt = 1;
Tnc = []; T = []; RP = []; AR = [];
% Tnc_avg = [];

% fixed labels
LB = ones(N,1);

% discretize spatial domain
% Lx = rowN * 2 * r * 4;      % length of the periodic domain
% Ly = columnN * 2 * r * 6;   % height of the periodic domain
Lx = 60;
Ly = 0;

% open file
filename = 'motion0_at1.dat';
fid = fopen(filename,'r');

figcell = figure;
% figTc = figure; 
% figTnc = figure;
time_count = 0; gen_num = 0;
t_block = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         initial state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read matrix nc
nc_array_init = fscanf(fid,'%d',N*N);
nc_init = reshape(nc_array_init,[N,N]);
% read cell label
LB_init = fscanf(fid,'%d',N);
vlen_init = N + pre_count(N,N+1,nc_init);
DATA_init = fscanf(fid,'%f',[2,vlen_init]);
DATA_init = DATA_init';

% plot initial state
gen_num = gen_num + 1;
figure(figcell);
% plot cells
for i = 1:N
    if (LB_init(i) == 0) % PST cell
        circle(DATA_init(i,1),DATA_init(i,2),r,red);
        %text(V(i,1),V(i,2),num2str(i));
    elseif (LB_init(i) == 1) % PSP cell
        circle(DATA_init(i,1),DATA_init(i,2),r,blue);
        %text(V(i,1),V(i,2),num2str(i));
    elseif (LB_init(i) == 2) % boundary cell
        circle(DATA_init(i,1),DATA_init(i,2),r,purple);
        %text(V(i,1),V(i,2),num2str(i));
    end
    hold on;
end
%  % plot c-sites
%  for i = 1:N
%      for j = (i+1):N
%          startP = N + pre_count(i,j,nc_init)+1;
%          endP = N + pre_count(i,j,nc_init)+nc_init(i,j);
%          colorC = [i/N,j/N,0];
%          plot(DATA_init(startP:endP,1),DATA_init(startP:endP,2),'*','Color',colorC,'MarkerSize',2);
%          hold on;
% %          for k = 1: nc_init(i,j)
% %              plot([DATA_init(i,1),DATA_init(N + pre_count(i,j,nc_init)+k,1)],...
% %                  [DATA_init(i,2),DATA_init(N + pre_count(i,j,nc_init)+k,2)],'Color','black');
% %              hold on;
% %              plot([DATA_init(j,1),DATA_init(N + pre_count(i,j,nc_init)+k,1)],...
% %                  [DATA_init(j,2),DATA_init(N + pre_count(i,j,nc_init)+k,2)],'Color','black');
% %              hold on;
% %          end
%      end
%  end
grid on;
xlim([0,Lx]);ylim([0,Ly]);
set(gca, 'FontSize', 28);
figname = [pwd '/2022/Aug/Aug11/','cell_', num2str(gen_num),'_T', num2str(time_count),'.png'];
saveas(figcell,figname);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      advance in time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (time_count < T_MAX)
    
    % read matrix nc
    nc_array = fscanf(fid,'%d',N*N);
    if isempty(nc_array)
        sprintf("the nc matrix is empty")
    end
    nc = reshape(nc_array,[N,N]);
    
    % read #of c-sites Tnc
    Tnc = [Tnc; fscanf(fid,'%d',1)];
    
%     % read #of cadherins of each cell
%     Tnc_single = fscanf(fid,'%d',N);
    
    % read matrix Tc
    % Tc_array = fscanf(fid,'%f',Nc * (N * (N-1))/2);
%     Tc_2d = reshape(Tc_array,[Nc,N*(N-1)/2]);
%     Tc_2d = Tc_2d';
%     Tc = zeros(N,N,Nc);
%     for i = 1:N
%         for j = (i+1):N
%             for k = 1:Nc
%                 row = N*(i-1) - i*(i-1)/2 + j - i;
%                 Tc(i,j,k) = Tc_2d(row,k);
%                 Tc(j,i,k) = Tc(i,j,k);
%             end
%         end
%     end
%     figure(figTc);
% %     S = size(Tc);
% %     [X,Y,Z] = ndgrid(1:S(1),1:S(2),1:S(3));
% %     scatter3(X(:),Y(:),Z(:),321,Tc(:),'filled')
%     H = heatmap(Tc(:,:,1));
%     figname = ['Tc_', num2str(gen_num), '_T', num2str(time_count), '.png'];
%     saveas(figTc,figname);

    % read time t0
    t0 = fscanf(fid,'%f',1);
    
    % read cell label
    LB = fscanf(fid,'%d',N);

    if (time_count >= 20)
        % read all the variables left
        vlen = N + pre_count(N,N+1,nc); % number of all the variables
        DATA = fscanf(fid,'%f',[2,vlen * 100]);
        DATA = DATA';

        % output information
        % sprintf("number of c-sites: %d", pre_count(N,N+1,nc))
        % sprintf("running time for this generation: %f", t0)
    
        ti = 100;
        % read the last block of data
        startR = (ti - 1) * vlen + 1;
        endR = startR + vlen - 1;
        V = DATA(startR:endR,:);


        % plotting ------------------------------------------------------
        if t_block < uni_dt
            t_block = t_block + t0;
        else
            t_block = 0;

            % plot cell motion
            figure(figcell);
            clf

            % plot cells
            for i = 1:N
                if (LB(i) == 0) % PST cell
                    circle(V(i,1),V(i,2),r,red);
                    %text(V(i,1),V(i,2),num2str(i));
                elseif (LB(i) == 1) % PSP cell
                    circle(V(i,1),V(i,2),r,blue);
                    %text(V(i,1),V(i,2),num2str(i));
                elseif (LB(i) == 2) % boundary cells
                    circle(V(i,1),V(i,2),r,purple);
                end
                hold on;
            end
%             % plot c-sites
%              for i = 1:N
%                  for j = (i+1):N
%                      startP = N + pre_count(i,j,nc)+1;
%                      endP = N + pre_count(i,j,nc)+nc(i,j);
%                      colorC = [i/N,j/N,0];
%                      plot(V(startP:endP,1),V(startP:endP,2),'*','Color',colorC,'MarkerSize',2);
%                      hold on;
%                  end
%              end
            grid on;
            xlim([0,Lx]);ylim([0,Ly]);
            set(gca, 'FontSize', 28);
            % figname = ['cell_', num2str(gen_num), '_T', num2str(time_count), '.png'];
            figname = [pwd '/2022/Aug/Aug11/', ...
                'cell_', num2str(gen_num), '_T', num2str(time_count), '.png'];
            saveas(figcell,figname);
            % generation number
            gen_num = gen_num + 1;
        end


        % relative position of red cells -----------------------------------
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
        T = [T; time_count];
        % Tnc_avg = [Tnc_avg, sum(Tnc_single) / N];
        RP = [RP; rp];
        AR = [AR; ar];
    else
        time_count = time_count + t0;
    end
end


% figure;
% plot(T, Tnc_avg, '-o');
% xlim([20,T_MAX]);
% ylim([0,10]);
% xlabel('time','FontSize',32)
% ylabel('average number of cadherins of each cell','FontSize',32)

figure;
subplot(2,1,1)
plot(T, RP, '-o');
xlim([20,T_MAX]);
ylim([0,1]);
xlabel('time','FontSize',32)
ylabel('relative position','FontSize',32)

subplot(2,1,2)
plot(T, AR, '-o');
xlim([20,T_MAX]);
ylim([0,1]);
xlabel('time','FontSize',32)
ylabel('Aspect Ratio','FontSize',32)