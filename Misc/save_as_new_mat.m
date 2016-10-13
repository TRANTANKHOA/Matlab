%% Initialization
clear; close all; clc;  % Clear and close all
% ord = 1;                % From 1 to 5
ERR = 'G';              % Could be 'G' or 'U'
for k=2:5
    ord = k;
    %% Metropolis
    OPT.sampler='metropolis';
    filename2 = [ERR num2str(ord) '_' OPT.sampler];
    load([filename2 '_0' '.mat'],'OPT')
    for i=1:OPT.div
        load([OPT.filename '_' num2str(i) '.mat'],'TH');
        save([OPT.filename '_' num2str(i) '.mat'],'TH','-v7.3');
    end
    clear OPT
    %% VA
    OPT.sampler='va';
    filename2 = [ERR num2str(ord) '_' OPT.sampler];
    load([filename2 '_0' '.mat'],'OPT')
    for i=1:OPT.div
        load([OPT.filename '_' num2str(i) '.mat'],'TH');
        save([OPT.filename '_' num2str(i) '.mat'],'TH','-v7.3');
    end
    clear OPT
    %% Slice
    OPT.sampler='slice';
    filename2 = [ERR num2str(ord) '_' OPT.sampler];
    load([filename2 '_0' '.mat'],'OPT')
    for i=1:OPT.div
        load([OPT.filename '_' num2str(i) '.mat'],'TH');
        save([OPT.filename '_' num2str(i) '.mat'],'TH','-v7.3');
    end
    clear OPT
    %% Parslice
    % OPT.sampler='parslice';
    % filename2 = [ERR num2str(ord) '_' OPT.sampler];
    % load([filename2 '_0' '.mat'],'OPT')
    % for i=1:OPT.div
    %     load([OPT.filename '_' num2str(i) '.mat'],'TH');
    %     save([OPT.filename '_' num2str(i) '.mat'],'TH','-v7.3');
    % end
    % clear OPT
end
ERR = 'U';              % Could be 'G' or 'U'
for k=2:5
    ord = k;
    %% Metropolis
    OPT.sampler='metropolis';
    filename2 = [ERR num2str(ord) '_' OPT.sampler];
    load([filename2 '_0' '.mat'],'OPT')
    for i=1:OPT.div
        load([OPT.filename '_' num2str(i) '.mat'],'TH');
        save([OPT.filename '_' num2str(i) '.mat'],'TH','-v7.3');
    end
    clear OPT
    %% VA
    OPT.sampler='va';
    filename2 = [ERR num2str(ord) '_' OPT.sampler];
    load([filename2 '_0' '.mat'],'OPT')
    for i=1:OPT.div
        load([OPT.filename '_' num2str(i) '.mat'],'TH');
        save([OPT.filename '_' num2str(i) '.mat'],'TH','-v7.3');
    end
    clear OPT
    %% Slice
    OPT.sampler='slice';
    filename2 = [ERR num2str(ord) '_' OPT.sampler];
    load([filename2 '_0' '.mat'],'OPT')
    for i=1:OPT.div
        load([OPT.filename '_' num2str(i) '.mat'],'TH');
        save([OPT.filename '_' num2str(i) '.mat'],'TH','-v7.3');
    end
    clear OPT
    %% Parslice
    % OPT.sampler='parslice';
    % filename2 = [ERR num2str(ord) '_' OPT.sampler];
    % load([filename2 '_0' '.mat'],'OPT')
    % for i=1:OPT.div
    %     load([OPT.filename '_' num2str(i) '.mat'],'TH');
    %     save([OPT.filename '_' num2str(i) '.mat'],'TH','-v7.3');
    % end
    % clear OPT
end