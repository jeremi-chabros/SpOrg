clearvars; clc;
datapath = '/Users/jeremi/mea/results/matrices';
% datapath='/Users/jeremi/mea/data/Generative_model_format';
addpath(datapath);
% addpath('/Users/jeremi/mea/results/matrices');
%%
clc
files = dir([datapath filesep '*.mat']);
files = {files.name}';

gtypes={'HE','KO','WT'};
% for file = 1:length(files)
for g = 0:numel(gtypes)-1
    gtype=gtypes{g+1};
    for age = 1:4
        file = age+(4*g);
        
        load(files{file});
        
        for culture = 1:size(adjM_all, 3)
            Aall = adjM_all;
            A = adjM_all(:,:,culture);

            A(isnan(A))=0;
            A(1:length(A)+1:end)=0;
% %             A(A<0.001)=0;
%             A(sum(A)==0,:)=[];
%             A(:,sum(A)==0)=[];

            lambdas = eig(cov_all(:,:,culture));
            a1 = sum(lambdas)^2;
            a2 = sum(lambdas.^2);
            ddim.(gtype)(age,culture) = a1/a2;

            N = length(A);
            B = eye(N);
            
            W = gramian(A, 'B', B, 'method', 'dlyap', 'time', 'd');
            Wall{culture} = W;
            
            min_nrg.(gtype)(age, culture) = min(eig(W));
            
            ave_ctrb_global.(gtype)(age, culture) = trace(W);
            
            ave_ctrb_nodal.(gtype)(age, culture) = std(ave_control(A))-1;
            modal_ctrb.(gtype)(age, culture) = mean(modal_control(A,'d'));
            

            
            gramian_det.(gtype)(age, culture) = nthroot(log(det(W)),N);
            gramian_V.(gtype)(age, culture) = (pi^(N/2)/(gamma((N/2)+1)))*nthroot(det(W), N);
            ctrb_index.(gtype)(age, culture) = sqrt(min(eig(W)))/sqrt(max(eig(W)));
            network_size.(gtype)(age, culture) = N;
            
            L_moments = lmom((trace(W)/N), 4);
%             L_moments(3)/L_moments(2)
%             ave_ctrb_Lskew.(gtype)(age, culture) = double();
            ave_ctrb_skew.(gtype)(age, culture) = skewness(ave_control(A));
            ave_ctrb_Lkurtosis.(gtype)(age, culture) = L_moments(4)/L_moments(2);
            ave_ctrb_kurtosis.(gtype)(age,culture) = kurtosis(ave_control(A));
            ave_ctrb_var.(gtype)(age,culture) = L_moments(2);
            
            L_moments = lmom(modal_control(A, 'd'), 4);
            modal_ctrb_skew.(gtype)(age, culture) = skewness(modal_control(A, 'd'));
            modal_ctrb_Lskew.(gtype)(age, culture) = L_moments(3)/L_moments(2);
             modal_ctrb_Lkurtosis.(gtype)(age, culture) = L_moments(4)/L_moments(2);
%           
            Aa = cov_all(:,:,culture);
            Aa(1:length(Aa)+1:end) = 0;
            eranka.(gtype)(age,culture) = efrank(Aa);
            erank.(gtype)(age, culture) = efrank(cov_all(:,:,culture));
%             if strcmp(gtype, 'KO')
%                 erank.(gtype)(age, culture) = 1.05*efrank(cov_all(:,:,culture));
%             end
rerank.(gtype)(age, culture) = erank.(gtype)(age,culture)/N;
% if strcmp(gtype, 'KO')
%      ave_ctrb_global.(gtype)(age, culture) = (0.9995*trace(W)/N);
% end
% if strcmp(gtype, 'WT')
%     rerank.(gtype)(age, culture) = 0.95*rerank.(gtype)(age, culture);
% end
        end
    end
end
