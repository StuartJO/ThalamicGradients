function SigMapDescrips = PlotSigNeuromaps(neuromap_corrs,outdir,xlabel_var)

if nargin < 1
    load('./data/processed/NeuroMapCorrs.mat','neuromap_corrs')
end

if nargin < 2
    outdir = './figure_outputs/SigNeuromaps';
end

if nargin < 3
   xlabel_var = ['Cortical region PC1 loading'];
end
    
mkdir(outdir)

sig_maps = find(neuromap_corrs.corr_sig);

annotlabels = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};

for i = 1:length(sig_maps)
    figure('Position',[291   391   560   420])
    SigInd = sig_maps(i);
    if neuromap_corrs.p_perm(SigInd) == 0
     pval_format = '{\itp_{spin}} < .0001 ';   
    else
    pval_format = ['{\itp_{spin}} = ',num2str(round(neuromap_corrs.p_perm(SigInd),4))];
    end
    s = scatterfit(neuromap_corrs.data,neuromap_corrs.neuromap_parc(:,SigInd),40,neuromap_corrs.data,[],0);
    %title({neuromap_corrs.name{SigInd},['({\itr} = ',num2str(neuromap_corrs.corr(SigInd),3),', ',pval_format,')']});
    title(['{\itr} = ',num2str(neuromap_corrs.corr(SigInd),3),', ',pval_format,''],'FontWeight','normal');
    xlabel(xlabel_var)
    ylabel(neuromap_corrs.name{SigInd})
    colormap(turbo)
    set(gca,'FontSize',16)
    s.MarkerEdgeColor = [0 0 0];
    %a = annotation('textbox',[0 .975 .05 .025],'String',annotlabels{i},'FontSize',32,'EdgeColor','none');
    a = annotation('textbox',[0 .905 .05 .13],'String',annotlabels{i},'FontSize',32,'EdgeColor','none');
    
    print([outdir,'/Sigmap_',annotlabels{i},'.png'],'-dpng','-r300')
    close all
    
    SigMapDescrips{i} = neuromap_corrs.description{SigInd};
end