%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MASTER SCRIPT FOR GENERATING FIGURE PANELS:
%
%     5A, 5B, 5C, 5E, 5F, 5G, 5H, 5I, 5J, 5K, 7A, S5A, S5B
%
% for the manuscript:
%
%     IFNg-producing Tfh cells control the differentiation of lung
%     resident memory B cells.
%
%     Nicole Arroyo-Diaz
%     Holly Bachus
%     Amber Papillion
%     Alex Rosenberg
%     Beatriz Leon-Ruiz
%     Andre Ballesteros-Tato
%
%     Immunity. Accepted (June 6, 2023)
%
% RNA-Seq Data available from GEO: GSE208322
%
%     https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208322
%
% Versions:
%
%     Matlab version 9.14 (2023a)
%     Matlab Statistics and Machine Learning Toolbox 12.5 (2023a)
%
% File dependencies:
%
%     nd23_gseafig.m - renders a pub-quality GSEA enrichment plot
%     gsea.m         - computes data for GSEA enrichment plot
%     usa.m          - red-white-blue color map
%     purplewhite.m  - purple-white color map
%
% Environment
%
%     Retain directory structure from this repo and change the variable
%     ROOT to be the root directory in your environment.
%
% Alex Rosenberg
% University of Alabama at Birmingham
% 10-July-2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % root directory
    ROOT = '/Users/afr/Desktop/ArroyoDiaz2023/';
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% READ IN GENE EXPRESSION DATA INTO STRUCTURE FOR FIGURES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % edgeR stats files for four A vs. B comparisons
    comps = [...
        "d12_Cxcr3_Hi_over_Lo";...
        "d30_Bmem_over_GCB";...
        "41_WT_over_IFNgRko";...
        "41_WT_over_d30_GCB_unpaired"];
    files = ROOT + "gex/" + [...
        "nd23_gex_d12_Cxcr3_Hi_over_Lo.txt";...
        "nd23_gex_d30_Bmem_over_GCB.txt";...
        "nd23_gex_41_WT_over_IFNg.txt";...
        "nd23_gex_41_WT_over_d30_GCB.txt"];
    
    % read in stats and expression levels for four A vs. B comparisons
    S = struct;
    for i = 1:4
        AA = readtable(files(i));
            S(i).comparison = comps(i);
            S(i).statsfile = files(i);
            S(i).rnkfile = regexprep(files(i), {'\.txt', '/gex/'}, {'_GSEA.rnk', '/gsea/'});
            S(i).ipafile = "";
            S(i).table = AA(:, {'gene', 'logfc', 'p','fdr'});
            S(i).table.score = sign(AA.logfc) .* -log10(AA.p);
            S(i).cpm = AA(:, 9:end);
            S(i).cpm.Properties.RowNames = AA.gene2;
            if comps(i) == "d30_Bmem_over_GCB"
                S(i).ipafile = regexprep(files(i), {'\.txt', '/gex/'}, {'_IPA_fdr05_2fold_INPUT.txt', '/ipa/'});
            end 
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CREATE INPUT FILES FOR GSEA AND IPA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % flag to create input files for GSEA or IPA or gene set files
    writeflag = 0;
        
    % if writing output files...
    if writeflag == 1
    
        % write ranked gene files for GSEA analysis input
        for i = 1:4
            temp = S(i).table(:, {'gene', 'score'});
            temp.gene = upper(temp.gene);
            writetable(temp, S(i).rnkfile,...
                'FileType', 'text',...
                'WriteVariableNames', false,...
                'Delimiter', '\t');
        end
    
        % make IPA input file for Bmem over GCB
        for i = 1:4
            if S(i).ipafile ~= ""
                writetable(...
                    S(i).table(S(i).table.fdr < .05 & abs(S(i).table.logfc) > 1, {'gene', 'logfc', 'fdr'}),...
                    S(i).ipafile,...
                    'FileType', 'text',...
                    'WriteVariableNames', false,...
                    'Delimiter', '\t');
            end
        end
    
        % make gene set (for GSEA): WT over IFNgRko fdr < .05 up
        k = [S.comparison]' == "41_WT_over_IFNgRko";
        gup = upper(S(k).table.gene(S(k).table.fdr < .05 & S(k).table.logfc > 0));
        writecell(gup, ROOT + "gsea/nd23_gene_set_WT_OVER_IFNRRKO_FDR05_UP.txt");
    
        % make gene set (for GSEA): Bmem over GCB fdr < .05 and 2-fold up
        k = [S.comparison]' == "d30_Bmem_over_GCB";
        gup = upper(S(k).table.gene(S(k).table.fdr < .05 & S(k).table.logfc > 1));
        writecell(gup, ROOT + "gsea/nd23_gene_set_D30_BMEM_OVER_GCB_FDR05_FC2_UP.txt");
    
    end
    
    % clean up workspace for figures
    clearvars -except S ROOT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PCA PLOT FOR WT AND IFNgR-/- (FIGURE 5A)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    data  = S(3).cpm{:, :}';
    data1 = S(3).cpm{S(3).table.fdr < .05, :}';
    n  = size(data, 2);
    n1 = size(data1, 2);
    [~,  score] = pca(data);
    [~, score1] = pca(data1);
    iwt = find(strncmp(S(3).cpm.Properties.VariableNames, 'WT',  2));
    iko = find(strncmp(S(3).cpm.Properties.VariableNames, 'IFN', 3));
    pw = 250;
    left = 45;
    bot = 45;
    top = 100;
    right = 10;
    fw = left + pw + right;
    fh = bot + pw + top;
    axmain   = [left / fw, bot / fh, pw / fw, pw / fh];
    axlegend = [left / fw, (bot + pw) / fh, top / fw, top / fh];
    figure('position', [20 20 fw fh]);
    % axes for legend
    subplot('position', axlegend, 'nextplot', 'add');
        plot(.20, .37, 'ko', 'markersize', 9, 'linewidth', 1);
        plot(.37, .37, 'mv', 'markersize', 9, 'linewidth', 1);
        plot(.20, .20, 'ko', 'markersize', 9, 'linewidth', 1, 'markerfacecolor', 'k');
        plot(.37, .20, 'mv', 'markersize', 9, 'linewidth', 1, 'markerfacecolor', 'm');
        text(.50, .37, ['All (n = ' num2str(n) ')'], 'fontname', 'arial', 'fontsize', 14);
        text(.50, .20, ['FDR < .05 (n = ' num2str(n1) ')'], 'fontname', 'arial', 'fontsize', 14);
        text(.20, .50, 'WT', 'fontname', 'arial', 'fontsize', 14, 'rotation', 60);
        text(.37, .50, 'IFNgR^-^/^-', 'fontname', 'arial', 'fontsize', 14, 'rotation', 60);
        set(gca,...
            'xlim', [0 1],...
            'ylim', [0 1],...
            'xcolor', 'w',...
            'ycolor', 'w',...
            'xtick', [],...
            'ytick', [],...
            'box', 'off');
    % axes for main plot
    subplot('position', axmain, 'nextplot', 'add');
        plot([-24 24], [0 0], '-', 'color', [1 1 1] * .6, 'linewidth', .5);
        plot([0 0], [-24 24], '-', 'color', [1 1 1] * .6, 'linewidth', .5);
        plot(score(iwt,  1), score(iwt,  2), 'ko', 'markersize', 9, 'linewidth', 1);
        plot(score(iko,  1), score(iko,  2), 'mv', 'markersize', 9, 'linewidth', 1);
        plot(score1(iwt, 1), score1(iwt, 2), 'ko', 'markersize', 9, 'linewidth', 1, 'markerfacecolor', 'k');
        plot(score1(iko, 1), score1(iko, 2), 'mv', 'markersize', 9, 'linewidth', 1, 'markerfacecolor', 'm');
        set(gca,...
            'linewidth', 1,...
            'box', 'on',...
            'ticklength', [0 0],...
            'fontname', 'arial',...
            'fontsize', 12,...
            'xlim', [-24 24],...
            'ylim', [-24 24],...
            'xtick', -20:5:20,...
            'ytick', -20:5:20);
        xlabel('PC1', 'fontname', 'arial', 'fontsize', 14);
        ylabel('PC2', 'fontname', 'arial', 'fontsize', 14);
    set(gcf, 'inverthardcopy', 'off', 'color', 'w');     
    clearvars -except S ROOT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% VOLCANO PLOT FOR WT AND IFNgR-/- (FIGURE 5C, 7A)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
    x = S(3).table.logfc;
    y = -log10(S(3).table.p);
    g = S(3).table.gene;
    q05 = -log10(max(S(3).table.p(S(3).table.fdr < .05)));
    isig = abs(x) > 1 & S(3).table.fdr < .05;
    igoi = [find(strcmp(g, 'Cxcr3')); find(strcmp(g, 'Tbx21'))];
    pw = 350;
    ph = 300;
    bot = 50;
    left = 50;
    top = 10;
    right = 15;
    fw = left + pw + right;
    fh = bot + ph + top;
    axmain = [left / fw, bot / fh, pw / fw, ph / fh];
    figure('position', [20 20 fw fh]);
    subplot('position', axmain, 'nextplot', 'add');
    plot([1 1], [0 14], 'r-', 'linewidth', .75);
    plot([-1 -1], [0 14], 'r-', 'linewidth', .75);
    plot([-7 7], [1 1] * q05, 'r--', 'linewidth', .75);
    plot(x(~isig), y(~isig), '.', 'markersize', 9, 'color', [1 1 1] * .65);
    plot(x(isig), y(isig), '.', 'markersize', 9, 'color', [0 0 1]);
    plot(x(igoi), y(igoi), 'ro', 'markersize', 10, 'linewidth', 1);
    text(x(igoi) + .3, y(igoi), g(igoi), 'fontname', 'arial', 'fontsize', 11, 'color', 'r');
    set(gca,...
        'ticklength', [0 0],...
        'linewidth', 1,...
        'xlim', [-7 7],...
        'ylim', [0 14],...
        'fontname', 'arial',...
        'fontsize', 12,...
        'box', 'on');
    xlabel('log_2(WT / IFNgR^-^/^-)', 'fontname', 'arial', 'fontsize', 14);
    ylabel('-log_1_0(p-Value)', 'fontname', 'arial', 'fontsize', 14);
    set(gcf, 'inverthardcopy', 'off', 'color', 'w');     
    clearvars -except S ROOT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SELECTED GENE HEATMAP (FIGURE S5A)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    goi = ["Zbtb32"; "Klf2"; "Hhex"; "Ski"; "Bcl2"; "Aicda"; "Bcl6"];
    m = S(2).cpm(goi, :);
    z = linkage(pdist(m{:, :}, 'euclidean'), 'average');
    f = figure; [~, ~, isort] = dendrogram(z, 0); close(f)
    m1 = m(isort, :);
    g1 = m1.Properties.RowNames;
    s1 = regexprep(m1.Properties.VariableNames, '_.*', '');
    data = zscore(m1{:, :}, [], 2);
    zrange = [-2 2];
    cmap = usa(256);
    [ng, ns] = size(data);
    pix = 25;
    pw = pix * ns;
    ph = pix * ng;
    left = 10;
    bot = 20;
    right = 60;
    top = 50;
    sp = 40;
    cw = 15;
    ch = 120;
    border = 25;
    fw = left + pw + right + sp + cw + border;
    fh = bot + ph + top;
    axmain = [left / fw, bot / fh, pw / fw, ph / fh];
    axgene = [(left + pw) / fw, bot / fh, right / fw, ph / fh];
    axsamp = [left / fw, (bot + ph) / fh, pw / fw, top / fh];
    axcbar = [(left + pw + right + sp) / fw, bot / fh, cw / fw, ch / fh];
    pref = struct(...
        'xcolor', 'w',...
        'ycolor', 'w',...
        'ticklength', [0 0],...
        'xtick', [],...
        'ytick', [],...
        'box', 'off');        
    figure('position', [20 20 fw fh]);   
    % axes for sample labels
    subplot('position', axsamp);
        text(1:ns, ones(ns, 1) * .15, s1,...
            'fontname', 'arial',...
            'fontsize', 14,...
            'rotation', 45);
        set(gca, pref,...
            'xlim', [.5 ns + .5],...
            'ylim', [0 1]);  
    % axes for gene names
    subplot('position', axgene);
        text(ones(ng, 1) * .05, 1:ng, g1,...
            'fontname', 'arial',...
            'fontsize', 14,...
            'fontangle', 'italic');
        set(gca, pref,...
            'ydir', 'reverse',...
            'xlim', [0 1],...
            'ylim', [.5 ng + .5]);        
    % axes for main heatmap
    subplot('position', axmain, 'nextplot', 'add');
        imagesc(data, zrange);
        colormap(cmap);
        set(gca,...
            'ydir', 'reverse',...
            'xlim', [.5 ns + .5],...
            'ylim', [.5 ng + .5],...
            'ticklength', [0 0],...
            'xtick', [],...
            'ytick', []);
        plot(repmat((0:ns) + .5, 2, 1), [.5 ng + .5],...
            '-', 'color',...
            [1 1 1] * .5,...
            'linewidth', .5);
        plot([.5 ns + .5], repmat((0:ng) + .5, 2, 1),...
            '-', 'color',...
            [1 1 1] * .5,...
            'linewidth', .5);
        cbar = colorbar;
        cbar.Position = axcbar;
        set(cbar,...
            'yaxislocation', 'left',...
            'ticklength', 0,...
            'box', 'on',...
            'fontname', 'arial',...
            'fontsize', 12);
        cbar.Label.String = 'z-Score';
        cbar.Label.FontSize = 14;
        cbar.Label.FontName = 'arial';
        cbar.Label.Position = [0, (diff(zrange) * .1) + max(zrange), 0];
        cbar.Label.Rotation = 0;
    set(gcf, 'inverthardcopy', 'off', 'color', 'w');            
    clearvars -except S ROOT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% IPA UPSTREAM REGULATORS: BMEM OVER GCB FDR < .05, 2-FOLD (FIGURE 5H)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    IPA = readtable([ROOT '/ipa/nd23_gex_d30_Bmem_over_GCB_IPA_fdr05_2fold_REGULATORS.txt']);
    IPA = IPA(:, [1 2 5 7]);
    IPA.Properties.VariableNames = {'reg', 'logfc', 'z', 'bhp'};
    IPA.neglogbhp = -log10(IPA.bhp);
    IPA = IPA(IPA.bhp < .01 & abs(IPA.z) > 3, :);
    IPA = sortrows(IPA, 'z', 'descend');
    for i = 1:size(IPA, 1)
        temp = lower(IPA.reg{i});
        temp(1) = upper(temp(1));
        IPA.reg(i) = cellstr(temp);
    end
    ncol = 256;
    cmap = flipud(purplewhite(ncol));
    left = 15;
    bot = 50;
    top = 10;
    right = 40;
    cw = 15;
    ch = 100;
    border = 30;
    pix = 18;
    pw = 150;
    lw = 60;
    ph = size(IPA, 1) * pix;
    fw = left + pw + lw + pw + right + cw + border;
    fh = bot + ph + top;
    axgene = [(left + pw) / fw, bot / fh, lw / fw, ph / fh];
    axpos  = [(left + pw + lw) / fw, bot / fh, pw / fw, ph / fh];
    axneg  = [left / fw, bot / fh, pw / fw, ph / fh];
    axcbar = [(left + pw + lw + pw + right) / fw, bot / fh, cw / fw, ch / fh];
    pref = struct(...
        'ticklength', [0 0],...
        'xcolor', 'w',...
        'ycolor', 'w',...
        'xtick', [],...
        'ytick', [],...
        'box', 'off');
    figure('position', [20 20 fw fh]);
    % axes for upstream regulator names
    subplot('position', axgene);
        text(ones(size(IPA, 1), 1) * .5, 1:size(IPA, 1), IPA.reg,...
            'fontname', 'arial',...
            'fontsize', 14,...
            'fontangle', 'italic',...
            'horizontalalignment', 'center');
        set(gca, pref,...
            'xlim', [0 1],...
            'ylim', [.5 size(IPA, 1) + .5],...
            'ydir', 'reverse');
    % axes for bar plot of negative z-score regulators
    an = subplot('position', axneg, 'nextplot', 'add');
    % axes for bar plot of positive z-score regulators
    ap = subplot('position', axpos, 'nextplot', 'add');
        for i = 1:size(IPA, 1)
            v = round((IPA.neglogbhp(i) / 10) * (ncol - 1)) + 1;
            if v > ncol, v = ncol; end
            if IPA.z(i) > 0, ax = ap; else, ax = an; end
            barh(ax, i, IPA.z(i), .7, 'facecolor', cmap(v, :));
        end
        set([ap an],...
            'ylim', [.5 size(IPA, 1) + .5],...
            'ytick', [],...
            'ycolor', 'w',...
            'linewidth', 1,...
            'tickdir', 'out',...
            'fontname', 'arial',...
            'fontsize', 12,...
            'ydir', 'reverse');
        set(an, 'yaxislocation', 'right', 'xlim', [-10 0]);
        set(ap, 'yaxislocation', 'right', 'xlim', [0 10]);
        xlabel([an ap], 'Activation z-Score', 'fontname', 'arial', 'fontsize', 14);
    % axes for colorbar
    subplot('position', axcbar);
        imagesc((0:.1:10)');
        colormap(cmap);
        set(gca,...
            'ydir', 'normal',...
            'ticklength', [0 0],...
            'ytick', [min(ylim) mean(ylim) max(ylim)],...
            'yticklabel', {'0', '5', '\geq10'},...
            'xtick', 0,...
            'box', 'on',...
            'linewidth', .5,...
            'fontname', 'arial',...
            'fontsize', 12);
        text(1, 125, '-log_1_0 BHP',...
            'fontname', 'arial',...
            'fontsize', 14,...
            'horizontalalignment', 'center');
    set(gcf, 'inverthardcopy', 'off', 'color', 'w');            
    clearvars -except S ROOT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GSEA SIGNIFICANT HALLMARKS BUBBLE PLOT - BMEM OVER GCB (FIGURE 5I)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    gfile = ROOT + "gsea/nd23_gsea_d30_Bmem_over_GCB_vs_HALLMARK_1656454322517/" + ...
        "gsea_report_for_na_pos_1656454322517.tsv";    
    G0 = readtable(gfile, 'FileType', 'text', 'Delimiter', '\t');
    G = G0(:, {'NAME', 'NES', 'NOMP_val', 'FDRQ_val'});
    G = G(G.FDRQ_val < .05, :);
    G.FDRQ_val(G.FDRQ_val == 0) = .0001;    
    nh = size(G, 1);
    pw = 180;
    ph = 280 * nh / 11; % to keep height per row same as other bubble plot
    left = 10;
    bot = 50;
    top = 10;
    right = 350;
    fw = left + pw + right;
    fh = bot + ph + top;
    axgene = [(left + pw) / fw, bot / fh, right / fw, ph / fh];
    axmain = [left / fw, bot / fh, pw / fw, ph / fh];
    figure('position', [20 20 fw fh]);
    % axes for gene names
    subplot('position', axgene, 'nextplot', 'add');
        text(.015 * ones(nh, 1), 1:nh, G.NAME,...
            'fontname', 'arial',...
            'fontsize', 12,...
            'interpreter', 'none');
        set(gca,...
            'xcolor', 'w',...
            'ycolor', 'w',...
            'xtick', [],...
            'ytick', [],...
            'ticklength', [0 0],...
            'box', 'off',...
            'xlim', [0 1],...
            'ylim', [0 nh + 1],...
            'ydir', 'reverse');
    % axes for main plot
    subplot('position', axmain, 'nextplot', 'add');
        for i = 1:nh
            msize = 5 + -log10(G.FDRQ_val(i)) * 6;
            plot([0 3], [1 1] * i, '-', 'color', [1 1 1] * .6, 'linewidth', .5);
            if G.FDRQ_val(i) < .05, col = 'b'; else, col = [1 1 1] * .5; end
            plot(G.NES(i), i, 'o', 'color', col, 'markerfacecolor', col, 'markersize', msize);
        end
        set(gca,...
            'tickdir', 'out',...
            'ytick', [],...
            'ycolor', 'w',...
            'box', 'off',...
            'linewidth', 1,...
            'ylim', [0 nh + 1],...
            'fontname', 'arial',...
            'fontsize', 12,...
            'ydir', 'reverse');
        xlabel('Normalized Enrichment', 'fontname', 'arial', 'fontsize', 14);
    set(gcf, 'inverthardcopy', 'off', 'color', 'w'); 
    clearvars -except S ROOT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GSEA IMMUNE SIG C7 (SELECTED) BUBBLE PLOT - BMEM OVER GCB (FIGURE S5B)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    imm = [...
        "GSE11961_GERMINAL_CENTER_BCELL_DAY7_VS_MEMORY_BCELL_DAY40_DN", "GC",    "GSE11961";...
        "GSE11961_MEMORY_BCELL_DAY7_VS_GERMINAL_CENTER_BCELL_DAY7_UP",  "GC",    "GSE11961";...
        "GSE12366_GC_VS_MEMORY_BCELL_DN",                               "GC",    "GSE12366";...
        "GSE12366_NAIVE_VS_MEMORY_BCELL_DN",                            "Naive", "GSE12366";...
        "GSE13411_NAIVE_VS_SWITCHED_MEMORY_BCELL_DN",                   "Naive", "GSE13411";...
        "GSE12366_PLASMA_CELL_VS_MEMORY_BCELL_DN",                      "PC",    "GSE12366";...
        "GSE13411_SWITCHED_MEMORY_BCELL_VS_PLASMA_CELL_UP",             "PC",    "GSE13411";...
        "GSE13411_PLASMA_CELL_VS_MEMORY_BCELL_DN",                      "PC",    "GSE13411";...
        "GSE11961_MEMORY_BCELL_DAY7_VS_PLASMA_CELL_DAY7_UP",            "PC",    "GSE11961"];
    nimm = size(imm, 1);
    gfile = ROOT + "gsea/nd23_gsea_d30_Bmem_over_GCB_vs_IMMSIG_C7_1656450499616_REPORT_positive.tsv";
    G0 = readtable(gfile, 'FileType', 'text');
    k = nan(nimm, 1);
    for i = 1:nimm, k(i) = find(G0.NAME == imm(i, 1)); end
    G = G0(k, {'NAME', 'NES', 'NOMP_val', 'FDRQ_val'});
    G.FDRQ_val(G.FDRQ_val == 0) = .0001;
    y = [1 2 3 5 6 8 9 10 11];
    pw = 180;
    ph = 280;
    left = 10;
    bot = 50;
    top = 10;
    right = 500;
    fw = left + pw + right;
    fh = bot + ph + top;
    axlegend = [(left + pw) / fw, 0, right / fw, bot / fh];
    axgeneset = [(left + pw) / fw, bot / fh, right / fw, ph / fh];
    axmain = [left / fw, bot / fh, pw / fw, ph / fh];
    figure('position', [20 20 fw fh]);
    % axes for gene set (GSE) number
    subplot('position', axlegend, 'nextplot', 'add');
        text(.1,  .5, 'GSE11961', 'fontname', 'arial', 'fontsize', 12, 'color', [0.0 0.6 0.0]);
        text(.25, .5, 'GSE12366', 'fontname', 'arial', 'fontsize', 12, 'color', [0.6 0.0 0.6]);
        text(.4,  .5, 'GSE13411', 'fontname', 'arial', 'fontsize', 12, 'color', [0.8 0.4 0.0]);
        set(gca,...
            'xcolor', 'w',...
            'ycolor', 'w',...
            'xtick', [],...
            'ytick', [],...
            'ticklength', [0 0],...
            'box', 'off',...
            'xlim', [0 1],...
            'ylim', [0 1]);
    % axes for gene set name and set grouping labels
    subplot('position', axgeneset, 'nextplot', 'add');
        for i = 1:nimm
            switch imm(i, 3)
                case "GSE11961", col = [0.0 0.6 0.0];
                case "GSE12366", col = [0.6 0.0 0.6];
                case "GSE13411", col = [0.8 0.4 0.0];
            end        
            text(.015, y(i), regexprep(imm(i, 1), 'GSE....._', ''),...
                'color', col,...
                'fontname', 'arial',...
                'fontsize', 12,...
                'interpreter', 'none');
        end
        plot(.83 * [1 1], [.5 3.5],   '-', 'color', [1 1 1] * .4, 'linewidth', 2);
        plot(.83 * [1 1], [4.5 6.5],  '-', 'color', [1 1 1] * .4, 'linewidth', 2);
        plot(.83 * [1 1], [7.5 11.5], '-', 'color', [1 1 1] * .4, 'linewidth', 2);
        text(.845, 2,   {'Up in Mem'; 'vs GC'},    'fontname', 'arial', 'fontsize', 14);
        text(.845, 5.5, {'Up in Mem'; 'vs Naive'}, 'fontname', 'arial', 'fontsize', 14);
        text(.845, 9.5, {'Up in Mem'; 'vs PC'},    'fontname', 'arial', 'fontsize', 14);
        set(gca,...
            'xcolor', 'w',...
            'ycolor', 'w',...
            'xtick', [],...
            'ytick', [],...
            'ticklength', [0 0],...
            'box', 'off',...
            'xlim', [0 1],...
            'ylim', [0 max(y) + 1],...
            'ydir', 'reverse');
    % axes to plot data
    subplot('position', axmain, 'nextplot', 'add');
        for i = 1:nimm
            msize = 5 + -log10(G.FDRQ_val(i)) * 6;
            plot([0 3], [1 1] * y(i), '-', 'color', [1 1 1] * .6, 'linewidth', .5);
            if G.FDRQ_val(i) < .05, col = 'b'; else, col = [1 1 1] * .5; end
            plot(G.NES(i), y(i), 'o', 'color', col, 'markerfacecolor', col, 'markersize', msize);
        end
        set(gca,...
            'tickdir', 'out',...
            'ytick', [],...
            'ycolor', 'w',...
            'box', 'off',...
            'linewidth', 1,...
            'ylim', [0 max(y) + 1],...
            'fontname', 'arial',...
            'fontsize', 12,...
            'ydir', 'reverse');
        xlabel('Normalized Enrichment', 'fontname', 'arial', 'fontsize', 14);
    set(gcf, 'inverthardcopy', 'off', 'color', 'w');
    clearvars -except S ROOT
 
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INDIVIDUAL GSEA ENRICHMENT PLOTS (FIGURE 5B, 5E, 5F, 5G, 5J, 5K)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % FIGURE 5B: ranked list = WT over IFNgR-/-, set = Mem B cell sig
    nd23_gseafig(...
        [ROOT 'gsea/nd23_gsea_41_WT_over_IFNg_vs_Custom_Sets_1656533192856/'],...
        'D30_BMEM_OVER_GCB_FDR05_FC2_UP',...
        'WT over IFNgR^-^/^-', ...
        [-.25 .6]);

    % FIGURE 5E: ranked list = Cxcr3 hi over lo, set = IFNg-induced prog
    nd23_gseafig(...
        [ROOT 'gsea/nd23_gsea_d12_Cxcr3_Hi_over_Lo_vs_Custom_Sets_1657226055952/'],...
        'WT_OVER_IFNGRKO_FDR05_UP',...
        'Cxcr3 Hi over Cxcr3 Lo',...
        [-.1 .8]);

    % FIGURE 5F: ranked list = Cxcr3 hi over lo, set = Mem B cell sig
    nd23_gseafig(...
        [ROOT 'gsea/nd23_gsea_d12_Cxcr3_Hi_over_Lo_vs_Custom_Sets_1657226055952/'],...
        'D30_BMEM_OVER_GCB_FDR05_FC2_UP',...
        'Cxcr3 Hi over Cxcr3 Lo',...
        [-.1 .8]);

    % FIGURE 5G: ranked list = Cxcr3 hi over lo, set = Tan et al 2022
    nd23_gseafig(...
        [ROOT 'gsea/nd23_gsea_d12_Cxcr3_Hi_over_Lo_vs_TAN_1682547648588/'],...
        'TAN_2022_LUNG_BRM_GREATER_THAN_AVERAGE',...
        'Cxcr3 Hi over Cxcr3 Lo',...
        [-.2 .6]);

    % FIGURE 5J: ranked list = Bmem over GCB, set = Hallmark IFNg
    nd23_gseafig(...
        [ROOT 'gsea/nd23_gsea_d30_Bmem_over_GCB_vs_HALLMARK_1656454322517/'],...
        'HALLMARK_INTERFERON_GAMMA_RESPONSE',...
        'D30 Bmem over D30 GCB',...
        [-.1 .6]);

    % FIGURE 5K: ranked list = Bmem over GCB, set = IFNg-induced prog
    nd23_gseafig(...
        [ROOT 'gsea/nd23_gsea_d30_Bmem_over_GCB_vs_41_WT_over_IFNgRko.1657225051104/'],...
        'WT_OVER_IFNGRKO_FDR05_UP',...
        'D30 Bmem over D30 GCB',...
        [-.1 .6]);


