function nd23_gseafig(adir, geneset, rankname, yrange)

% GSEAFIG make nice version of gsea enrichment plot
%
%    GSEAFIG(adir, setname) makes a gsea plot.  The first argument is the
%    directory of a gsea analysis (full path) and the second argument is
%    the name of the gene set to plot.

% structure for analysis directory to get file names in directory

    d = dir(adir);
    files = {d.name}';

% summary path file names for pos and neg enrichment (either .xls or .tsv)

    posfile = [adir char(files(~cellfun('isempty', regexp(files, 'gsea_report_for_na_pos.*[sv]$'))))];
    negfile = [adir char(files(~cellfun('isempty', regexp(files, 'gsea_report_for_na_neg.*[sv]$'))))];

% ranked gene list full path and file name

    rankfile = [adir char(files(~cellfun('isempty', regexp(files, 'ranked_gene_list'))))];

% specified gene set full path and file name

    setfile = [adir char(files(~cellfun('isempty', regexp(files, [geneset '\...[sv]$']))))];

% read in summary files, rank file and specified set file

    P = readtable(posfile,  'FileType', 'text', 'Delimiter', '\t');
    N = readtable(negfile,  'FileType', 'text', 'Delimiter', '\t');
    R = readtable(rankfile, 'FileType', 'text', 'Delimiter', '\t');
    S = readtable(setfile,  'FileType', 'text', 'Delimiter', '\t');

% concatenate summary tables for pos and neg

    if isempty(P)
        E = N;
    elseif isempty(N)
        E = P;
    else
        E = [P; N];
    end

% ranks of specified gene set in ranked list

    hits = S.RANKINGENELIST + 1;

% scores from ranked gene lists

    scores = R.SCORE;
    n = length(scores);

% running enrichment score

    es = gsea(scores, hits);

% identify list in summary table, get summary stats

    k = find(strcmp(E.NAME, geneset));
    esmax = E.ES(k);
    nes = E.NES(k);
    p = E.NOMP_val(k);
    q = E.FDRQ_val(k);
    rmax = E.RANKATMAX(k);

% figure subplot and axis dimensions in pixels

    pw = 280;
    ph = 140;
    hh = 20;
    sp = 5;
    bot = 25;
    left = 65;
    right = 25;
    top = 50;
    fw = left + pw + right;
    fh = bot + hh + sp + ph + top;

% axis positions in normalized coords

    posmain = [left / fw, (bot + hh + sp) / fh, pw / fw, ph / fh];
    poshits = [left / fw, bot / fh, pw / fw, hh / fh];
    posname = [left / fw, (bot + hh + sp + ph) / fh, pw / fw, top / fh];

% preferences

    pref = struct(...
        'xtick', [],...
        'ytick', [],...
        'xcolor', 'w',...
        'ycolor', 'w',...
        'box', 'off');
    fpref = struct(...
        'fontname', 'arial',...
        'fontsize', 14,...
        'horizontalalignment', 'center',...
        'verticalalignment', 'top');
    escolor = [0 .6 0];

% draw figure window

    figure('position', [20 20 fw fh]);

% axes for title & stats

    subplot('position', posname);
    text(.5, .8, geneset,...
        'fontname', 'arial',...
        'fontsize', 14,...
        'horizontalalignment', 'center',...
        'interpreter', 'none');
    if p == 0, pnom = 'p_n_o_m < 0.001'; else, pnom = ['p_n_o_m = ' sprintf('%.3f', p)]; end
    if q == 0, fdrq = 'q < 0.001'; else, fdrq = ['q = ' sprintf('%.3f', q)]; end
%     text(.15, .5, pnom, fpref);
%     text(.5,  .5, fdrq, fpref);
%     text(.85, .5, ['NES = ' sprintf('%.3f', nes)], fpref);
    text(.25, .5, pnom, fpref);
    text(.75, .5, ['NES = ' sprintf('%.3f', nes)], fpref);
    set(gca, pref, 'xlim', [0 1], 'ylim', [0 1]);

% running ES plot

    subplot('position', posmain, 'nextplot', 'add');
    plot([1 n], [0 0], 'k-', 'linewidth', .75);
    plot(1:n, es, '-', 'color', escolor, 'linewidth', 3);
    if sign(esmax) == 1, xm = rmax + 1; else, xm = 1 + n - rmax; end
    plot([1 1] * xm, [0 esmax], '-', 'color', escolor, 'linewidth', .75);
    set(gca,...
        'xlim', [1 n],...
        'xtick', [],...
        'xcolor', 'w',...
        'linewidth', 1,...
        'tickdir', 'out',...
        'fontname', 'arial',...
        'fontsize', 12,...
        'ylim', yrange);
    hy = ylabel('Enrichment Score', 'fontname', 'arial', 'fontsize', 14);
    hy.Position(1) = 1 - (.12 * n);

% hits plot

    subplot('position', poshits);
    plot([hits hits]', [0 1], 'k-', 'linewidth', .5);
    set(gca, pref, 'xlim', [1 n], 'ylim', [0 1]);
    xlabel(['Rank in ' rankname], 'fontname', 'arial', 'fontsize', 14, 'color', 'k');

% white background

    set(gcf, 'inverthardcopy', 'off', 'color', 'w');
    
return

    
