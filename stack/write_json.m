flight=40030;
inst=2;
masklim = false;
mypaths=get_paths(flight);
dx = 1200;
nbins = 25;
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));


for ifield = 4:8
    dt=get_dark_times(flight,inst,ifield);
    load(sprintf('%s/stackdat_%s',loaddir,dt.name),'stackdatall');
    load(sprintf('%s/excessdat_%s',loaddir,dt.name),'excessdatall');

    for im = 1:4
        stackdat = stackdatall(im).stackdat;
        excessdat = excessdatall(im).excessdat;
        data(im).m_min = stackdat.m_min;
        data(im).m_max = stackdat.m_max;
        data(im).r_arr = stackdat.rsub_arr;
        data(im).rfull_arr = stackdat.r_arr;

        data(im).profg = stackdat.all.profcbgsub - stackdat.bg.profcbgsub;
        data(im).profgfull = stackdat.all.profcbg - stackdat.bg.profcbg;
        data(im).profex = excessdat.excess.profcbgsub;
        data(im).profexfull = excessdat.excess.profcbg;
        data(im).cov = excessdat.excov.covcbsub;
        data(im).covfull = excessdat.excov.covcb;

        data(im).profpsfcbfull = excessdat.psf.profcb;
        data(im).profpsfpsfull = excessdat.psf.profps;
        data(im).profpsfcb = excessdat.psf.profcbsub;
        data(im).profpsfps = excessdat.psf.profpssub;
        data(im).profpsfcb100 = excessdat.psf.profcb100;
        data(im).profpsfps100 = excessdat.psf.profps100;

        data(im).covpsfcbfull = excessdat.psfcov.covcb;
        data(im).covpsfpsfull = excessdat.psfcov.covps;
        data(im).covpsfcb = excessdat.psfcov.covcbsub;
        data(im).covpsfps = excessdat.psfcov.covpssub;
        data(im).covpsfcb100 = excessdat.psfcov.covcb100;
        data(im).covpsfps100 = excessdat.psfcov.covps100;
        
        data(im).r_weight = stackdat.all.profhitg;

        hscclusdat = get_hsc_clus_prof(flight, inst, ...
            masklim, data(im).r_weight);
        clus = hscclusdat(4).linefitmodelsub;
        clusfull = hscclusdat(4).linefitmodel;
        data(im).profclus = clus;
        data(im).profclusfull = clusfull;


    end

    jsontxt = jsonencode(data);
    fname = strcat(mypaths.alldat, 'TM',num2str(inst),'/',...
        dt.name, '_datafit.json');
    fid = fopen(fname,'wt');
    fprintf(fid, jsontxt);
    fclose(fid);
end