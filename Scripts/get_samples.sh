sdir="/mnt/isilon/projects/ecosystem_biology/prospectomicsCourse2024/imp"
tdir="/mnt/lscratch/users/ohickl/Thesis/run_001/samples"
files2cp=(Analysis/annotation/bakta/annotation.faa.gz Analysis/annotation/bakta/annotation.filt.contig2ID.tsv.gz Binning/MetaWrap/bins.contigs.gz Binning/MetaWrap/bins.stats.gz)
for i in M11-01 M08-01 M08-04 M06-02 M06-01 M05-04 M05-03; do
    mkdir -p "$tdir/$i"
    for e in ${files2cp[@]}; do
    cp ${sdir}/${i}/${e} "$tdir/$i"
    done
done
