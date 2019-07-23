library(RJSONIO)

args = commandArgs(trailingOnly=TRUE)

indir = paste0(args[2], "_in")
outdir = args[2]

dir.create(indir)
dir.create(outdir)

unzip(args[1], exdir=indir) 
url <- paste0('http://proteomics2.ucsd.edu/ProteoSAFe/DownloadResultFile?task=', args[1], '&block=m=main&file=req/')
download.file(url, paste0(indir, "/req.zip"))
unzip(paste0(indir, "/req.zip"), exdir=indir)

url <- paste0('http://proteomics2.ucsd.edu/ProteoSAFe/DownloadResultFile?task=', args[1], '&block=m=main&file=final_out/node_attributes_table.tsv')
download.file(url, paste0(outdir, "/node_attributes_table.tsv"))

load(paste0(indir, "/split_data/tabgnps.rda"))
write.table(tabgnps1, paste0(outdir, "/tabgnps.tsv"), row.names=FALSE, sep="\t", quote=FALSE)

load(paste0(indir, "/split_data/net.rda"))
write.table(net, paste0(outdir, "/net.tsv"), row.names=FALSE, sep="\t", quote=FALSE)

load(paste0(indir, "/split_data/allspectra.rda"))
write.table(allspectra2, paste0(outdir, "/allspectra.mgf"), row.names=FALSE, col.names=FALSE, quote=FALSE)

url <- paste0('http://proteomics2.ucsd.edu/ProteoSAFe/DownloadResultFile?task=', args[1], '&block=m=main&file=lid_res/lid.rda')
download.file(url, paste0(indir, "/lid.rda"))

load(paste0(indir, "/lid.rda"))
exportJson <- RJSONIO::toJSON(lid, pretty = TRUE, auto_unbox = TRUE)
write(exportJson, paste0(outdir, "/lid.json"))

load(paste0(indir, "/merge_fusion/fusion.rda"))
targJson <- RJSONIO::toJSON(targL, pretty = TRUE, auto_unbox = TRUE)
write(targJson, paste0(outdir, "/fusion.json"))

load(paste0(indir, "/merge_consensus/consensus.rda"))
targ2Json <- RJSONIO::toJSON(targL2, pretty = TRUE, auto_unbox = TRUE)
write(targ2Json, paste0(outdir, "/consensus.json"))

url <- paste0('http://proteomics2.ucsd.edu/ProteoSAFe/DownloadResultFile?task=', args[1], '&block=m=main&file=metfrag_out/mlist.rda')
download.file(url, paste0(indir, "/mlist.rda"))

load(paste0(indir, "/mlist.rda"))
mlistJson <- RJSONIO::toJSON(mlist, pretty = TRUE, auto_unbox = TRUE)
write(mlistJson, paste0(outdir, "/mlist.json"))

unlink(indir, recursive=TRUE)
