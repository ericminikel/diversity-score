#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_3.0.2/bin/Rscript

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-p", "--pcs"), action="store", default=NA, type='character',
              help="Path to CSV file of principal components"),
  make_option(c("-s", "--samples"), action="store", default=NA, type='character',
              help="Comma-separated list of samples to highlight"),
  make_option(c("-o", "--outpng"), action="store", default='variant.png', type='character',
              help="PNG file to save plot [default %default]"),
  make_option(c("-t", "--title"), action="store", default='Principal components', type='character',
              help="Main title of plot [default %default]"),
  make_option(c("-u", "--subtitle"), action="store", default='', type='character',
              help="Subtitle of plot [default %default]")
)
opt = parse_args(OptionParser(option_list=option_list))

pcs = read.table(opt$pcs,sep=',',row.names=1,header=TRUE)
samples = strsplit(opt$samples,',')[[1]]

png(opt$outpng,width=600,height=600)
# first column of pcs is sample
plot(pcs[,1],pcs[,2],pch='.',col='#CCCCCC',
    xlab='PC1',ylab='PC2',main=opt$title,sub=opt$subtitle)
points(pcs[samples,1],pcs[samples,2],pch=19,col='#000000')
dev.off()

