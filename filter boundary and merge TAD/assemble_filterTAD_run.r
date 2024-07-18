library("optparse")

option_list = list(

    make_option(c("-i", "--input_matrix"),, type="character", help="Matrix file"),
    make_option(c("-p", "--input_predict"), ,type="character", help="Predict file"),
    make_option(c("-o", "--output"),type="character", help="Output file"),
    make_option(c("-w", "--window_size"),  type="integer", help="Window size")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)



mat_in <- opt$input_matrix
predict_file <- opt$input_predict
domain_out<- opt$output
window <- opt$window_size


source('./assemble_filterTAD.R')
Process_matrix(matrix.file=mat_in,predict_file=predict_file,outFile=domain_out,window.size = window)
