
Process_matrix <- function(matrix.file,predict_file,outFile=NULL,window.size,statFilter=T) {
  bins <- NULL
  matrix.data <- NULL
  if (is.character(matrix.file)) {
    print("Step 0 : File Read ")
    matdf <- read.table(matrix.file, header = F)
    if (ncol(matdf) - nrow(matdf) == 3) {
      colnames(matdf) <- c("chr", "from.coord", "to.coord")
    } else if (ncol(matdf) - nrow(matdf) == 4) {
      colnames(matdf) <- c("id", "chr", "from.coord", "to.coord")
    } else {
      print("Unknown Type of matrix file")
      return(NULL)
    }
    n_bins <- nrow(matdf)
    bins <- data.frame(
      id = 1:n_bins,
      chr = matdf[, "chr"],
      from.coord = matdf[, "from.coord"],
      to.coord = matdf[, "to.coord"]
    )
    matrix.data <- as.matrix(matdf[, (ncol(matdf) - nrow(matdf) + 1):ncol(matdf)])
    print("Step 0 : Done !!")
    pvalue <- rep(1, n_bins)
    mean.cf <- rep(0, times = n_bins)
    local.ext = rep(-0.5, n_bins)
    data <- read.table(predict_file,header=F)
    index <- data[, 1]
    local.ext[index] <- -1

    print("Step 1: find gap domains")
    gap.idx = Which.Gap.Region2(matrix.data=matrix.data, window.size) 

    proc.regions = Which.process.region(rmv.idx=gap.idx, n_bins=n_bins, min.size=3)

  }
  if(statFilter)
  {

    print("Step 2 :  Filtering of false positive TAD boundaries")
    scale.matrix.data = matrix.data
    for( i in 1:(2*window.size) )
    {
      scale.matrix.data[ seq(1+(n_bins*i), n_bins*n_bins, 1+n_bins) ] = scale( matrix.data[ seq(1+(n_bins*i), n_bins*n_bins, 1+n_bins) ] )
    }
    
    print("-- Compute p-values by Wilcox Ranksum Test")
    for( i in 1:nrow(proc.regions))
    {
      start = proc.regions[i, "start"]
      end = proc.regions[i, "end"]
      pvalue[start:end] <- Get.Pvalue(matrix.data=scale.matrix.data[start:end, start:end], size=window.size, scale=1)
    }
    local.ext[intersect( union(which( local.ext==-1), which(local.ext==-1)), which(pvalue<0.05))] = -2
    local.ext[which(local.ext==-1)] = 0
    local.ext[which(local.ext==-2)] = -1
    now_local <- which(local.ext == -1)  
    print(now_local)
    if (length(now_local) >= 2) {
      adjacent_indices <- now_local[-length(now_local)] + 1 == now_local[-1]  
      adjacent_a <- now_local[adjacent_indices]  
      min_pvalue_indices <- c()
      for (i in adjacent_a) {
        pvalues <- pvalue[c(i, i + 1)]  
        max_pvalue_index <- i + which.max(pvalues) - 1
        now_local <- now_local[now_local != max_pvalue_index]
      }
    }
    a_plus_1 <- now_local + 1  

    data_to_save <- cbind(now_local, a_plus_1) 

    print("Step 2 : Done!")
  } else pvalue = 0
    domains = Convert.Bin.To.Domain.TMP(bins = bins, 
                                  signal.idx = now_local, 
                                  gap.idx = gap.idx , 
                                  pvalues = pvalue, 
                                  pvalue.cut = 0.05)
    bedform = domains[, c("chr", "from.coord", "to.coord", "tag")]
    colnames(bedform) = c("chrom", "chromStart", "chromEnd", "name")

    if( !is.null(outFile) ) {
      print("Writing Files")
      outDomain = paste(outFile, ".domain", sep="")
      write.table( domains, file=outDomain, quote=F, row.names=F, col.names=T, sep="\t")
      outBed = paste(outFile, ".bed", sep="")
      write.table( bedform, file=outBed, quote=F, row.names=F, col.names=F, sep="\t")
      
  }
  print("Done!!")
  return(list(domain=domains,bed=bedform))
}

Get.Diamond.Matrix <- function(mat.data, i, size)
{
  n_bins <- nrow( mat.data )
  if(i == n_bins) return(NA)
  
  lowerbound = max( 1, i-size+1 )
  upperbound = min( i+size, n_bins)
  
  return(mat.data[lowerbound:i, (i+1):upperbound])
}

Which.Gap.Region2 <- function(matrix.data, w)
{
  n_bins = nrow(matrix.data)
  gap = rep(0, n_bins)
  
  for(i in 1:n_bins)
  {
    if( sum( matrix.data[i, max(1, i-w):min(i+w, n_bins)] ) == 0 ) gap[i]=-0.5
  }
  
  idx = which(gap == -0.5)
  return(idx)
}
Which.process.region <- function(rmv.idx, n_bins, min.size=3)
{
  gap.idx = rmv.idx
  
  proc.regions = data.frame(start=numeric(0), end=numeric(0))
  proc.set = setdiff(1:n_bins, gap.idx)
  n_proc.set = length(proc.set)
  
  i=1
  while(i < n_proc.set )
  {
    start = proc.set[i]
    j = i+1
    
    while(j <= n_proc.set)
    {
      if( proc.set[j] - proc.set[j-1] <= 1) j = j + 1
      else {
        proc.regions = rbind(proc.regions, c(start=start, end=proc.set[j-1]) )
        i = j
        break
      }
    }
    
    if(j >= n_proc.set ) {
      proc.regions = rbind(proc.regions, c(start=start, end=proc.set[j-1]) )
      break
    }
  }
  
  colnames(proc.regions) = c("start", "end")
  proc.regions <- proc.regions[ which( abs(proc.regions[,"end"] - proc.regions[, "start"]) >= min.size ), ]
  
  return(proc.regions)
}

Get.Upstream.Triangle <- function(mat.data, i, size)
{
  n_bins = nrow(mat.data)
  
  lower = max(1, i-size)
  tmp.mat = mat.data[lower:i, lower:i]
  return(tmp.mat[upper.tri(tmp.mat, diag=F)])
}

Get.Downstream.Triangle <- function(mat.data, i, size)
{
  n_bins = nrow(mat.data)
  if(i==n_bins) return(NA)
  
  upperbound = min(i+size, n_bins)
  tmp.mat = mat.data[(i+1):upperbound, (i+1):upperbound]
  return( tmp.mat[ upper.tri( tmp.mat, diag=F ) ] )
}


Get.Diamond.Matrix2 <- function(mat.data, i, size){
  n_bins = nrow(mat.data)
  new.mat = matrix(rep(NA, size*size), nrow=size, ncol=size)
  for(k in 1:size)
  {
    if(i-(k-1) >= 1 && i < n_bins)
    {
      lower = min(i+1, n_bins)
      upper = min(i+size, n_bins)
      new.mat[size-(k-1), 1:(upper-lower+1)] = mat.data[i-(k-1), lower:upper]
    }
  } 
  return(new.mat)
}


Get.Pvalue <- function(matrix.data, size, scale=1 )
{
  n_bins = nrow(matrix.data)
  pvalue <- rep(1, n_bins)

  for( i in 1:(n_bins-1) )
  {
    dia = as.vector( Get.Diamond.Matrix2(matrix.data, i, size=size) )
    ups = as.vector( Get.Upstream.Triangle(matrix.data, i, size=size) )
    downs = as.vector( Get.Downstream.Triangle(matrix.data, i, size=size))
    x_values <- dia*scale
    y_values <- c(ups, downs)

    wil.test <- wilcox.test(x = x_values, y = y_values, alternative = "less", exact = FALSE)

    pvalue[i] <- wil.test$p.value
  }
  pvalue[ is.na(pvalue) ] = 1
  return(pvalue)
}

Convert.Bin.To.Domain.TMP <- function(bins, signal.idx, gap.idx, pvalues=NULL, pvalue.cut=NULL)
{
  n_bins = nrow(bins)
  ret = data.frame(chr=character(0), from.id=numeric(0), from.coord=numeric(0), to.id=numeric(0), to.coord=numeric(0), tag=character(0), size=numeric(0))
  levels( x=ret[, "tag"] ) = c("domain", "gap", "boundary")
  
  rmv.idx = setdiff(1:n_bins, gap.idx)
  proc.region = Which.process.region(rmv.idx, n_bins, min.size=0)
  from.coord = bins[proc.region[, "start"], "from.coord"]
  n_procs = nrow(proc.region)
  gap = data.frame(chr=rep( bins[1, "chr"], n_procs), from.id=rep(0, n_procs), from.coord=from.coord, to.id=rep(0, n_procs), to.coord=rep(0, n_procs), tag=rep("gap", n_procs), size=rep(0, n_procs), stringsAsFactors=F)
  
  rmv.idx = union(signal.idx, gap.idx)
  proc.region = Which.process.region(rmv.idx, n_bins, min.size=0)
  n_procs = nrow(proc.region)
  from.coord = bins[proc.region[, "start"], "from.coord"]
  domain = data.frame(chr=rep( bins[1, "chr"], n_procs), from.id=rep(0, n_procs), from.coord=from.coord, to.id=rep(0, n_procs), to.coord=rep(0, n_procs), tag=rep("domain", n_procs), size=rep(0, n_procs), stringsAsFactors=F)
  
  rmv.idx = setdiff(1:n_bins, signal.idx)
  proc.region = as.data.frame( Which.process.region(rmv.idx, n_bins, min.size=1) )
  n_procs = nrow(proc.region)
  if(n_procs>0)
  {
    from.coord = bins[proc.region[, "start"]+1, "from.coord"]  
    boundary = data.frame(chr=rep( bins[1, "chr"], n_procs), from.id=rep(0, n_procs), from.coord=from.coord, to.id=rep(0, n_procs), to.coord=rep(0, n_procs), tag=rep("boundary", n_procs), size=rep(0, n_procs), stringsAsFactors=F)
    ret = rbind(ret, boundary)
  }
  
  ret = rbind(gap, domain)
  ret = ret[order(ret[,3]), ]
  
  ret[, "to.coord"] = c(ret[2:nrow(ret), "from.coord"], bins[n_bins, "to.coord"])
  ret[, "from.id"] = match( ret[, "from.coord"], bins[, "from.coord"] )
  ret[, "to.id"] = match(ret[, "to.coord"], bins[, "to.coord"])
  ret[, "size"] = ret[,"to.coord"]-ret[,"from.coord"]
  
  
  return(ret)
}
