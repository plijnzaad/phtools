
read.binding.site.list <- function(binding.site.list){
    if( sum(sapply(binding.site.list, function(x) "weight" %in% colnames(x))) != length(binding.site.list)){
        #some binding site don't have weight, set to 1
        binding.site.list <- lapply(binding.site.list, function(x) {x$weight<-1; return(x)})
    }
    return(lapply(binding.site.list, function(x) split(x[,setdiff(colnames(x), "chr")], x$chr, drop = TRUE) ) )
}

#return union of all chromosomes in all binding site lists, i.e, all chromosomes that are needed
get.binding.chr <- function(bs.list){
    chr.list <- lapply(bs.list, function(x) names(x))
    chr.binding.vec <- unique(unlist(chr.list, use.names=FALSE))
    return(chr.binding.vec)
}

#return chromosomes in ChIP and/or control samples
get.chr <- function(chip.list, input.list, chr.exclusion=NULL){
    if(!is.null(chr.exclusion)){
        temp <- lapply(chip.list, function(x) setdiff(names(x), chr.exclusion))
    }else{
        temp <- lapply(chip.list, function(x) names(x))
    }
    len <- sapply(temp, length)
    chr.chip.unique <- unique(unlist(temp, use.names=FALSE))
    len.unique <- length(chr.chip.unique)
    if(sum(len==len.unique)!=length(len)){
        cat("ChIP samples have different set of chromosome.\n")
        cat("Please check the data.\n")
        cat("One option is to filter the chromosomes that are not shared by all samples, possibly through chr.exclusion parameter.\n")
        cat("Existing chromosomes for each ChIP sample are:\n")
        print(temp)
        stop("ChIP samples have different set of chromosomes.")
    }

    if(!is.null(input.list)){
        if(!is.null(chr.exclusion)){
            temp <- lapply(input.list, function(x) setdiff(names(x), chr.exclusion))
        }else{
            temp <- lapply(input.list, function(x) names(x))
        }
        len <- sapply(temp, length)
        chr.input.unique <- unique(unlist(temp, use.names=FALSE))
        len.unique <- length(chr.input.unique)
        if(sum(len==len.unique)!=length(len)){
            cat("Control samples have different set of chromosome.\n")
            cat("Please check the data.\n")
            cat("One option is to filter the chromosomes that are not shared by all samples, possibly through chr.exclusion parameter.\n")
            cat("Existing chromosomes for each control sample are:\n")
            print(temp)
            stop("Control samples have different set of chromosomes.")
        }
        if(!setequal(chr.chip.unique, chr.input.unique)){
            cat("ChIP and control samples have different set of chromosomes.\n")
            cat("ChIP sample chromosomes are\n")
            print(chr.chip.unique)
            cat("Control sample chromosomes are\n")
            print(chr.input.unique)
            stop("ChIP and control samples have different set of chromosomes.")
        }
    }

    return(chr.chip.unique)
}


get.chr.len <- function(chip.list, input.list=NULL, chr.vec){
    chr.len.vec <- sapply( chr.vec, function(x) max(sapply( chip.list, 
        function(y) max(max(y[[x]][["+"]]), max(y[[x]][["-"]])) )) )
    if(!is.null(input.list)){
        chr.len.vec.input <- sapply( chr.vec, function(x) max(sapply( input.list, 
            function(y) max(max(y[[x]][["+"]]), max(y[[x]][["-"]])) )) )
        chr.len.vec <- pmax(chr.len.vec, chr.len.vec.input)
    }
    return(chr.len.vec)
}

get.NCIS.norm.factor <- function(dat){
    norm.factor.vec <- sapply(names(dat$chip.list), function(x) NCIS.internal(
        dat$chip.list[[x]], dat$input.list[[dat$matching.input.names[x]]], 
        chr.vec=dat$chr.vec, chr.len.vec=dat$chr.len.vec)$est)
    return(norm.factor.vec)
}

DBChIP <- function(binding.site.list, chip.data.list, conds, input.data.list=NULL, 
    data.type=c("MCS", "AlignedRead", "BED"), frag.len=200,
    chr.vec=NULL, chr.exclusion=NULL,
    chr.len.vec=NULL, subtract.input=FALSE, norm.factor.vec=NULL,
    in.distance=100, out.distance=250, window.size=250,
    dispersion=NULL, common.disp=TRUE, prior.n=10,
    two.sample.method="composite.null", allowable.FC=1.5, collapsed.quant=0.5){

    bs.list <- read.binding.site.list(binding.site.list)

    #if no binding site found in one condition, which is unlikely, skip the next command
    #stopifnot(setequal(names(bs.list), unique(as.character(conds))))

    ## compute consensus site
    consensus.site <- site.merge(bs.list, in.distance=in.distance, out.distance=out.distance)
    
    dat <- load.data(chip.data.list=chip.data.list, conds=conds, consensus.site=consensus.site, 
        input.data.list=input.data.list, data.type=data.type, chr.vec=chr.vec, chr.exclusion=chr.exclusion, chr.len.vec=chr.len.vec, 
        norm.factor.vec=norm.factor.vec, frag.len=frag.len)
    
    dat <- get.site.count(dat, subtract.input=subtract.input, window.size=window.size)

    dat <- test.diff.binding(dat, dispersion=dispersion, common.disp=common.disp, prior.n=prior.n, 
        two.sample.method=two.sample.method, allowable.FC=allowable.FC, collapsed.quant=collapsed.quant)

    return(dat)
}

median.ratio <- function (counts) {
    #lib.size <- colSums(counts)
    gm <- exp(rowMeans(log(counts)))
    apply(counts, 2, function(x) median((x/gm)[gm > 0]))
}

#compute sample background read depth by excluding binding reads
#chop genome into 500bp bins, exclude the bin the sites in and their left and right neighbors.
comp.background.size <- function(dat, consensus.site, binsize=500, shift.size=100){
    return(comp.bkg.size(dat$chip.list, consensus.site, dat$chr.vec, dat$chr.len.vec, dat$input.list, binsize=binsize, shift.size=shift.size))
}

comp.bkg.size <- function(chip.list, consensus.site, chr.vec, chr.len.vec, input.list=NULL, binsize=500, shift.size=100){
    res <- NULL
    res.input <- NULL
    if(class(consensus.site)=="data.frame"){
        consensus.site <- split(consensus.site[,setdiff(colnames(consensus.site), "chr")], consensus.site$chr, drop = TRUE)
    }
    for(i in 1:length(chr.vec)){
        chr <- chr.vec[i]
        bk.f <- c(0, seq(from=binsize-shift.size, by=binsize, length.out=ceiling(chr.len.vec[i]/binsize)-1), chr.len.vec[i])
        bk.r <- c(0, seq(from=binsize+shift.size, by=binsize, length.out=ceiling(chr.len.vec[i]/binsize)))
        bin.count <- sapply(chip.list, function(x) hist(x[[chr]][["+"]], breaks=bk.f, plot = FALSE)$counts) + sapply(chip.list, function(x) hist(x[[chr]][["-"]], breaks=bk.r, plot = FALSE)$counts)
        #bin.count: #bins * #chip samples
        if(!is.null(input.list)) bin.count.input <- sapply(input.list, function(x) hist(x[[chr]][["+"]], breaks=bk.f, plot = FALSE)$counts) + sapply(input.list, function(x) hist(x[[chr]][["-"]], breaks=bk.r, plot = FALSE)$counts)
        if(!is.null(consensus.site[[chr]])){
            index <- ceiling(consensus.site[[chr]]$pos/binsize)
            index.all <- cbind(index, index+1, index-1)
            index.all <- index.all[index.all>0 & index.all<=nrow(bin.count)]
            index.all <- sort(unique(index.all))
            bin.count <- bin.count[-index.all, ]
            if(!is.null(input.list)) bin.count.input <- bin.count.input[-index.all, ]
        }
        #return bin*sample matrix
        res <- rbind(res, colSums(bin.count))
        if(!is.null(input.list)) res.input <- rbind(res.input, colSums(bin.count.input))
        cat(".")
    }
    if(!is.null(input.list)){
        return(list(chip.background.size=colSums(res), input.background.size=colSums(res.input)))
    }else{
        return(list(chip.background.size=colSums(res), input.background.size=NULL))
    }
}

read.data.list <- function(data.list, data.type){
    if(data.type=="MCS"){
        res <- lapply(data.list, read.MCS)
    }else{
        if(data.type=="AlignedRead"){
            require(ShortRead)
            res <- lapply(data.list, read.AlignedRead)
        }else{
            if(data.type=="BED"){
                res <- lapply(data.list, read.BED)
            }else{
                stop("unknown format")
            }
        }
    }
    return(res)
}

#cluster close sites together for further test
#input: binding site list for each rep and chr, may come with weight
#output: consensus site locations and their corresponding ChIP counts for each replicate and p-value
site.merge <- function(bs.list, in.distance=100, out.distance=250){
    
    cat("merging sites from different conditions to consensus sites")
    #for each chr, weighted clustering using centroid method
    cond.vec <- names(bs.list)
    consensus.site <- list()
    chr.binding.vec <- get.binding.chr(bs.list)
    for(chr in chr.binding.vec){
        bs.pos <- lapply(cond.vec, function(x) bs.list[[x]][[chr]]$pos)
        if(sum(sapply(bs.pos, length))==0) next
        cond.origin <- rep(cond.vec, times=sapply(bs.pos, length))
        bs.pos <- unlist(bs.pos, use.names = FALSE)
        weight <- unlist(lapply(cond.vec, function(x) bs.list[[x]][[chr]]$weight), use.names = FALSE)
        bs <- data.frame(pos=bs.pos, cond=cond.origin, weight=weight)
        
        od <- order(bs$pos, decreasing = FALSE)
        bs <- bs[od, ]
        
        diff.pos <- diff(bs$pos)
        bk.pts <- which(diff.pos >= out.distance)
        cluster.size <- c(bk.pts, length(bs$pos)) - c(0, bk.pts)
        cum.size <- cumsum(cluster.size)
        res <- data.frame()
        #singleton
        if(sum(cluster.size==1)>0){
            index <- cum.size[cluster.size==1]
            res <- data.frame(pos=bs$pos[index], nsig=1, origin=bs$cond[index], ori.pos=bs$pos[index])
        }
        #pair
        if(sum(cluster.size==2)>0){
            index <- cum.size[cluster.size==2]
            dist42 <- bs$pos[index]-bs$pos[(index-1)]
            if(sum(dist42<=in.distance)>0){
                index2 <- index[dist42<=in.distance]
                res <- rbind(res, data.frame(pos=sapply(index2, function(x) round(weighted.mean(bs$pos[(x-1):x], w=bs$weight[(x-1):x]))), 
                    nsig=2, 
                    origin=sapply(index2, function(x) paste(bs$cond[(x-1)], bs$cond[x], sep=",")),
                    ori.pos=sapply(index2, function(x) paste(bs$pos[(x-1)], bs$pos[x], sep=","))))
            }
            if(sum(dist42>in.distance)>0){
                index2 <- index[dist42>in.distance]
                ind <- bs$weight[index2] > bs$weight[(index2-1)]
                index3 <- c(index2[ind], (index2-1)[!ind])
                res <- rbind(res, data.frame(pos=bs$pos[index3], 
                    nsig=1, 
                    origin=bs$cond[index3],
                    ori.pos=bs$pos[index3]))
            }
        }
        if(sum(cluster.size>2)>0){
            end.index <- cum.size[cluster.size>2]
            start.index <- end.index - cluster.size[cluster.size>2]+1
            for(i in 1:length(start.index)){
                res <- rbind(res, get.cluster(bs[start.index[i]:end.index[i],], in.distance=in.distance, out.distance=out.distance))
            }
        }
        res <- res[order(res$pos, decreasing = FALSE),]
        consensus.site[[chr]] <- res
        cat(".")
    }
    cat("done\n")
    return(consensus.site)
}


get.cluster <- function(bsc, in.distance=100, out.distance=250){
    D <- dist(bsc$pos)
    out <- hclust(D, method="centroid")
    #out <- hclust(D,method="average")
    cid <- cutree(out,h=in.distance)
    #for each cluster compute the weighted mean
    clusters=unique(cid)
    K=length(clusters)
    if(K==1){
        return(data.frame(pos=round(mean(bsc$pos)), nsig=nrow(bsc), origin=paste(bsc$cond, collapse=","), ori.pos=paste(bsc$pos, collapse=",")))
    }
    pos <- rep(-1, K)
    pos.weight <- rep(-1, K)
    nsig <- rep(-1, K)
    origin <- rep("", K)
    ori.pos <- rep("", K)
    for(i in 1:K){
        nsig[i] <- sum(cid==i)
        pos[i] <- round(weighted.mean(bsc$pos[cid==i], w=bsc$weight[cid==i]))
        pos.weight[i] <- max(bsc$weight[cid==i])
        origin[i] <- paste(bsc$cond[cid==i], collapse=",")
        ori.pos[i] <- paste(bsc$pos[cid==i], collapse=",")
    }
    od <- order(pos, decreasing = FALSE)
    pos <- pos[od]
    binding.info <- data.frame(cbind(nsig, origin, ori.pos))
    binding.info <- binding.info[od, ]
    pos.weight <- pos.weight[od]
    index <- which(diff(pos) < out.distance)
    if(length(index)>0){
        ind <- pos.weight[index] < pos.weight[(index+1)]
        index <- c(index[ind], (index+1)[!ind])
        pos <- pos[-index]
        binding.info <- binding.info[-index, ]
    }
    return(data.frame(cbind(pos, binding.info)))
}

get.site.count.hist <- function(pos, window.size=250){
    len <- length(pos)
    if(len==0) stop("no pos in get.site.count.hist")
    #additonal -1 is because hist is right-inclusive
    left.end <- pos - window.size-1
    bk.l <- unique(sort(c(0, left.end, pos-1, 3e9)))
    if(len>1){
        #in case the pos to the left is too close
        left.end.final <- pmax(left.end, c(0, pos[1:(len-1)]-1))
    }else{
        left.end.final <- max(left.end, 0)
    }
    index.l <- sapply(left.end.final, function(x) which(x==bk.l))
    
    right.end <- pos + window.size
    bk.r <- unique(sort(c(0, pos, right.end, 3e9)))
    if(len>1){
        right.end.final <- pmin(right.end, c(pos[-1], 3e9))
    }else{
        right.end.final <- right.end
    }
    index.r <- sapply(right.end.final, function(x) which(x==bk.r))-1

    return(list(bk.l=bk.l, index.l=index.l, bk.r=bk.r, index.r=index.r))

}

#for each binding site s, count the the 5' ends on the positive strand within the upstream window [s-w, s-1] 
#and the number of 5' ends on the negative strand within the downstream window [s+1, s+w]
#input:
#chr.data, read data for a chromosome in the form of a list, one element is for positive strand (+), and the other for negative strand (-)
#pos, vector of binding site positiions
#window.size, parameter w that determines the size of the upstream/downstrem window, 
#should be set slightly larger than the estimated average fragment length, default 250 bp
get.chr.site.count <- function(chr.data, hi){
    return(hist(chr.data[["+"]], breaks=hi$bk.l, plot = FALSE)$counts[hi$index.l]+
        hist(chr.data[["-"]], breaks=hi$bk.r, plot = FALSE)$counts[hi$index.r])
}

get.site.count <- function(dat, subtract.input=FALSE, window.size=250){
    cat("count ChIP reads around each binding site")
    consensus.site <- dat$consensus.site
    res <- c()
    for(chr in names(consensus.site)){
        hi <- get.site.count.hist(consensus.site[[chr]]$pos, window.size=window.size)
        chip.count.vec <- sapply(dat$chip.list, function(x) get.chr.site.count(x[[chr]], hi))
        
        if(length(consensus.site[[chr]]$pos)==1) chip.count.vec <- t(as.matrix(chip.count.vec))
        if(!is.null(dat$input.list) && subtract.input){
            input.count.vec <- sapply(names(dat$chip.list), function(x) get.chr.site.count(dat$input.list[[dat$matching.input.names[x]]][[chr]], hi))
            chip.count.vec <- chip.count.vec - t(t(input.count.vec)*dat$norm.factor.vec)
            chip.count.vec <- pmax(round(chip.count.vec), 0)
        }
        rownames(chip.count.vec) <- paste(chr, "_", consensus.site[[chr]]$pos, sep="")
        res <- rbind(res, chip.count.vec)
        cat(".")
    }
    dat$site.count <- res
    cat("done\n")
    return(dat)
}

#############################################################
get.unique.read.count.hist <- function(pos, window.size=250, buffer.size=50){

    temp <- sapply(pos, function(x) (x - window.size-buffer.size):(x+buffer.size))
    temp <- unique(sort(c(as.vector(temp), 3e9)))
    if(temp[1]>0){
        bk.l <- c(0, temp)
        bk.r <- c(0, temp+window.size)
        left.end <- pos - window.size-buffer.size+1
        index.l <- sapply(left.end, function(x) which(x==bk.l))
        index.r <- index.l
    }else{
        #the only way temp[1] can be < 0 is by including 0
        bk.l <- temp
        left.end <- pos - window.size-buffer.size+1
        index.l <- sapply(left.end, function(x) which(x==bk.l))
        bk.r <- temp+window.size
        if(bk.r[1]>0){
            bk.r <- c(0, bk.r)
            index.r <- index.l+1
        }else{
            index.r <- index.l
        }
    }
    
    return(list(bk.l=bk.l, index.l=index.l, bk.r=bk.r, index.r=index.r, len=window.size+2*buffer.size))

}

get.chr.unique.read.count <- function(chr.data, hi){
    his <- hist(chr.data[["+"]], breaks=hi$bk.l, plot = FALSE)$counts
    uc <- sapply(hi$index.l, function(x) sum(his[x:(x+hi$len-1)]>0))
    his <- hist(chr.data[["-"]], breaks=hi$bk.r, plot = FALSE)$counts
    uc <- uc + sapply(hi$index.r, function(x) sum(his[x:(x+hi$len-1)]>0))
    return(uc)
}

get.unique.read.count <- function(dat, window.size=250, buffer.size=50){
    cat("count unique read count around each binding site")
    consensus.site <- dat$consensus.site
    res <- c()
    for(chr in names(consensus.site)){
        hi <- get.unique.read.count.hist(consensus.site[[chr]]$pos, window.size=window.size, buffer.size=buffer.size)
        chip.count.vec <- sapply(dat$chip.list, function(x) get.chr.unique.read.count(x[[chr]], hi))
        if(length(consensus.site[[chr]]$pos)==1) chip.count.vec <- t(as.matrix(chip.count.vec))
        rownames(chip.count.vec) <- paste(chr, "_", consensus.site[[chr]]$pos, sep="")
        res <- rbind(res, chip.count.vec)
        cat(".")
    }
    dat$site.unique.read.count <- res
    cat("done\n")
    return(dat)
}
#####################################################

test.diff.binding  <- function(dat, lib.size=NULL, dispersion=NULL, common.disp=TRUE, prior.n=10, two.sample.method="composite.null", allowable.FC=1.5, collapsed.quant=0.5){
    site.count <- dat$site.count
    conds <- dat$conds
    if(is.null(lib.size)) lib.size <- median.ratio(site.count)
    #normalize lib.size such that its product is 1
    lib.size <- lib.size/exp(mean(log(lib.size)))
    dat$lib.size <- lib.size
    if(length(unique(conds)) == ncol(site.count) & is.null(dispersion)){ #no replicate
        if(ncol(site.count)==2){
            if(two.sample.method=="composite.null"){
                cat("Testing two sample without replicates, resort to testing composite null with allowable fold change", allowable.FC, "\n")
                dat$test.stat <- test.diff.binding.2sample(site.count, conds, lib.size=lib.size, allowable.FC=allowable.FC)
                return(dat)
            }else{
                if(!is.numeric(dispersion)) stop("Two sample with no replicates, you need to use default method or provide a dispersion parameter.")
            }
        }else{ #more than 2 conditions
            if(!is.numeric(dispersion)){
                cat("Estimating common dispersion by collapsed quantile method.\n")
                dispersion <- est.common.disp.by.collapse(site.count, lib.size, quant=collapsed.quant)
                cat("Common dispersion at", collapsed.quant, "quantile:", dispersion, "\n")
            }
        }
    }
    dat$test.stat <- test.diff.binding.edgeR(site.count, conds, lib.size=lib.size, dispersion=dispersion, common.disp=common.disp, prior.n=prior.n)
    return(dat)
}

test.diff.binding.2sample <- function(site.count, conds, lib.size, allowable.FC=1.5){
    if(setequal(names(lib.size), colnames(site.count))) lib.size <- lib.size[colnames(site.count)]
    total <- rowSums(site.count)
    prop <- site.count[,1]/total
    ratio <- lib.size[1]/lib.size[2]
    lp <- ratio/allowable.FC
    lp <- lp/(1+lp)
    up <- ratio*allowable.FC
    up <- up/(1+up)
    prop[prop < lp] <- lp
    prop[prop > up] <- up
    data <- cbind(site.count[,1], total, prop)
    pval <- apply(data, 1, function(x) binom.test(x[1], n=x[2], p=x[3], alternative="two.sided")$p.value)
    ratio <- (site.count[,2]/lib.size[2]+1)/(site.count[,1]/lib.size[1]+1)
    res <- data.frame(ratio, pval)
    colnames(res) <- c(paste("FC.", conds[2], sep=""), "pval")
    return(res)
}

est.common.disp.by.collapse <- function(site.count, lib.size, quant=0.5){
    n <- ncol(site.count)
    clp <- combn(n ,2)
    disp <- rep(-1, ncol(clp))
    for(i in 1:ncol(clp)){
        conds <- 1:n
        conds[clp[2,i]] <- conds[clp[1,i]]
        d <- DGEList(counts=site.count, lib.size=lib.size, group=conds)
        n.cond <- length(unique(conds))
        if(n.cond==2){
            d <- estimateCommonDisp(d)
        }else{#multi group
            conds <- factor(conds)
            design <- model.matrix(~conds)
            d <- estimateGLMCommonDisp(d, design)
        }
        disp[i] <- d$common.dispersion
    }
    return(quantile(disp, quant))
}

test.diff.binding.edgeR <- function(site.count, conds, lib.size, dispersion=NULL, common.disp=TRUE, prior.n=10){
    require(edgeR)
    d <- DGEList(counts=site.count, lib.size=lib.size, group=conds)
    #d <- calcNormFactors(d)
    n.cond <- length(unique(conds))
    if(n.cond==2){
        if(!is.null(dispersion)) {
            res <- exactTest(d, dispersion=dispersion)
        }else{
            d <- estimateCommonDisp(d)
            if(common.disp){
                cat("Common dispersion:", d$common.dispersion, "\n")
            }else{
                d <- estimateTagwiseDisp(d, prior.n = prior.n, verbose=FALSE)
            }
            #
            res <- exactTest(d)
        }
    }else{#multi group
        conds <- factor(conds)
        design <- model.matrix(~conds)
        colnames(design)[2:n.cond] <- levels(conds)[2:n.cond]
        if(!is.null(dispersion)) {
            glmfit <- glmFit(d, design, dispersion = dispersion)
        }else{
            if(common.disp){
                d <- estimateGLMCommonDisp(d, design)
                cat("Common dispersion:", d$common.dispersion, "\n")
                glmfit <- glmFit(d, design, dispersion = d$common.dispersion)
            }else{
                d <- estimateGLMTrendedDisp(d, design)
                d <- estimateGLMTagwiseDisp(d, design)
                glmfit <- glmFit(d, design, dispersion = d$tagwise.dispersion)
            }
        }
        res <- glmLRT(d, glmfit, coef = 2:length(unique(conds)))
    }
    norm.count <- t(t(site.count)/lib.size)
    #add one to avoid dividing by zero or tiny count
    mean.count <- sapply(levels(conds), function(x){
        if(sum(conds==x)==1){
            norm.count[,conds==x]+1
        }else{
            rowMeans(norm.count[,conds==x])+1
        }
    })
    FC <- mean.count[,-1]/mean.count[,1]
    res <- data.frame(FC, pval=res$table$PValue)
    colnames(res) <- c(paste("FC.", levels(conds)[-1], sep=""), "pval")
    rownames(res) <- rownames(site.count)
    return(res)
}

#only work for two conditions
test.diff.binding.DESeq <- function(dat, lib.size=NULL){
    site.count <- dat$site.count
    conds <- dat$conds
    if(length(levels(dat$conds))!=2) stop("DESeq only handles two conditions comparison.")
    if(is.null(lib.size)) lib.size <- median.ratio(site.count)
    lib.size <- lib.size/exp(mean(log(lib.size)))
    dat$lib.size <- lib.size

    require(DESeq)
    cds <- newCountDataSet(countData=site.count, conditions=conds)
    pData(cds)$sizeFactor <- lib.size
    if(ncol(site.count)==2){#no rep
        cds <- estimateDispersions(cds, method="blind")
    }else{
        cds <- estimateDispersions(cds)
    }
    pval <- nbinomTest( cds, levels(dat$conds)[1], levels(dat$conds)[2], pvals_only=TRUE)
    pval[is.na(pval)] <- 1

    norm.count <- t(t(site.count)/lib.size)
    #add one to avoid dividing by zero or tiny count
    mean.count <- sapply(levels(conds), function(x){
        if(sum(conds==x)==1){
            norm.count[,conds==x]+1
        }else{
            rowMeans(norm.count[,conds==x])+1
        }
    })
    FC <- mean.count[,2]/mean.count[,1]
    res <- data.frame(FC, pval)
    colnames(res) <- c(paste("FC.", levels(conds)[2], sep=""), "pval")
    dat$test.stat <- res
    return(dat)
}

unlist.data.frame <-
function(sth.list, add.row.name=TRUE, list.name="chr", use.list.name = TRUE){
    res <- c()
    for(index in names(sth.list)){
        if(use.list.name){
            res <- rbind(res, cbind(data.frame(index), sth.list[[index]]))
        }else{
            res <- rbind(res, sth.list[[index]])
        }
    }
    if(use.list.name) colnames(res)[1] <- list.name
    if(add.row.name) rownames(res) <- paste(res[,1], "_", res[,2], sep="")
    return(res)
}

#report peak by FDR (by BH or adaptive) or top n
#return chr, loc, pval, FDR and FC
report.peak <- function(test.res, FDR=NULL, FDR.method="BH", n=10, add.origin=TRUE, adaptive.threshold=c(0.05, 0.95)){
    qval <- fdr(test.res$test.stat$pval, method=FDR.method, adaptive.threshold=adaptive.threshold)
    if(add.origin){
        res <- data.frame(cbind(unlist.data.frame(test.res$consensus.site), test.res$test.stat), FDR=qval)
    }else{
        res <- data.frame(cbind(unlist.data.frame(test.res$consensus.site)[, c("chr", "pos")], test.res$test.stat), FDR=qval)
    }
    res <- res[order(res$pval),]
    rownames(res) <- NULL
    if(is.null(FDR)){
        return(res[1:min(n, nrow(res)), ])
    }else{
        if(sum(res$FDR<=FDR)==0) cat("FDR too strigent, no result is returned. The minimum FDR obtained is ", min(res$FDR))
        return(res[res$FDR<=FDR, ])
    }
}


load.data <- function(chip.data.list, conds, consensus.site, input.data.list=NULL,
    data.type="MCS", chr.vec=NULL, chr.exclusion=NULL, chr.len.vec=NULL,
    norm.factor.vec=NULL, frag.len=200){

    conds <- as.factor(conds)
    if(length(conds)!=length(chip.data.list)) stop("conds need to have equal length as chip.data.list")
    if(!setequal(names(chip.data.list), names(conds))){
        names(conds) <- names(chip.data.list)
    }else{
        #put conds in the same order as chip
        conds <- conds[names(chip.data.list)]
    }
#check the name in chip.data.list, conds, input.data.list match
#if(sum(names(input.data.list) %in% union(names(chip.data.list), conds)) < length(input.data.list)) stop("input.data.list names need to match condition name (conds) or chip.data.list names.")
    if(!is.null(input.data.list)){
        matching.input.names <- sapply(names(chip.data.list), 
            function(x) {
                if(sum(x==names(input.data.list))==1) return(x)
                cname <- as.character(conds[[x]])
                if(sum(cname==names(input.data.list))==1){
                    return(cname)
                }else{
                    stop("The name of a input need to be uniquely matched to the name of a chip sample/replicate or a condition.")
                }
            }
        )
    }else{
        matching.input.names <- NULL
    }
    if(class(consensus.site)=="data.frame") consensus.site <- split(consensus.site[,setdiff(colnames(consensus.site), "chr")], consensus.site$chr, drop = TRUE)

    cat("reading data...")
    chip.list <- read.data.list(chip.data.list, data.type)
    if(!is.null(input.data.list)){
        input.list <- read.data.list(input.data.list, data.type)
    }else{
        input.list <- NULL
    }
    cat("done\n")
    
    if(is.null(chr.vec)){
        chr.vec <- get.chr(chip.list, input.list, chr.exclusion=chr.exclusion)
    }
    if(length(setdiff(names(consensus.site), chr.vec))>0) stop("chr.vec don't contain all chromosomes in binding sites")
    if(is.null(chr.len.vec)){
        chr.len.vec <- get.chr.len(chip.list, input.list, chr.vec)
    }

    res <- list(chip.list=chip.list, conds=conds, consensus.site=consensus.site, 
        input.list=input.list, matching.input.names=matching.input.names, 
        chr.vec=chr.vec, chr.len.vec=chr.len.vec, frag.len=frag.len)
    if(!is.null(input.list)){
        if(is.null(norm.factor.vec)){
            cat("computing normalization factor between ChIP and control samples")
            if(is.null(consensus.site)){
                res$norm.factor.vec <- get.NCIS.norm.factor(res)        
            }else{
                bkg <- comp.background.size(res, consensus.site, binsize=500, shift.size=round(frag.len/2))
                res$chip.background.size <- bkg$chip.background.size
                res$input.background.size <- bkg$input.background.size
                res$norm.factor.vec <- sapply(names(chip.list), function(x) res$chip.background.size[[x]]/res$input.background.size[[matching.input.names[x]]])
            }
            cat("done\n")
        }else{
            if(is.null(names(norm.factor.vec))) names(norm.factor.vec) <- names(chip.list)
            res$norm.factor.vec <- norm.factor.vec
        }
    }
    return(res)
}

plotPeak <- function(rept, dat, lib.size=NULL, w=400, ext=200, combine.rep=FALSE, cap=NULL, n.row.per.page=6, caption=NULL){
    if(!is.null(cap)) if(cap<=1) stop("cap need to be at least 1.")
    ext <- round(ext)
    if(ext<1) stop("ext need to be at least 1.")
    if(is.null(lib.size)){
        if(is.null(dat$lib.size)){
            stop("User need to provide library sizes for the ChIP data.")
        }else{
            lib.size <- dat$lib.size
        }
    }else{
        if(length(lib.size)!=length(dat$chip.list)) stop("lib.size need to have same length as ChIP libraries.")
    }
    if(length(w)!=nrow(rept)) w <- rep(w, length.out=nrow(rept))
    if(length(caption)!=nrow(rept)) caption <- rep(caption, length.out=nrow(rept))
    rep(1:3, length.out=5)
    if(combine.rep){
        dat <- combine.reps(dat, lib.size)
        lib.size <- dat$lib.size
    }
    for(i in 1:nrow(rept)){
        single.peak(as.character(rept$chr[i]), as.integer(rept$pos[i]), dat, lib.size, w=w[i], ext=ext, cap=cap, n.row.per.page=n.row.per.page, caption=caption[i])
    }
}

# plot peak for multiple chip samples
# assume following exist: chip.list, input.list, lib.size, norm.factor.vec
# input parameter:
# chr: chromosome
# center: binding site
# w: half window size
# ext: length of extension for the coverage
# center.line: whether draw center line
# ind: indicator whether to draw center line in red
single.peak <- function(chr, center, dat, lib.size, w=300, ext=200, center.line=TRUE, ind=FALSE, cap=NULL, n.row.per.page=6, caption=NULL){
    ext <- round(ext)
    if(ext<1) stop("ext need to be at least 1.")
    start <- center-w
    end <- center+w
    scale <- max(lib.size)/lib.size
    coverage <- list()
    for(s in names(dat$chip.list)){
        fp <- scale[[s]]*get.coverage(start+1, end+1, dat$chip.list[[s]][[chr]][["+"]], ext=ext, cap=cap)
        rp <- scale[[s]]*get.coverage(start-1, end-1, dat$chip.list[[s]][[chr]][["-"]], "D", ext=ext, cap=cap)
        temp <- cbind(fp, rp)
        if(!is.null(dat$input.list)){
            fp <- scale[[s]]*dat$norm.factor.vec[[s]]*get.coverage(start+1, end+1, dat$input.list[[dat$matching.input.names[s]]][[chr]][["+"]], ext=ext, cap=cap)
            rp <- scale[[s]]*dat$norm.factor.vec[[s]]*get.coverage(start-1, end-1, dat$input.list[[dat$matching.input.names[s]]][[chr]][["-"]], direction="D", ext=ext, cap=cap)
            temp <- cbind(temp, fp, rp)
        }
        coverage[[s]] <- temp
    }
    ymax <- max(sapply(coverage, max))+10
    ymin <- max(min(sapply(coverage, min))-10, 0)
    par(mfrow=c(min(length(dat$chip.list), n.row.per.page),1))
    par(mar=c(1, 4, 1.2, 1))
    par(oma=c(3,2.2,2.2,0))
    for(s in names(dat$chip.list)){
        if(!is.null(dat$input.list)){
            matplot(start:end, coverage[[s]], col=c("blue", "red", "green", "orange"), ylab=s, ylim=c(ymin, ymax), type = "l", lty=1)
        }else{
            matplot(start:end, coverage[[s]], col=c("blue", "red"), ylab=s, ylim=c(ymin, ymax), type = "l", lty=1)
        }
        if(center.line) abline(v=center, lty=2, col=ifelse(ind, "red", "black"))
    }
    if(ext>1) mtext("Coverage", side=2, line=1, outer=TRUE, at=.5)
    mtext("Location", side=1, line=2, outer=TRUE, at=.5)
    mtext(paste(chr, ", at: ", center, " ", caption, sep=""), side=3, line=1, outer=TRUE, at=.5)
}

get.coverage <- function(start, end, loc.vec, direction="U", ext=200, cap=NULL){
    if(start > end) stop("start > end")
    if(direction=="U"){
        rstart <- start-ext+1
        rend <- end
    }else{
        rstart <- start
        rend <- end+ext-1
    }

    bk <- c(0, (rstart-1):rend, 3e9)
    phis <- hist(loc.vec, breaks=bk, plot = FALSE)$counts
    phis <- phis[-1]
    phis <- phis[-length(phis)]
    if(!is.null(cap)) phis <- pmin(phis, cap)
    if(ext==1) return(phis)
    return(convolve(phis, rep(1,ext), type="f"))
}

#combine replicates within conditions
combine.reps <- function(dat, lib.size){
    lev.vec <- levels(dat$conds)
    chip.list <- list()
    if(is.null(dat$input.list)){
        input.list <- NULL
        matching.input.names <- NULL
    }else{
        input.list <- list()
        matching.input.names <- lev.vec
        names(matching.input.names) <- lev.vec
    }
    lib.size.all <- rep(-1, length(lev.vec))
    names(lib.size.all) <- lev.vec
    chip.background.size <- lib.size.all
    input.background.size <- lib.size.all
    for(lev in lev.vec){
        chip.name <- names(dat$conds)[as.character(dat$conds)==lev]
        chip.list[[lev]] <- MCS.list.merge(dat$chip.list[chip.name])
        chip.background.size[lev] <- sum(dat$chip.background.size[chip.name])
        if(!is.null(dat$input.list)){
            input.name <- unique(dat$matching.input.names[chip.name])
            input.list[[lev]] <- MCS.list.merge(dat$input.list[input.name])
            input.background.size[lev] <- sum(dat$input.background.size[input.name])
        }
        lib.size.all[lev] <- sum(lib.size[chip.name])
    }
    conds <- factor(lev.vec, levels=lev.vec)
    names(conds) <- lev.vec
    if(is.null(dat$input.list)){
        norm.factor.vec <- NULL
    }else{
        norm.factor.vec <- chip.background.size[lev.vec]/input.background.size[lev.vec]
    }
    
    return(list(chip.list=chip.list, conds=conds, input.list=input.list, matching.input.names=matching.input.names,
        chr.vec=dat$chr.vec, chr.len.vec=dat$chr.len.vec, 
        consensus.site=dat$consensus.site, norm.factor.vec=norm.factor.vec, lib.size=lib.size.all, frag.len=dat$frag.len, 
        chip.background.size=chip.background.size, input.background.size=input.background.size))
}

MCS.list.merge <- function(lst){
    for(i in 1:length(lst)){
        if(i==1){
            temp <- lst[[1]]
        }else{
            temp <- MCS.merge(temp, lst[[i]])
        }
    }
    return(temp)
}

MCS.merge <- function(rep1, rep2){
    for(chr in names(rep1)){
        rep1[[chr]][["+"]] <- c(rep1[[chr]][["+"]], rep2[[chr]][["+"]])
        rep1[[chr]][["-"]] <- c(rep1[[chr]][["-"]], rep2[[chr]][["-"]])
    }
    return(rep1)
}

read.MCS.list <- function(mnames, mfiles, dir="", merge.by.name=TRUE, chr.exclusion=NULL, suffix.filter="random", suffix.stripper=".fa", sort.chr=TRUE){
    if(length(mnames)!=length(mfiles)) stop("The lengths of names and files are not the same.")
    if(length(mnames) > length(unique(mnames)) & !merge.by.name) stop("Duplicate names. Name need to be unique if merge.by.name is set to FALSE.")
    res <- list()
    if(nchar(dir)>0) if(substring(dir, nchar(dir))!="/") dir <- paste(dir, "/", sep="")
    for(i in 1:length(mnames)){
        load(paste(dir, mfiles[i], sep=""))
        if(!is.null(suffix.filter)) mcs <- mcs[!(substring(names(mcs), nchar(names(mcs))-nchar(suffix.filter)+1) == suffix.filter)]
        if(!is.null(suffix.stripper)){
            ind <- substring(names(mcs), nchar(names(mcs))-nchar(suffix.stripper)+1) == suffix.stripper
            if(sum(ind)>0) names(mcs)[ind] <- substr(names(mcs)[ind], 1, nchar(names(mcs))-nchar(suffix.stripper))
        }
        #do it again, in case suffix is added before stripper
        if(!is.null(suffix.filter)) mcs <- mcs[!(substring(names(mcs), nchar(names(mcs))-nchar(suffix.filter)+1) == suffix.filter)]
        if(!is.null(chr.exclusion)) mcs <- mcs[!(names(mcs) %in% chr.exclusion)]
        if(is.null(res[[mnames[i]]])){
            res[[mnames[i]]] <- mcs
        }else{
            cat("merging replicates for condition", mnames[i], "\n")
            res[[mnames[i]]] <- MCS.merge(res[[mnames[i]]], mcs)
        }
    }
    if(sort.chr) res <- lapply(res, function(x) x[sort(names(x))])
    return(res)
}

#method: BH or adaptive
fdr <- function (p, method, n = length(p), adaptive.threshold=c(0.05, 0.95)) 
{
    if (method == "fdr") 
        method <- "BH"
    p0 <- p
    if (all(nna <- !is.na(p))) 
        nna <- TRUE
    p <- as.vector(p[nna])
    lp <- length(p)
    stopifnot(n >= lp)
    if (n <= 1) 
        return(p0)
    p0[nna] <- switch(method, 
    BH = {
        i <- lp:1
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        pmin(1, cummin(n/i * p[o]))[ro]
    }, adaptive = {
        i <- lp:1
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        n0 <- min((sum(p > adaptive.threshold[1] & p <= adaptive.threshold[2])+1)/(adaptive.threshold[2]-adaptive.threshold[1]), lp)
        #n0 <- min((sum(p>lambda)+1)/(1-lambda), lp)
        pmin(1, cummin(n0/i * p[o]))[ro]
    }, BY = {
        i <- lp:1
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        q <- sum(1/(1L:n))
        pmin(1, cummin(q * n/i * p[o]))[ro]
    }, none = p)
    p0
}
