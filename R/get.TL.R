get.TL <-
function(rast.stack){
    p=nlayers(rast.stack)
    M=ncell(rast.stack[[1]])
    B=list()
    for(i in 1:p){
        B[[i]]=Matrix(0,nrow=M,ncol=M)
    }
    n.x=ncol(rast.stack[[1]])
    n.y=nrow(rast.stack[[1]])
    N=n.x*n.y

    xy=rowColFromCell(rast.stack[[1]],1:M)

    ## left cell
    idx=which(xy[,2]>1)
    left.cells=cellFromRowCol(rast.stack[[1]],xy[idx,1],xy[idx,2]-1)
    B.idx=idx+N*(left.cells-1)
    for(i in 1:p){
        vals=.5*(extract(rast.stack[[i]],idx)+extract(rast.stack[[i]],
            left.cells))
        idx.nonzero=which(vals!=0)
        B[[i]]=B[[i]]+sparseMatrix(i=idx[idx.nonzero],
             j=left.cells[idx.nonzero],x=vals[idx.nonzero],dims=c(N,N))
    }

    ## right cell
    idx=which(xy[,2]<n.x)
    right.cells=cellFromRowCol(rast.stack[[1]],xy[idx,1],xy[idx,2]+1)
    B.idx=idx+N*(right.cells-1)
    for(i in 1:p){
        vals=.5*(extract(rast.stack[[i]],idx)+extract(rast.stack[[i]],right.cells))
        idx.nonzero=which(vals!=0)
        B[[i]]=B[[i]]+sparseMatrix(i=idx[idx.nonzero],j=right.cells[idx.nonzero],x=vals[idx.nonzero],dims=c(N,N))
    }

    ## up cell
    idx=which(xy[,1]>1)
    up.cells=cellFromRowCol(rast.stack[[1]],xy[idx,1]-1,xy[idx,2])
    B.idx=idx+N*(up.cells-1)
    for(i in 1:p){
        vals=.5*(extract(rast.stack[[i]],idx)+extract(rast.stack[[i]],up.cells))
        idx.nonzero=which(vals!=0)
        B[[i]]=B[[i]]+sparseMatrix(i=idx[idx.nonzero],j=up.cells[idx.nonzero],x=vals[idx.nonzero],dims=c(N,N))
    }

    ## down cell
    idx=which(xy[,1]<n.y)
    down.cells=cellFromRowCol(rast.stack[[1]],xy[idx,1]+1,xy[idx,2])
    B.idx=idx+N*(down.cells-1)
    for(i in 1:p){
        vals=.5*(extract(rast.stack[[i]],idx)+extract(rast.stack[[i]],down.cells))
        idx.nonzero=which(vals!=0)
        B[[i]]=B[[i]]+sparseMatrix(i=idx[idx.nonzero],j=down.cells[idx.nonzero],x=vals[idx.nonzero],dims=c(N,N))
    }
    B
}
