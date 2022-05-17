if (F) {
    #Author: zhuzhiyong
    #mail: zhuzhiyong@genomics.cn
    #create time   : Tue 17 May 2022 02:38:27 PM CST
    # object: input your matrix with seurat object
    # organize: with column in meta.date(cell-level) you want group by
    # cell_n: cell number define to condense
    # method: aggregate expression by average
    # chunk_size: split datasets to chunk, default 10000 cells a chunk 
}
#define a pseudo_cell func
pseudo_cell <- function(object,organize="orig.ident",cell_n=10,method = 'average',chunk_size=100000){
    average_expr <- function(myobject,categorie,chunk.no=1){
        if( ncol(myobject) >1){
            ##generate category.matrix
            category <- myobject@meta.data[,'pseudo_label']
            if (length(table(category)) >1 ) {
                category.matrix <- Matrix::sparse.model.matrix( ~ -1 + category )
                colnames(category.matrix) <- gsub("category","",colnames(category.matrix))
                # colnames(category.matrix) <- paste0(colnames(category.matrix),".",chunk.no)
            }else{
                category.matrix <- matrix(data = 1, nrow = ncol(myobject),dimnames = list(Cells(myobject), paste0(categorie,"_N.1.",chunk.no)))
            }
            ##aggregate cell by average expr
            if (method == 'average') {category.matrix <- sweep(category.matrix, 2,Matrix::colSums(category.matrix), FUN = "/") }
            data.use <- GetAssayData(myobject, slot = "counts")
            aggreate.counts = data.use %*% category.matrix
        }else{
            aggreate.counts <- myobject[["RNA"]]@counts
            colnames(aggreate.counts) <- paste0(categorie,"_N.1.",chunk.no);
        }
        return(aggreate.counts)
    }
    categories=names(table(object@meta.data[,organize]))
    data.list <- list()
    for ( idx in 1:length(categories)) {
        categorie = categories[idx]
        subObj <- object[, object@meta.data[, organize] == categorie]
        ####如果子集太大，cut chunks to process
        if( ncol(subObj) > chunk_size){
            subObj@meta.data <- subObj@meta.data %>% mutate(cellN=colnames(subObj)) %>% group_by(get(organize)) %>% 
                mutate( chunk_label = rep(1:ceiling(n()/chunk_size), each=chunk_size, length.out=n()) ) %>% 
                tibble::column_to_rownames('cellN') %>% data.frame()
            data.list.tmp <- list() #用于存储chunk处理某一类细胞类型的返回数据
            for (chunk in 1:length(table(subObj$chunk_label))) {
                subObj.chunk <- subset(subObj,subset=chunk_label==chunk)
                ###labeling group number on cell-level meta.data
                subObj.chunk@meta.data <- subObj.chunk@meta.data %>% mutate(cellN=colnames(subObj.chunk)) %>% group_by(get(organize)) %>% 
                    mutate( bins.no = rep(1:ceiling(n()/cell_n), each=cell_n, length.out=n()) ) %>% 
                    mutate(pseudo_label=paste0(get(organize),"_N.",bins.no,".",chunk)) %>% tibble::column_to_rownames('cellN') %>% data.frame()
                
                agg.counts <- average_expr(subObj.chunk,categorie=categorie,chunk.no = chunk)
                data.list.tmp[[chunk]] <- CreateSeuratObject(counts = agg.counts)
                data.list.tmp[[chunk]]$orig.ident <- categorie
            }
            toRet.chunk = data.list.tmp[[chunk]]
            if (length(data.list.tmp) > 1) {toRet.chunk <- merge(data.list.tmp[[1]],data.list.tmp[2:length(data.list.tmp)])}
            data.list[[idx]] <- DietSeurat(toRet.chunk)
            
        }else{#小于chunk size的数据集
            subObj@meta.data <- subObj@meta.data %>% mutate(cellN=colnames(subObj)) %>% group_by(get(organize)) %>% 
                mutate( bins.no = rep(1:ceiling(n()/cell_n), each=cell_n, length.out=n()) ) %>% 
                mutate(pseudo_label=paste0(get(organize),"_N.",bins.no)) %>% tibble::column_to_rownames('cellN') %>% data.frame()
            agg.counts <- average_expr(subObj,categorie=categorie,chunk.no = 1)
            data.list[[idx]] <- CreateSeuratObject(counts = agg.counts)
            data.list[[idx]]$orig.ident <- categorie #返回小chunk的obj
        }
    }
    toRet =  data.list[[idx]]
    if (length(data.list) > 1) {toRet <- merge(data.list[[1]],data.list[2:length(data.list)])}
    toRet <- DietSeurat(toRet)
    return(toRet)
}
