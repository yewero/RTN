##This function is used for argument checking
tnai.checks <- function(name, para) {
  if(name=="alpha") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>1 || para<0)
      stop("'alpha' should be an integer or numeric value <=1 and >=0 !\n",call.=FALSE)
  }
 else if(name=="maxgap") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<0 || as.integer(para)<para)
      stop("'maxgap' should be an integer value >=0!\n",call.=FALSE)
  }
 else if(name=="tna.what"){
    opts<-c("summary","status","para",
            "pheno","hits", "regulatoryElements", 
            "tnet", "refnet",
            "regulons","refregulons",
            "regulons.and.mode", "refregulons.and.mode",
            "rowAnnotation", "colAnnotation",
            "mra","gsea1", "gsea2",
            "overlap","synergy","shadow")
    supps <- c("nondpiregulons.and.mode","regulons.and.pheno","refregulons.and.pheno")
    if(!is.character(para) || length(para)!=1 || !(para %in% c(opts,supps)))
      stop(paste("'what' should be any one of the options: \n", paste(opts,collapse = ", ") ),call.=FALSE )
  }
 else if(name=="tni.what"){
    opts<-c("summary","status", "para",
            "gexp","regulatoryElements",
            "tnet","refnet",
            "regulons","refregulons",
            "regulons.and.mode", "refregulons.and.mode",
            "rowAnnotation", "colAnnotation", "cdt","cdtrev")
    supps <- c("regulons.and.mode.gmm","refregulons.and.mode.gmm")
    if(!is.character(para) || length(para)!=1 || !(para %in% c(opts,supps)))
      stop(paste("'what' should be any one of the options: \n", paste(opts,collapse = ", ") ) ,call.=FALSE )
  }
 else if(name=="avs.what"){
    opts<-c("markers","validatedMarkers","variantSet","randomSet","summary",
            "status","linkedMarkers","randomMarkers","vse","evse","pevse",
            "annotation.vse","annotation.evse","annotation.pevse")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'what' should be any one of the options: \n", paste(opts,collapse = ", ") ),call.=FALSE )
  }
 else if(name=="avs.plot.what"){
    opts<-c("vse","evse","pevse")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'avs.plot.what' should be any one of the options: \n", paste(opts,collapse = ", ") ),call.=FALSE )
  }
 else if(name=="reldata"){
    opts<-c("RTNdata.LDrel27","RTNdata.LDHapMapRel27","RTNdata.LD1000gRel20130502")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("available 'reldata':", paste(opts,collapse = ", ") ),call.=FALSE)
  }
 else if(name=="tni.gtype"){
    opts<-c("rmap","amap","amapDend","mmap","mmapDetailed")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'gtype' should be any one of the options:", paste(opts,collapse = ", ") ),call.=FALSE )
  }
 else if(name=="hcl"){
    if( !is.null(para) && !class(para)=="hclust" )
      stop("'hcl' should be an 'hclust' object whith TF's IDs!",call.=FALSE)
  }
 else if(name=="tna.gtype"){
    opts<-c("rmap","amap")
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'gtype' should be any one of the options:", paste(opts,collapse = ", ") ),call.=FALSE )
  }
 else if(name=="order") {
    if(!is.logical(para) || length(para)!=1)
      stop("'order' should be a logical value!",call.=FALSE)
  }
 else if(name=="regulon.order") {
    opts<-c('name','score','size', 'pvalue','adj.pvalue','none')
    if(!is.character(para) || length(para)!=1 || !(para %in% opts))
      stop(paste("'regulon.order' should be any one of the options:", paste(opts,collapse = ", ") ),call.=FALSE )
    para<-switch(which(opts%in%para),'Regulon','Observed.Score','Regulon.Size','Pvalue','Adjusted.Pvalue','none')
    return(para)
  }
 else if(name=="ntop") {
    if(!is.null(para) && ( !(is.numeric(para) || is.integer(para)) || length(para)!=1 || round(para,0)!=para) )
      stop("'ntop' should be an integer value!",call.=FALSE)
  }
 else if(name=="estimator"){
    if(!is.character(para) || length(para)!=1 || !(para %in% c("pearson", "kendall", "spearman")))
      stop("'estimator' should be any one of 'pearson', 'kendall' and 'spearman'!",call.=FALSE)
  }
 else if(name=="pwtransform"){
    if(!is.logical(para) || length(para)!=1)
      stop("'pwtransform' should be a logical value!",call.=FALSE)
  }
 else if(name=="medianEffect"){
    if(!is.logical(para) || length(para)!=1)
      stop("'medianEffect' should be a logical value!",call.=FALSE)
  }
 else if(name=="mdStability.custom"){ 
    if( (!is.numeric(para) && !is.integer(para)) || length(para)<1 || length(para)>2){
      stop("custom 'mdStability' should be a numeric value, or a vector of length = 2!",call.=FALSE)
    }
    if(length(para)==1){
      para<-c(-abs(para),abs(para))
    }else {
      para<-sort(para)
      if(sum(para>0)!=1)stop("custom 'mdStability' upper and lower bounds should have different signals!")
    }
    return(para)
  }
 else if(name=="reportNames"){
    if(!is.logical(para) || length(para)!=1)
      stop("'reportNames' should be a logical value!",call.=FALSE)
  }
 else if(name=="autoformat"){
    if(!is.logical(para) || length(para)!=1)
      stop("'autoformat' should be a logical value!",call.=FALSE)
  }
 else if(name=="labpair") {
    if( (!is.matrix(para) || ncol(para)!=2) )
      stop("'labpair' should be a two-column matrix with regulon pair labels!",call.=FALSE)
  }
 else if(name=="tnet"){
    if(!is.character(para) || length(para)!=1 || !(para %in% c("dpi", "ref")))
      stop("'tnet' should be any one of 'dpi' and 'ref'!",call.=FALSE)
  }  
 else if(name=="gsea.tnet"){
    if(!is.character(para) || length(para)!=1 || !(para %in% c("dpi","ref")))
      stop("'tnet' should be any one of 'dpi' and 'ref'!",call.=FALSE)
  }
 else if(name=="area.tnet"){
    if(!is.character(para) || length(para)!=1 || !(para %in% c("dpi","ref")))
      stop("'tnet' should be any one of 'dpi' and 'ref'!",call.=FALSE)
  }
 else if(name=="TNI"){
    if(!class(para)=="TNI")
      stop("'object' should be an object of class 'TNI'!",call.=FALSE)
    if( sum(names(para@results)%in%c("tn.ref", "tn.dpi") )<2 )
      stop("'object' should contain results in the slot 'result'!",call.=FALSE)
    if(!is.character(para@regulatoryElements) || 
      any(is.na(para@regulatoryElements)) || any(para@regulatoryElements==""))
      stop("'object' should contain names of regulatory elements in the slot 'regulatoryElements'!",call.=FALSE)
  }
 else if(name=="transcriptionalNetwork") {
    if( !is.matrix(para) )
      stop(" 'transcriptionalNetwork' should be a matrix with targets on rows and TFs on cols!",call.=FALSE)
    if( is.null(colnames(para)) || is.null(rownames(para)) || 
        any(duplicated(rownames(para))) || any(duplicated(colnames(para))) )
      stop("'transcriptionalNetwork' should be a matrix of named rows and cols, witg unique names!",call.=FALSE)  
    if(sum(is.na(para))>0)
      stop("'transcriptionalNetwork' matrix should have no NAs!",call.=FALSE)
  }
 else if(name=="referenceNetwork") {
    if( !is.matrix(para) )
      stop("'referenceNetwork' should be a matrix with targets on rows and TFs on cols!",call.=FALSE)
    if(  is.null(colnames(para)) || is.null(rownames(para)) || 
         any(duplicated(rownames(para))) || any(duplicated(colnames(para))) )
      stop("'referenceNetwork' should be a matrix of named rows and cols (unique names)!",call.=FALSE)  
    if(sum(is.na(para))>0)
      stop("'referenceNetwork' matrix should have no NAs!",call.=FALSE)
  }  
 else if(name=="gexp") {
    if( !is.matrix(para) || !is.numeric(para[1,]))
      stop("'gexp' should be a numeric matrix with genes on rows and samples on cols!",call.=FALSE)
    if(  is.null(colnames(para)) || is.null(rownames(para)) || 
         any(duplicated(rownames(para))) || any(duplicated(colnames(para))) )
      stop("the 'gexp' matrix should be named on rows and cols (unique names)!",call.=FALSE)
  } else if(name=="gxdata") {
    if( !is.matrix(para) || !is.numeric(para[1,]))
      stop("'gxdata' should be a numeric matrix with genes on rows and samples on cols!",call.=FALSE)
    if(  is.null(colnames(para)) || is.null(rownames(para)) || 
         any(duplicated(rownames(para))) || any(duplicated(colnames(para))) )
      stop("the 'gxdata' matrix should be named on rows and cols (unique names)!",call.=FALSE)
  }
 else if(name=="regulatoryElements") {
    if(!is.character(para) || any(is.na(para)) || any(para==""))
      stop("'regulatoryElements' should be a character vector, without 'NA' or empty names!",call.=FALSE)
    para <- para[!duplicated(para)]
    return(para)
  } else if(name=="samples") {
    if(!is.null(para)){
      if(!all.characterValues(para)){
        stop("'samples' should be a character vector, without 'NA'!",call.=FALSE)
      }
    }
  }else if(name=="features") {
    if(!is.null(para)){
      if(!all.characterValues(para)){
        stop("'features' should be a character vector, without 'NA'!",call.=FALSE)
      }
    }
  }else if(name == "refsamp"){
    if(!is.null(para)){
      if(!all.characterValues(para)){
        stop("'refsamp' should be a character vector, without 'NA'!",call.=FALSE)
      }
    }
  }else if(name=="tfs") {
    if(!is.null(para)){
    if( !is.character(para) || any(is.na(para)) || any(para==""))
      stop("'tfs' should be a character vector, without 'NA' or empty names!",call.=FALSE)
    }
  }else if(name=="mds") {
    if(!is.null(para)){
      if( !is.character(para) || any(is.na(para)) || any(para==""))
        stop("'mds' should be a character vector, without 'NA' or empty names!",call.=FALSE)
    }
  }else if(name=="modulators") {
    if(!is.null(para)){
      if( !is.character(para) || any(is.na(para)) || any(para==""))
        stop("'modulators' should be a character vector, without 'NA' or empty names!",call.=FALSE) 
    }
  }else if(name=="sampling") {
    if( !(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>100 || para<1 || round(para,0)!=para )
      stop("'sampling' should be an integer value >= 1 and <= 100 !",call.=FALSE)
  }else if(name=="probs") {
    if( !(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>1 || para<0)
      stop("'probs' should be a numeric value >= 0 and <= 1 !",call.=FALSE)
  }else if(name=="minRegCor") {
    if( !(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>1 || para<0)
        stop("'minRegCor' should be a numeric value >= 0 and <= 1 !",call.=FALSE)
  }else if(name=="prob") {
    if( !(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>1 || para<0)
      stop("'prob' should be a numeric value >= 0 and <= 1 !",call.=FALSE)
  }else if(name=="amapCutoff") {
    if(!is.null(para) && (!is.numeric(para) || length(para)!=1 || para>1 || para<0))
      stop("'amapCutoff' should be a numeric value >=0 and <=1!",call.=FALSE)
  }else if(name=="amapFilter"){
    if(!is.character(para) || length(para)!=1 || !(para %in% c("phyper","quantile","custom")))
      stop("'amapFilter' should be any one of 'quantile', 'phyper' and 'custom'!",call.=FALSE)
  }else if(name=="listOfRegulons") {
    if(!is.list(para))
      stop("'listOfRegulons' should be a list of gene sets!\n")
    if(  is.null(names(para)) || length(unique(names(para)))<length(names(para))  )
      stop("'listOfRegulons' should be a named list (unique names)!",call.=FALSE)
    junk<-lapply(para,function(reg){
      if(!is.character(reg))
        stop("'listOfRegulons' should include character vectors!",call.=FALSE)
    })
  }else if(name=="listOfRegulonsAndMode") {
    if(!is.list(para))
      stop("'listOfRegulonsAndMode' should be a list of gene sets!\n")
    if(  is.null(names(para)) || length(unique(names(para)))<length(names(para))  )
      stop("'listOfRegulonsAndMode' should be a named list (unique names)!",call.=FALSE)
    junk<-lapply(para,function(reg){
      if( !( is.numeric(reg) || is.integer(reg) ) )
        stop("'listOfRegulonsAndMode' should include named numeric or integer vectors!",call.=FALSE)
      if(is.null(names(reg)) || any(is.na(names(reg))) || any(names(reg)==""))
        stop("'listOfRegulonsAndMode' should include named vectors, without 'NA' or empty names!",call.=FALSE)
    })
  }else if(name=="glist") {
    if(!is.null(para)){
      if(!is.list(para))
        stop("'glist' should be a list of gene sets or regulons!\n")
      if(  is.null(names(para)) || any(duplicated(para)) )
        stop("'glist' should be a named list with unique names!",call.=FALSE)
      if(any(!unlist(lapply(para,class))=="character")){
        para<-lapply(para,function(gl){
          as.character(gl)
        })
      }
      return(para)
    }
  }else if(name=="pAdjustMethod") {
    if(!is.character(para) || length(para)!=1 || 
      !(para %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")))
      stop("'pAdjustMethod' should be any one of 'holm','hochberg','hommel','bonferroni','BH','BY' and 'fdr'!",call.=FALSE)
  }else if(name=="coverage") {
    if(!is.character(para) || length(para)!=1 || 
         !(para %in% c("tnet", "all")))
      stop("'coverage' should be any one of 'tnet' and 'all'!",call.=FALSE)
  }else if(name=="pValueCutoff") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>1 || para<0)
      stop("'pValueCutoff' should be an integer or numeric value >=0 and <=1  !",call.=FALSE)
  }else if(name=="miThreshold") {
    if(is.character(para)){
      if(!is.character(para) || length(para)!=1 || !(para %in% c("md","md.tf")))
        stop("'miThreshold' should be any one of 'md' and 'md.tf''!",call.=FALSE)
    }else {
      if( !(is.integer(para) || is.numeric(para)) || length(para)<1 || length(para)>2)
        stop("custom 'miThreshold' should be a numeric value, or a vector of length = 2!",call.=FALSE)
    }
  }else if(name=="consensus") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para>100 || para<1)
      stop("'consensus' should be an integer or numeric value <=100 and >=1 !",call.=FALSE)
  }else if(name=="pooledNullDistribution") {
    if(!is.singleLogical(para))
      stop("'pooledNullDistribution' should be a logical value!",call.=FALSE)
  }else if(name=="doSizeFilter") {
    if(!is.singleLogical(para))
      stop("'doSizeFilter' should be a logical value!",call.=FALSE)
  }else if(name=="targetContribution") {
    if(!is.singleLogical(para))
      stop("'targetContribution' should be a logical value!",call.=FALSE)
  }else if(name=="additionalData") {
    if(!is.singleLogical(para))
      stop("'additionalData' should be a logical value!",call.=FALSE)
  }else if(name=="scale") {
    if(!is.singleLogical(para))
      stop("'scale' should be a logical value!",call.=FALSE)
  }else if(name=="nPermutations") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || round(para,0)!=para)
      stop("'nPermutations' should be an integer >=1 !",call.=FALSE)
  }else if(name=="parChunks") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<2 || round(para,0)!=para)
      stop("'parChunks' should be an integer >=2 !",call.=FALSE)
  }else if(name=="nBootstraps") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || round(para,0)!=para)
      stop("'nBootstraps' should be an integer >=1 !",call.=FALSE)
  }else if(name=="minRegulonSize") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || round(para,0)!=para)
      stop("'minRegulonSize' should be an integer >=1 !",call.=FALSE) 
  }else if(name=="minPrunedSize") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || round(para,0)!=para)
        stop("'minPrunedSize' should be an integer >=1 !",call.=FALSE) 
  } else if(name=="activityMethod") {
  if(!is.character(para) || length(para)!=1 || !(para %in% c("gsea2", "GSEA2", "area3", "AREA3")))
    stop("'activityMethod' should be 'GSEA2' or 'AREA3'!",call.=FALSE)
  } else if(name=="minSize") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || round(para,0)!=para)
      stop("'minSize' should be an integer >=1 !",call.=FALSE) 
  }else if(name=="evse.minSize") {
    if(!(is.integer(para) || is.numeric(para)) || any(para<1) || any(round(para,0)!=para) )
      stop("'minSize' should be integer >=1 !",call.=FALSE)
    if( length(para)<1 || length(para)>2)
      stop("'minSize' should either be a single value or a vector of length = 2 !",call.=FALSE)
    if(length(para)==1)para<-c(para,para)
    return(para)
  }else if(name=="minIntersectSize") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<0 || para>100)
      stop("'minIntersectSize' should be an integer >=0 and <=100 !",call.=FALSE)
  }else if(name=="minAgreement") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1 || para>100)
      stop("'minAgreement' should be an integer >=1 and <=100 !",call.=FALSE)
  }else if(name=="statFilter") {
    if(!is.character(para) || length(para)!=1 || !(para %in% c("phyper","pchisq")))
      stop("'statFilter' should be any one of 'phyper' and 'pchisq'!",call.=FALSE)
  }else if(name=="statUniverse") {
    if(!is.character(para) || length(para)!=1 || !(para %in% c("all","tnet","reg")))
      stop("'statUniverse' should be any one of 'all' and 'tnet'!",call.=FALSE)
  }else if(name=="eps") {
    if(!is.na(para)){
      if(!is.numeric(para) || length(para)!=1 || para<0)
        stop("'eps' should be an numeric value >=0 or NA!",call.=FALSE)   
    }
  }else if(name=="exponent") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<1)
      stop("'exponent' should be an integer or numeric value >=1 !",call.=FALSE)
  }else if(name=="cvfilter") {
    if(!is.logical(para) || length(para)!=1)
      stop("'cvfilter' should be a logical value!",call.=FALSE)
  }else if(name=="verbose") {
    if(!is.logical(para) || length(para)!=1)
      stop("'verbose' should be a logical value!",call.=FALSE)
  }else if(name=="log") {
    if(!is.logical(para) || length(para)!=1)
      stop("'log' should be a logical value!",call.=FALSE)
  }else if(name=="iConstraint") {
    if(!is.logical(para) || length(para)!=1)
      stop("'iConstraint' should be a logical value!",call.=FALSE)
  }else if(name=="mask") {
    if(!is.logical(para) || length(para)!=1)
      stop("'mask' should be a logical value!",call.=FALSE)
  } else if(name=="decreasing") {
    if(!is.logical(para) || length(para)!=1)
      stop("'decreasing' should be a logical value!",call.=FALSE)
  }else if(name=="fineMapping") {
    if(!is.logical(para) || length(para)!=1)
      stop("'fineMapping' should be a logical value!",call.=FALSE)
  }else if(name=="boxcox") {
    if(!is.logical(para) || length(para)!=1)
      stop("'boxcox' should be a logical value!",call.=FALSE)
  }else if(name=="tnet") {
    if(!is.character(para) || length(para)!=1 || !(para %in% c("dpi", "ref")))
      stop("'tnet' should be any one of 'dpi' and 'ref'!",call.=FALSE)
  }else if(name=="stepFilter") {
    if(!is.logical(para) || length(para)!=1)
      stop("'stepFilter' should be a logical value!",call.=FALSE)
  } else if(name=="orderAbsValue") {
    if(!is.logical(para) || length(para)!=1)
      stop("'orderAbsValue' should be a logical value!",call.=FALSE)
  }else if(name=="phenotype") {
    if(!is.null(para)){
      if( !( is.numeric(para) || is.integer(para) ) || length(para)==0 )
        stop("'phenotype' should be a named numeric or integer vector with length >0 !",call.=FALSE)
      if(is.null(names(para)))
        stop("'phenotype' should be a named vector, without 'NA' or empty names!",call.=FALSE)
    }
  }else if(name=="hits") {
    if(!is.null(para)){
      if(!is.character(para) || length(para)==0 || !is.vector(para))
        stop("'hits' should be a character vector with length >0 !",call.=FALSE)      
    }
  }else if(name=="phenoIDs") {
    if(is.null(para)){
      return(para)
    }else {
      if( (!is.matrix(para) && !is.data.frame(para) ) || ncol(para)<2 ){
        stop("'phenoIDs' should be a dataframe (or a matrix of characters) with ncol >=2 !",call.=FALSE)
      }
      junk<-sapply(1:ncol(para),function(i){
        tp<-para[,i]
        if(is.list(tp))tp<-unlist(tp)
        if(i<=2 || is.factor(tp))tp<-as.character(tp)
        para[,i]<<-tp
      })
      rownames(para)<-para[,1]
      if( (any(is.na(para[,1])) || any(para[,1]=="") ) ){
        stop("Col 1 in 'phenoIDs' matrix should have no 'NA' or empty values!",call.=FALSE)      
      }
      para<-as.data.frame(para,stringsAsFactors=FALSE)
      return(para)
    }
  } else if (name == "tarPriorityMethod") {
      if(!(para %in% c("TC", "EC", "MI"))) {
          stop("`tarPriorityMethod` must be one of 'TC', 'EC' or 'MI'.")
      }
  } else if (name == "acitivityMethod") {
    if(!(para %in% c("GSEA2", "gsea2", "AREA3", "area3"))) {
        stop("`tarPriorityMethod` must be 'GSEA2' or 'AREA3'.")
    }
  } else if(name=="annotation.vse"){
    if( !is.data.frame(para) || ncol(para)<3 ){
      stop("'annotation' should be a dataframe with ncol>=3!",call.=FALSE)
    }
    colnames(para)<-toupper(colnames(para))
    cnames<-c("CHROM","START","END")
    if( !all(cnames%in%colnames(para)) ){
      tp1<-"\nCol names in the 'annotation' data do not match valid labels!\n"
      tp2<-"Please, revise the 'annotation' format:\n"
      tp3<-"Name col1: <CHROM> chromosome name\n"
      tp4<-"Name col2: <START> start position\n"
      tp5<-"Name col3: <END> end position\n"
      tp6<-"Name col4: <ID> any genomic id or name of the line\n"
      stop(tp1,tp2,tp3,tp4,tp5,tp6,call.=FALSE)
    }
    if(ncol(para)==3){
      para$ID<-paste(para$CHROM,para$START,para$END,sep="_")
    }else {
      para<-para[,1:4]
    }
    idx<-c(which(!colnames(para)%in%cnames),match(cnames,colnames(para)))
    para<-para[,idx,drop=FALSE]
    #---
    idlab<-colnames(para)[1]
    if(any(duplicated(para[,1]))){
      stop("Input data should have no duplicated annotation!",call.=FALSE)
    }
    colnames(para)[1]<-"ID"
    #---
    sapply(1:ncol(para),function(i){
      tp<-para[,i]
      if( any( is.na(tp) || any(tp=="") ) ){
        stop("'annotation' matrix should have no 'NA' or empty values!",call.=FALSE)      
      }
      if(is.list(tp))para[,i]<<-unlist(tp)
    })
    #---
    para$ID<-as.character(para$ID)
    para$CHROM<-as.character(para$CHROM)
    if(is.numeric(para$CHROM) || is.integer(para$CHROM)){
      para$CHROM <- paste("chr",para$CHROM,sep="")
    }
    chrs<-c(paste("chr",1:22,sep=""),"chrX")
    chrChecks<-!para$CHROM%in%chrs
    if( any(chrChecks) ){
      n<-sum(chrChecks)/length(chrChecks)
      tp1<-paste("chromosome names in 'annotation' should be listed in ", sep="")
      tp2<-"[chr1, chr2, chr3, ..., chr22, chrX]!"
      if(n>0.95){
        stop(tp1,tp2,call.=FALSE)
      }else {
        nonvalid <- unique(para$CHROM[chrChecks])
        if(length(nonvalid)>1) nonvalid <- c(nonvalid[1],"...",nonvalid[length(nonvalid)])
        nonvalid <- paste(nonvalid,collapse=", ")
        tp3<-"\n...the following values were removed: "
        warning(tp1,tp2,tp3,nonvalid,call.=FALSE)
        para<-para[!chrChecks,]
      }
    }
    #---
    b1<-is.integer(para$START) || is.numeric(para$START)
    b2<-is.integer(para$END) || is.numeric(para$END)
    if( !(b1 && b2) ){
      stop("Chromosome start/end positions in 'annotation' should be integer or numeric vectors!",call.=FALSE)
    }
    rownames(para)<-para$ID
    return(para)
  } else if(name=="annotation.evse"){
    if( !is.data.frame(para) || ncol(para)!=4 ){
      stop("'annotation' should be a dataframe with ncol=4 !",call.=FALSE)
    }
    colnames(para)<-toupper(colnames(para))
    cnames<-c("CHROM","START","END")
    if( !all(cnames%in%colnames(para)) ){
      tp1<-"\nCol names in the 'annotation' data do not match valid labels!\n"
      tp2<-"Please, revise the 'annotation' format:\n"
      tp3<-"Name col1: <CHROM> chromosome name\n"
      tp4<-"Name col2: <START> start position\n"
      tp5<-"Name col3: <END> end position\n"
      tp6<-"Name col4: <ID> any genomic id or name of the line\n"
      stop(tp1,tp2,tp3,tp4,tp5,tp6,call.=FALSE)
    }
    idx<-c(which(!colnames(para)%in%cnames),match(cnames,colnames(para)))
    para<-para[,idx,drop=FALSE]
    #---
    idlab<-colnames(para)[1]
    if(any(duplicated(para[,1]))){
      stop(paste("Col '",idlab,"' in the 'annotation' data should have non-duplicated IDs!",sep=""),call.=FALSE)
    }
    colnames(para)[1]<-"ID"
    #---
    sapply(1:ncol(para),function(i){
      tp<-para[,i]
      if( any( is.na(tp) || any(tp=="") ) ){
        stop("'annotation' matrix should have no 'NA' or empty values!",call.=FALSE)      
      }
      if(is.list(tp))para[,i]<<-unlist(tp)
    })
    #---
    para$ID<-as.character(para$ID)
    para$CHROM<-as.character(para$CHROM)
    if(is.numeric(para$CHROM) || is.integer(para$CHROM)){
      para$CHROM <- paste("chr",para$CHROM,sep="")
    }
    chrs<-c(paste("chr",1:22,sep=""),"chrX")
    chrChecks<-!para$CHROM%in%chrs
    if( any(chrChecks) ){
      n<-sum(chrChecks)/length(chrChecks)
      tp1<-paste("chromosome values in 'annotation' should be listed in ", sep="")
      tp2<-"[chr1, chr2, chr3, ..., chr22, chrX]!"
      if(n>0.95){
        stop(tp1,tp2,call.=FALSE)
      }else {
        nonvalid <- unique(para$CHROM[chrChecks])
        if(length(nonvalid)>1) nonvalid <- c(nonvalid[1],"...",nonvalid[length(nonvalid)])
        nonvalid <- paste(nonvalid,collapse=", ")
        tp3<-"\n...the following values were removed: "
        warning(tp1,tp2,tp3,nonvalid,call.=FALSE)
        para<-para[!chrChecks,]
      }
    }
    #---
    b1<-is.integer(para$START) || is.numeric(para$START)
    b2<-is.integer(para$END) || is.numeric(para$END)
    if( !(b1 && b2) ){
      stop("Chromosome start/end positions in 'annotation' should be integer or numeric vectors!",call.=FALSE)
    }
    rownames(para)<-para$ID
    return(para)
  } else if(name=="eqtls"){
    if( !is.data.frame(para) || ncol(para)<2 ){
      stop("'eqtls' should be a dataframe with ncol>=2 !",call.=FALSE)
    }
    colnames(para)<-toupper(colnames(para))
    if(!all(c("RSID","GENEID") %in% colnames(para))){
      stop("'eqtls' should include two columns named <RSID> and <GENEID> !",call.=FALSE)
    }
    para$RSID <- as.character(para$RSID)
    para$GENEID <- as.character(para$GENEID)
    return(para)
  } else if(name=="rowAnnotation"){
    if(is.null(para)){
      return(para)
    }else {
      if( (!is.matrix(para) && !is.data.frame(para) ) || ncol(para)<2){
        stop("'rowAnnotation' should be a data frame (or a matrix of characters) with ncol >=2 !",call.=FALSE)
      }
      if(is.matrix(para))para<-data.frame(para,stringsAsFactors=FALSE, check.names=FALSE)
      junk<-sapply(1:ncol(para),function(i){
        tp<-para[,i]
        if(is.list(tp))tp<-unlist(tp)
        tp<-as.character(tp)
        para[,i]<<-tp
      })
      #check colnames
      if(is.null(colnames(para))){
        stop("'rowAnnotation' must have colnames!")
      }
      colnames(para)<-toupper(colnames(para))
      if(any(duplicated(colnames(para)))){
        stop("'rowAnnotation' should have unique colnames (not case sensitive)!")
      }
      #check rownames
      if(is.null(rownames(para))){
        if(any(is.na(para[,1])) || any(para[,1]=="")){
          stop("Col 1 in 'rowAnnotation' should have no NA or empty value!",call.=FALSE)
        }
        if(any(duplicated(para[,1])))
          stop("Col 1 in 'rowAnnotation' should have unique ids!",call.=FALSE)
        rownames(para) <- para[,1]
      }else {
        if(!identical(rownames(para),para[,1])){
          para <- cbind(ID=rownames(para), para)
        }
      }
      return(para)
    }
  } else if(name=="colAnnotation"){
    if(is.null(para)){
      return(para)
    }else {
      if((!is.matrix(para) && !is.data.frame(para)) || ncol(para)<2){
        stop("'colAnnotation' should be a data frame (or a matrix of characters) with ncol >=2 !",call.=FALSE)
      }
      if(is.matrix(para)) {
        para <- data.frame(para,stringsAsFactors=FALSE, check.names=FALSE)
      }
      for (i in 1:ncol(para)) {
        column <- para[,i]
        if (is.list(column)) {
          column <- sapply(column, paste, collapse = ", ")
          para[,i] <- column
        }
      }
      #check colnames
      if(is.null(colnames(para))){
        stop("'colAnnotation' must have colnames!")
      }
      if(any(duplicated(colnames(para)))){
        stop("'colAnnotation' should have unique col names!")
      }
      #check rownames
      if(is.null(rownames(para))){
        if(any(is.na(para[,1])) || any(para[,1]=="")){
          stop("Col 1 in 'colAnnotation' should have no NA or empty value!",call.=FALSE)
        }
        if(any(duplicated(para[,1])))
          stop("Col 1 in 'colAnnotation' should have unique ids!",call.=FALSE)
        rownames(para) <- para[,1]
      }else {
        if(!identical(rownames(para),para[,1])){
          para <- cbind(ID=rownames(para), para)
        }
      }
      return(para)
    }
  } else if(name=="duplicateRemoverMethod") {
    if(!is.character(para) || length(para) != 1 || !(para %in% c("max","min","average")))
      stop("'duplicateRemoverMethod' should be one of the following character strings: 'max', 'min' and 'average'",call.=FALSE)
  } else if(name=="gs") {
    if(!is.character(para) || length(para)==0 || any(is.na(para)) || any(para==""))
      stop("'geneSet/regulon' should be a character vector with length > 0, without 'NA' or empty names!",call.=FALSE)
  } else if(name=="universe") {
    if(!is.character(para) || any(is.na(para)) || any(para=="")) 
      stop("'universe' should be a character vector without 'NA' or empty values!",call.=FALSE)
  }  else if(name=="globalAdjustment") {
    if(!is.logical(para) || length(para)!=1)
      stop("'globalAdjustment' should be a logical value!",call.=FALSE)
  } else if(name=="filepath"){
    if(!is.character(para) || length(para)!=1)
      stop("'filepath' should be a single character!",call.=FALSE)
  } else if(name=="file"){
    if(!is.character(para) || length(para)!=1)
      stop("'file' should be a single character!",call.=FALSE)
  } else if(name=="ylimPanels") {
    if(!(is.numeric(para) || is.integer(para)) || length(para)!=4 )
      stop("'ylimPanels' should be a numeric vector with length =4 !",call.=FALSE)
  } else if(name=="heightPanels") {
      if(!(is.numeric(para) || is.integer(para)) || length(para)!=3 )
      stop("'heightPanels' should be a numeric vector with length =3 !",call.=FALSE)
  } else if(name=="bxseq") {
    if(!(is.numeric(para) || is.integer(para)) || length(para)<2 )
      stop("'bxseq' should be a numeric vector with length >=2 !",call.=FALSE)
  } else if(name=="width") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<0 )
      stop("'width' should be an integer or numeric value >0 !",call.=FALSE)
  } else if(name=="height") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<0 )
      stop("'height' should be an integer or numeric value >0 !",call.=FALSE)
  } else if(name=="maxy") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<0)
      stop("'maxy' should be an integer or numeric value >0 !",call.=FALSE)
  } else if(name=="rmargin") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<0 )
      stop("'rmargin' should be an integer or numeric value >0 !",call.=FALSE)
  } else if(name=="bxsp") {
    if(!(is.integer(para) || is.numeric(para)) || length(para)!=1 || para<0)
      stop("'bxsp' should be an integer or numeric value >0 !",call.=FALSE)
  } else if(name=="ylabPanels"){
    if(!is.character(para) || length(para)!=3)
      stop("'ylabPanels' should be a character vector with length=3 !",call.=FALSE)
  } else if(name=="xlab"){
    if(!is.character(para) || length(para)!=1)
      stop("'xlab' should be a character value with length=1 !",call.=FALSE)
  } else if(name=="ylab"){
    if(!is.character(para) || length(para)!=1)
      stop("'ylab' should be a character value with length=1 !",call.=FALSE)
  } else if(name=="lab"){
    if(!is.character(para) || length(para)!=1)
      stop("'lab' should be a character value with length=1 !",call.=FALSE)
  } else if(name=="tlt"){
    if(!is.character(para) || length(para)!=1)
      stop("'tlt' should be a character value with length=1 !",call.=FALSE)
  } else if(name=="fname") {
    if(!is.character(para) || length(para)!=1)
      stop("'fname' should be a character value with length=1!",call.=FALSE)
  } else if(name=="lab"){
    if(!is.character(para) || length(para)!=1)
      stop("'lab' should be a character value with length=1 !",call.=FALSE)
  } else if(name=="labPheno") {
    if(!is.character(para) || length(para)!=1)
      stop("'labPheno' should be a character value with length=1 !",call.=FALSE)
  } else if(name=="idkey") {
    if(!is.null(para)){
      if(!is.character(para) || length(para)!=1)
        stop("'idkey' should be a character value with length=1 !",call.=FALSE)
    }
  } else if (name=="plotpdf") {
    if(!is.logical(para) || length(para)!=1)
      stop("'plotpdf' should be a logical value with length=1 !", call.=FALSE)
  } else {
    stop("...<",name,"> check is missing!!",call. = FALSE)
  }
}
##------------------------------------------------------------------------------
is.singleNumber <- function(para){
  (is.integer(para) || is.numeric(para)) && length(para)==1L && !is.na(para)
}
is.singleInteger <- function(para){
  lg <- (is.integer(para) || is.numeric(para)) && length(para)==1L && 
    !is.na(para)
  if(lg) lg <- (para / ceiling(para)) == 1
  return(lg)
}
is.singleString <- function(para){
  is.character(para) && length(para) == 1L && !is.na(para)
}
is.singleLogical <- function(para){
  is.logical(para) && length(para) == 1L && !is.na(para)
}
all.binaryValues <- function(para){
  all( para %in% c(0, 1, NA) )
}
all.integerValues <- function(para){
  lg <- ( all(is.integer(para)) || all(is.numeric(para)) ) && 
    !any(is.na(para))
  if(lg) lg <- all ( (para / ceiling(para)) == 1 )
  return(lg)
}
all.characterValues <- function(para){
  all(is.character(para)) && !any(is.na(para))
}
is.color <- function(x){
  res <- try(col2rgb(x),silent=TRUE)
  return(!"try-error"%in%class(res))
}
