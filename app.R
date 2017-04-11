library(shiny)
library(Biobase)
library(data.table)
library(rjson)
library(DT)
library(CBMRtools)

#get top and bottom matches for given eset of connectivity scores
summarize_eset<-function(mat, 
  summarize.func = c("mean", "median", "max", "min"),
  do.scorecutoff = TRUE, scorecutoff = c(-0.6, 0.6), 
  do.nmarkers = TRUE, nmarkers = c(100, 100)
  ){

  summarize.func<- match.arg(summarize.func)
  x<-apply(mat, 1, match.fun(summarize.func))
  x<-as.numeric(x)
  n<-length(x)
  
  if(do.nmarkers){
    ord<-order(x, decreasing = TRUE)
    x.ind.nmarkers<-c(ord[1:nmarkers[1]], ord[(n-nmarkers[2]+1):n])
  } else
    x.ind.nmarkers<-1:n

  if(do.scorecutoff)
    #TODO: rank by score here too
    x.ind.scorecutoff<-which(x > scorecutoff[2] | x < scorecutoff[1])
  else
    x.ind.scorecutoff<-1:n

  inds<-intersect(x.ind.nmarkers, x.ind.scorecutoff)
  inds<-inds[order(x[inds], decreasing = TRUE)]
  return(list(inds = inds, scores = x[inds]))
} 

#wrapper to retrieve top hits from an eset
get_connectivity<-function(input, tab, 
  pdat,fdat, 
  summarize.func = c("mean", "median", "max", "min"),
  do.scorecutoff = TRUE, scorecutoff = c(-0.6, 0.6), 
  do.nmarkers = TRUE, nmarkers = c(100, 100)){
  
  i<-get_BUID(input, tab)
  mat<-get_mat(i, datadir)
  #subset mat to rows that are in fdat
  inds.keep<-match(rownames(fdat), rownames(mat))
  mat<-mat[inds.keep, ]
  mat.cols<-round(pdat[match(colnames(mat), rownames(pdat)),"pert_dose"])
  colnames(mat)<-paste("WTCS_dose", mat.cols, sep = "")

  res<-summarize_eset(mat, summarize.func, do.scorecutoff, scorecutoff,
    do.nmarkers, nmarkers)

  res.ind<-res$inds
  res.scores<-res$scores
  tab<-cbind(fdat[res.ind,], summaryscore=res.scores, mat[res.ind,])
  tab<-data.table(tab)
}

get_chemical_description<-function(input, pdat, 
  datadir, tab,
  cols.keep = c("BUID", "Chemical.name", "CAS", "Broad_external_Id", "carc_liv_final")){
  i<-get_BUID(input, tab)
  mat<-get_mat(i, datadir)
  res<-pdat[which(rownames(pdat) %in% colnames(mat)),cols.keep]
  return(data.table(res[1,]))
}

get_mat<-function(i, datadir){
  return(readRDS(paste(datadir, "/", i, ".RDS", sep = "")))
}

subset_fdat<-function(fdat, keyword, x.split){
  if(keyword %in% "all")
    return(fdat)
  else 
    inds<-grep(keyword, x.split)
  return(fdat[inds,])
}

get_cellines<-function(fdat){
  x<-as.character(fdat$id)
  x.split<-unlist(lapply(x, function(i){strsplit(i, split = "_")[[1]][2]}))
  unique(as.character(x.split))
}

get_BUID<-function(input, tab){
  as.character(tab[which(apply(tab, 1, function(i) any(i %in% input)))[1], "BUID"])
}

get_ids_pdat<-function(pdat, cols = c("Chemical.name", "CAS", "BUID"), col.unique = "BUID",
  val.ignore = c("", " ", NA, "NA", "NOCAS")){
  tab<-unique(pdat[, cols])
  res<-lapply(cols, function(i){
  x<-as.character(tab[,i])
  x.uniq<- setdiff(unique(x), union(x[duplicated(x)], val.ignore))
  })
  names(res)<-cols
  return(res)
}

summarize_gsproj<-function(eset, order.col = "median"){
  res<-apply(exprs(eset), 1, summary)
  res<-t(res)
  colnames(res)<-c("min","Q1", "median", "mean","Q3","max")
  
  mat<-exprs(eset)
  colnames(mat)<-paste("gsscore",round(pData(eset)$pert_dose), sep = "")

  res<-cbind(rowIDs = rownames(res), fData(eset), summaryscore = res[,order.col], mat)
  res<-res[order(res$summaryscore, decreasing = TRUE),, drop = FALSE]
  res<-data.table(res)
  return(res)
}

get_gsproj<-function(input, gslist, tab, gsname, gsmethod){
  i<-get_BUID(input, tab)
  res<-gslist[[gsname]][[gsmethod]]
  res<-res[, res$BUID %in% i]
  return(res)
}

get_gsproj_list<-function(gsnames, gsmethods, gsdir){
  res<-lapply(names(gsnames), function(i){
    gsnameval<-gsnames[[i]]
    res2<-lapply(gsmethods, function(j){
      readRDS(paste(gsdir, "/", gsnameval, "_",j, ".RDS", sep = ""))
      })
    names(res2)<-gsmethods
    return(res2)
    })
  names(res)<-names(gsnames)
  return(res)
}

get_de<-function(input, tab, eset, landmark = FALSE, do.scorecutoff = TRUE, scorecutoff = c(-2, 2), 
  do.nmarkers = TRUE, nmarkers = c(100, 100),
  summarize.func = c("mean", "median", "max", "min")){
  i<-get_BUID(input, tab)
  eset<-eset[, eset$BUID %in% i]
  
  if(landmark)
    eset<-eset[fData(eset)$pr_is_lmark %in% "Y",]
  
  mat<-exprs(eset)
  fdat<-fData(eset)[, c("id", "pr_gene_symbol", "pr_is_lmark")]
  pdat<-pData(eset)

  colnames(mat)<-paste("modz_dose", round(pdat$pert_dose), sep = "")
  res<-summarize_eset(mat, summarize.func, do.scorecutoff, scorecutoff,
    do.nmarkers, nmarkers)

  res.ind<-res$inds
  res.scores<-res$scores
  tab<-cbind(fdat[res.ind,, drop = FALSE], summaryscore=res.scores, mat[res.ind,, drop = FALSE])
  tab<-data.table(tab)
}

filtereset<-function(eset, filteropt, tab){
  if (filteropt %in% "all")
    eset <- eset
  else if (filteropt %in% "carc_pos") 
    eset<-eset[, eset$carc_liv_final %in% "POSITIVE"]
  else if (filteropt %in% "carc_neg") 
    eset<-eset[, eset$carc_liv_final %in% "NEGATIVE"]
  else if (filteropt %in% "carc_pos_geno_pos") 
    eset<-eset[, eset$carc_liv_final %in% "POSITIVE" & eset$Genotoxicity %in% "POSITIVE"]
  else if (filteropt %in% "carc_pos_geno_neg") 
    eset<-eset[, eset$carc_liv_final %in% "POSITIVE" & eset$Genotoxicity %in% "NEGATIVE"]
  else if (filteropt %in% "carc_neg_geno_pos") 
    eset<-eset[, eset$carc_liv_final %in% "NEGATIVE" & eset$Genotoxicity %in% "POSITIVE"]
  else if (filteropt %in% "carc_neg_geno_neg") 
    eset<-eset[, eset$carc_liv_final %in% "NEGATIVE" & eset$Genotoxicity %in% "NEGATIVE"]
  else if (filteropt %in% "D.Sherr")
    eset<-eset[, grep("collaborator_David_Sherr", eset$source)]
  else if (filteropt %in% "ss4")
    eset<-eset[, eset$distil_ss > 4]
  else if (filteropt %in% "ss5")
    eset<-eset[, eset$distil_ss > 5]
  else if (filteropt %in% "ss6")
    eset<-eset[, eset$distil_ss > 6]
  else if (filteropt %in% "q75rep_0.2")
    eset<-eset[, eset$q75rep > 0.2]
  else if (filteropt %in% "q75rep_0.3")
    eset<-eset[, eset$q75rep > 0.3]
  else 
    eset<-eset[, eset$BUID %in% get_BUID(filteropt, tab)]
  return(eset)
}

subset_names<-function(x,n){
  if(nchar(x)> n) return(paste(substr(x, 1, n), "...", sep = ""))
  else return(x)
}

##load data dirs
dirs<-fromJSON(file = "datadirs.json")

##load data for chemical annotation tab
chemannot<-readRDS(file = dirs$chemannotation_filename)

##load data for Differential Expression tab
deeset<-readRDS(dirs$diffexp_filename)

##load data for GeneSetEnrichment tab
gsnames<-list(#c2.reactome = "gsscores_c2.cp.reactome.v5.0.annotated", 
  c2.reactome = "gsscores_c2.cp.reactome.v5.0", 
  hallmark = "gsscores_h.all.v5.0")
gsmethods<-c("gsproj", "gsva", "ssgsea", "zscore")
gsdir<-dirs$genesetenrich_dir
gslist<-get_gsproj_list(gsnames, gsmethods, gsdir)
gssort<-c("min","Q1", "median", "mean","Q3","max")

##load data for Connectivity tab
dir_connectivity<-dirs$connectivity_dir
datadir<-paste(dir_connectivity, "/expression", sep = "")
pdat<-readRDS(file = paste(dir_connectivity, "/", grep("pdat", list.files(dir_connectivity), value = T)[1], sep = ""))
fdat<-readRDS(file = paste(dir_connectivity, "/", grep("fdat", list.files(dir_connectivity), value = T)[1], sep = ""))

ids<-list.files(datadir)
ids<-gsub(".RDS", "", ids)
pdat<-pdat[pdat$BUID %in% ids,]

cols<-c("BUID", "Chemical.name", "CAS")
tab<-unique(pdat[, cols])
chemicals<-get_ids_pdat(pdat)


x<-as.character(fdat$id)
#get cell line associated with fdat$id
x.split<-unlist(lapply(x, function(i){strsplit(i, split = "_")[[1]][2]})) 

#unique cell lines
cellines<-c("all",get_cellines(fdat))
summarizefuncs<-c("max", "median", "mean", "min")

#directory of heatmap data files
heatmapdir<-dirs$heatmap_dir
heatmapfiles<-gsub(".RDS", "",list.files(heatmapdir))

filteropts<-c(list(premade_sets = c("all", "carc_pos", "carc_neg", "carc_pos_geno_pos",
  "carc_pos_geno_neg", "carc_neg_geno_pos", "carc_neg_geno_pos", "ss4", "ss5", "ss6", "q75rep_0.2",
  "q75rep_0.3", "D.Sherr")), chemicals)
##define app
app<-shinyApp(

ui = shinyUI(navbarPage("CRCGN Portal",

    tabPanel("About",
      titlePanel("CRCGN Liver Portal"),
      fluidRow(
        column(6, offset = 0.1,
          helpText("A portal for retrieving and visualising data related to the CRCGN liver project"))
      ),
      fluidRow(   
        column(6, offset = 0.1,
          helpText("Dataset details: 332 chemicals (128 liver carcinogens, 168 non-carcinogens, 36 others)"))
      )
    ),

    tabPanel("Chemical Annotation",
      DT::dataTableOutput("chemannot_result")
      ),

    tabPanel("Differential Expression",
      fluidPage(
        #dropdown menus
        fluidRow(
        column(3, selectInput("chemical_de", "CRCGN Chemical:", chemicals)),
        column(2, checkboxInput("landmark_de", "Landmark only", value =FALSE)),
        column(2, selectInput("summarizefunc_de", "Summarization:",summarizefuncs,
          selected = "median")),
        
        column(1,checkboxGroupInput("filterbyinput_de", "Filter by:",
                     c("score" = "score",
                       "number" = "number"),
                     selected = c("score", "number"))),
        column(2,sliderInput("range_de", "score threshold", min = -10, max = 10, 
          value = c(-2,2), step = 0.01)),
        column(1,sliderInput("numberthresleft_de", "Num +",
          min = 0, max = 1000, value = 10, ticks = FALSE, step = 10)),
        column(1,sliderInput("numberthresright_de", "Num -",
          min = 0, max = 1000, value = 10, ticks = FALSE, step = 10))
        ),
        # Table returned bshowing description for query chemical
        dataTableOutput("chemical_description_de"),
        tags$style(type="text/css", '#chemical_description_de tfoot {display:none;}'),

        ##insert empty space
        fluidRow(column(width = 1, offset = 0, style='padding:10px;')),

        # Table returned showing de markers
        dataTableOutput("result_de")
      )
    ),
		
		tabPanel("Gene set enrichment",
      fluidPage(
        fluidRow(
          column(3, selectInput("chemical_gs", "CRCGN Chemical:", chemicals)),
          column(2, selectInput("gsname_gs", "Gene set name:", names(gsnames))),
          column(2, selectInput("gsmethod_gs", "Projection method:", gsmethods)),
          column(2, selectInput("summarize_gs", "Sort by:", gssort,selected = "median"))
          )
        ),
      DT::dataTableOutput("gsproj_result")
    ),
    
    tabPanel("Connectivity",
      fluidPage(
        #dropdown menus
        fluidRow(
        column(3, selectInput("chemical", "CRCGN Chemical:", chemicals)),
        column(2, selectInput("celline", "CMAP cell line:", cellines)),
        column(2, selectInput("summarizefunc", "Summarization:",summarizefuncs,
          selected = "median")),
        
        column(1,checkboxGroupInput("filterbyinput", "Filter by:",
                     c("score" = "score",
                       "number" = "number"),
                     selected = c("score", "number"))),
        column(2,sliderInput("range", "score threshold", min = -1, max = 1, 
          value = c(-0.3,0.3), step = 0.01)),
        column(1,sliderInput("numberthresleft", "Num +",
          min = 0, max = 1000, value = 50, ticks = FALSE, step = 10)),
        column(1,sliderInput("numberthresright", "Num -",
          min = 0, max = 1000, value = 50, ticks = FALSE, step = 10))
        ),
        # Table returned bshowing description for query chemical
        dataTableOutput("chemical_description"),
        tags$style(type="text/css", '#chemical_description tfoot {display:none;}'),

        ##insert empty space
        fluidRow(column(width = 1, offset = 0, style='padding:10px;')),

        # Table returned showing connectivity scores for best matches
        dataTableOutput("result")
      )
    ),
    tabPanel("Heatmaps",
      fluidPage(
        fluidRow(
        column(3, selectInput("heatmapfile", "Dataset", heatmapfiles, 
          selected = "gsscores_h.all.v5.0_gsproj")),
        column(3, selectInput("filteropt", "Filter", filteropts, 
          selected = "2,6-Dinitrotoluene"))),
        plotOutput("heatmap_result")

      )
    )
	)),

server = shinyServer(function(input, output) {

  output$chemannot_result<-DT::renderDataTable({
    data.table(chemannot)
    },
    extensions = 'Buttons',
    server = FALSE,
    options = list (dom = 'T<"clear">Blfrtip',
      autoWidth = FALSE,
      columnDefs = list(list(width = '100px', targets = list(2),
          render = JS(
            "function(data, type, row, meta) {",
            "console.log(data);console.log(type);console.log(row);console.log(meta)",
            "return type === 'display' && data.length > 40 ?",
            "'<span title=\"' + data + '\">' + data.substr(0, 40) + '...</span>' : data;",
            "}")), 
      list(targets = list(1), visible = FALSE)
      ),
      deferRender=FALSE,
      scrollX=TRUE,scrollY=400,
      scrollCollapse=TRUE,
      pageLength = 1000, lengthMenu = c(10,50,100,200,1000),
      buttons=c('copy','csv','print')))

  output$chemical_description<-DT::renderDataTable({
    get_chemical_description(input$chemical, pdat, datadir, tab)
    }, options = list(dom = ''))

  output$result <- DT::renderDataTable({
      get_connectivity(input$chemical, tab, 
        pdat,
        subset_fdat(fdat, input$celline, x.split), 
        input$summarizefunc, "score" %in% input$filterbyinput,
        c(input$range[1], input$range[2]), 
        "number" %in% input$filterbyinput, 
        c(input$numberthresleft, input$numberthresright)
        )
    }, 
    extensions = 'Buttons', 
    server = FALSE,
    options = list(dom = 'T<"clear">Blfrtip', 
    pageLength = 10,
    deferRender=FALSE, 
    buttons=c('copy','csv','print')))

  output$gsproj_result <- DT::renderDataTable({
      eset<-get_gsproj(input$chemical_gs, gslist, tab, input$gsname_gs, input$gsmethod_gs)
      res<-summarize_gsproj(eset, order.col = input$summarize_gs)
      return(res)
    }, 
    extensions = 'Buttons',
    server = FALSE,
    options = list (dom = 'T<"clear">Blfrtip',
      autoWidth = FALSE,
      columnDefs = list(list(width = '100px', targets = list(2),
          render = JS(
            "function(data, type, row, meta) {",
            "console.log(data);console.log(type);console.log(row);console.log(meta)",
            "return type === 'display' && data.length > 40 ?",
            "'<span title=\"' + data + '\">' + data.substr(0, 40) + '...</span>' : data;",
            "}")), 
      list(targets = list(1), visible = FALSE)
      ),
      deferRender=FALSE,
      scrollX=TRUE,scrollY=400,
      scrollCollapse=TRUE,
      pageLength = 100, lengthMenu = c(10,50,100,200),
      buttons=c('copy','csv','print')
    )

  )

  output$chemical_description_de<-DT::renderDataTable({
    get_chemical_description(input$chemical_de, pdat, datadir, tab)
    }, options = list(dom = ''))

  output$result_de <- DT::renderDataTable({
      get_de(input$chemical_de, tab, 
        eset = deeset, 
        landmark = input$landmark_de, 
        do.scorecutoff = "score" %in% input$filterbyinput_de, 
        scorecutoff = c(input$range_de[1], input$range_de[2]), 
        do.nmarkers = "number" %in% input$filterbyinput_de, 
        nmarkers = c(input$numberthresleft_de, input$numberthresright_de),
        summarize.func = input$summarizefunc_de)

    }, 
    extensions = 'Buttons', 
    server = FALSE,
    options = list(dom = 'T<"clear">Blfrtip', 
    pageLength = 10,
    deferRender=FALSE, 
    buttons=c('copy','csv','print')))

  output$heatmap_result<-renderPlot({
      outfile<-paste(heatmapdir, "/", input$heatmapfile, ".RDS", sep = "")
      eset<-readRDS(outfile)
      eset<-filtereset(eset, input$filteropt, tab)

    
      cns<-as.character(sapply(as.character(eset$Chemical.name), 
        function(i) subset_names(i,30)))
      cns<-make.unique(cns) #in case truncation produces non-unique ids
      colnames(eset)<-paste(cns, "_", eset$pert_dose_RANK, sep = "")
      hc<-clust_eset(eset)
      col_legend<-list(carc_liv_final = list(col_breaks = c("POSITIVE", "NEGATIVE", ""), 
        col_values = sapply(c("orange", "green", "white"), to.hex),
        col_labels = c("POSITIVE", "NEGATIVE", "")),
      Genotoxicity = list(col_breaks = c("POSITIVE", "NEGATIVE", ""), 
        col_values = sapply(c("orange", "green", "white"), to.hex),
        col_labels = c("POSITIVE", "NEGATIVE", "")),
      pert_dose_RANK = list(col_breaks = 1:6, 
        col_values = rev(gray.colors(6)),
        col_labels = 1:6))
      hmcolors<-function(...) scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, ...)
      #xsize <- 50/ncol(eset)
      ncols<-ncol(eset)
      xsize<-3
      if (ncols<500) 
        xsize<-4
      if (ncols<100)
        xsize<-6
      if(ncols<10)
        xsize<-10
      if(ncols > 1000)
        xsize<-0 

      p<-ggheat.continuous.single(eset = eset, 
        hc = hc$hc, 
        hr = hc$hr, 
        hmcolors = hmcolors,
        hmtitle = "geneset score",
        col_lab = c("carc_liv_final", "Genotoxicity", "pert_dose_RANK"), 
        col_legend = col_legend,
        ylabstr = "",
        fout = NA, 
        p.heights = c(1.5, 0.5, 5),
        xsize = xsize,
        ysize = 3, 
        ysizelab = 7,
        xright = 0.18)
      return(p)
  })

})

)

