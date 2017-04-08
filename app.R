library(shiny)
library(Biobase)
library(data.table)
library(rjson)

#get top and bottom matches for given eset of connectivity scores
parse_conc_eset<-function(mat, 
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
  return(list(inds = inds, scores = x[inds]))
} 

#wrapper to retrieve top hits from an eset
get_result<-function(input, tab, fdat, summarize.func = c("mean", "median", "max", "min"),
  do.scorecutoff = TRUE, scorecutoff = c(-0.6, 0.6), 
  do.nmarkers = TRUE, nmarkers = c(100, 100)){
  
  i<-get_BUID(input, tab)
  mat<-get_mat(i, datadir)
  #subset mat to rows that are in fdat
  inds.keep<-match(rownames(fdat), rownames(mat))
  mat<-mat[inds.keep, ]

  res<-parse_conc_eset(mat, summarize.func, do.scorecutoff, scorecutoff,
    do.nmarkers, nmarkers)

  res.ind<-res$inds
  res.scores<-res$scores
  tab<-cbind(fdat[res.ind,], summaryscore=res.scores)
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

query<-"CRCGN_liver"
database<-"cmap"
titletext<-paste("Connectivity Summary: ",  "query: ", query, ", database: ", database, sep = "")


dirs<-fromJSON(file = "datadirs.json")
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

app<-shinyApp(

ui = shinyUI(navbarPage("CRCGN Portal",

		tabPanel("Connectivity",
			fluidPage(
			  # Application title
			  titlePanel(titletext),
        #dropdown menus
			  fluidRow(
        column(3, selectInput("chemical", "CRCGN Chemical:", chemicals)),
        column(2, selectInput("celline", "CMAP cell line:", cellines)),
        column(2, selectInput("summarizefunc", "Summarization:",summarizefuncs)),
        
        column(1,checkboxGroupInput("filterbyinput", "Filter by:",
                     c("score" = "score",
                       "number" = "number"),
                     selected = c("score", "number"))),
        column(2,sliderInput("range", "score threshold", min = -1, max = 1, value = c(-0.3,0.3), step = 0.01)),
        column(1,sliderInput("numberthresleft", "Num +",
          min = 0, max = 1000, value = 50, ticks = FALSE, step = 10)),
        column(1,sliderInput("numberthresright", "Num -",
          min = 0, max = 1000, value = 50, ticks = FALSE, step = 10))
        ),
        #textOutput("filterby"),
			  # Table returned bshowing description for query chemical
			  dataTableOutput("chemical_description"),
			  tags$style(type="text/css", '#chemical_description tfoot {display:none;}'),

			  # Table returned showing connectivity scores for best matches
			  dataTableOutput("result")
			)
		),
		
		tabPanel("Gene set enrichment",
			verbatimTextOutput("temp_message")),
		tabPanel("Differential Expression",
      verbatimTextOutput("temp_mesage2"))
		)),

server = shinyServer(function(input, output) {

  output$chemical_description<-renderDataTable({
    get_chemical_description(input$chemical, pdat, datadir, tab)
    }, options = list(dom = ''))

  output$result <- renderDataTable({
      get_result(input$chemical, tab, 
        subset_fdat(fdat, input$celline, x.split), 
        input$summarizefunc, "score" %in% input$filterbyinput,
        c(input$range[1], input$range[2]), 
        "number" %in% input$filterbyinput, 
        c(input$numberthresleft, input$numberthresright)
        )
    }, options = list(pageLength = 10))

  output$temp_message<-renderText({"Under construction" })
  output$temp_message2<-renderText({"Under construction" })

})

)

