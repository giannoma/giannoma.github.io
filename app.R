##
##      __________  _____
##     / ____/ __ \/ ___/
##    / /   / / / /\__ \
##   / /___/ /_/ /___/ /
##   \____/_____//____/
##
##  Creative Data Solutions
##  Vanderbilt University
##  https://cds.vanderbilt.edu
##
## Author(s): Matthew A. Cottam, Jean-Philippe Cartailler, Manuel A. Giannoni-Guzmán
## Date Created: 2023-06-27
## Date Last Updated: 2024-03-19
##
## ---------------------------

# Global libraries
library(shiny)
library(shinyjs)
library(shinythemes)
library(shinyhelper)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggforce)
library(knitr)
library(kableExtra)
library(colourpicker)
library(readr)
library(ggrepel)
library(stringr)
library(shinycssloaders)
library(polished)



# Password protection REMOVE ONCE PAPER IS ACCEPTED###
polished::polished_config(
  app_name = "RNA_seq_website",
  api_key = "z2WzC0IxQBKHe95aiVi2ICQ4B8czLw3Re2"
)

password_protect = FALSE # TRUE or FALSE

## -----------------------------------------------------------------------------
##
## CONFIGURATION
##
## -----------------------------------------------------------------------------

citation <- "https://doi.org/10.1101/2024.02.16.580759"

if(FALSE){
  # Testing, execute manually
  setwd("ShinyApp/")
}

load_gene_model <- "protein_coding" # "all" or "protein_coding"


options(spinner.type = 2, spinner.color = "#56B4E9", spinner.color.background = "#ffffff", spinner.size = 1)

color_green <- "#23a127"
color_red <- "#d73027"
color_blue <- "#4575b4"

# Global raw data
if(load_gene_model ==  "all"){

  # ALL GENES
  genemap <- read.csv("data/all_genes/genemap.csv")
  metadata <- read.csv("data/all_genes/metadata.csv")
  DCP_params <- read.csv("data/all_genes/DR_rhythm_params.csv")
  DCP_rhythm <- readRDS("data/all_genes/DCP_rhythm.rds")
  DEGs <- read.csv("data/all_genes/LONG__vs__SHORT-DESeq2-results-all-data.csv", row.names = 1)
  DEGs_normalized_counts_SHORT <- read.csv("data/all_genes/normalized_counts_all.SHORT.csv.gz", row.names = 1, check.names = FALSE)
  DEGs_normalized_counts_LONG <- read.csv("data/all_genes/normalized_counts_all.LONG.csv.gz", row.names = 1, check.names = FALSE)

}else if(load_gene_model ==  "protein_coding"){

  # PROTEIN_CODING only genes
  genemap <- read.csv("data/protein_coding/genemap.csv")
  metadata <- read.csv("data/protein_coding/metadata.csv")
  DCP_params <- read.csv("data/protein_coding/DR_rhythm_params.csv")

  DCP_rhythm <- readRDS("data/protein_coding/DCP_rhythm.rds")
  DEGs <- read.csv("data/protein_coding/LONG__vs__SHORT-DESeq2-results-all-data.csv", row.names = 1)
  DEGs_normalized_counts_SHORT <- read.csv("data/protein_coding/normalized_counts_all.SHORT.csv.gz", row.names = 1, check.names = FALSE)
  DEGs_normalized_counts_LONG <- read.csv("data/protein_coding/normalized_counts_all.LONG.csv.gz", row.names = 1, check.names = FALSE)

}else{
  
  stop("Cannot proceed, unknown gene model")
  
}




# -----------------------------------------------------------------------------
#
# Data remodeling
#
# -----------------------------------------------------------------------------


# Gene names available for selector
## 1. Prepare genemap

# genemap may have empty symbols, so let's substitute
genemap$mgi_symbol[genemap$mgi_symbol == ''] <- genemap$ensembl_gene_id

DCP_rhythm$rhythm.joint$ensemblid <- gsub("\\..*", "", DCP_rhythm$rhythm.joint$gname)
genelist <- genemap[genemap$ensembl_gene_id %in% DCP_rhythm$rhythm.joint$ensemblid,]$mgi_symbol # 23631 elements


# There are genes without symbols, so

# Sort gene list, so that letters come before numbers
genelist <- c(sort(genelist[!grepl("\\d", genelist)]), sort(genelist[grepl("\\d", genelist)]))
genelist <- genelist[nzchar(genelist)] # remove those without gene symbols, since empty string cannot be searched


## 2. Prepare raw values for DCP, SHORT and LONG
dcp_raw_SHORT <- DCP_rhythm$x1$rhythm
dcp_raw_LONG <- DCP_rhythm$x2$rhythm



# -----------------------------------------------------------------------------
#
# Shiny app
#
# -----------------------------------------------------------------------------


# Define UI for application that draws a histogram ----
ui <- fluidPage(

  useShinyjs(),

  # Set theme

  theme = shinytheme("yeti"),

  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "app.css")
  ),
  # Application title
  titlePanel("Circadian Photoperiod Seq: SCN Gene Rhythms and Expression Browser"),
  p("Please cite Cox et al., Journal of Biological Rhythms, 2024, Volume:pages"),
  #h4("Additional information..."), # todo change before launch

  # Sidebar Panel
  sidebarLayout(
    sidebarPanel(
      id = "well_light",
      width = 2,
      # Selector input box
      selectizeInput("gene",
                     "Select a gene:",
                     choices = genelist
                     #,selected = "Ciart"
                     ) %>%
        helper(type = "inline", size = "m", fade = TRUE,
               title = "How to use this app",
               content = "<b>Instructions:</b>
                 <p>Select a gene of interest using mouse gene symbols and then
                 click <strong>Update</strong>.</p>
                 <p><em>Optional:</em> Change the point size, line size, and
                 colors and then click <strong>Update</strong>.
                 To save high-resolution images of downloaded plots, click the
                 download buttons beneath each plot.</p>
                 <hr><b>Interpreting Plots:</b>
                 <p>The left plot displays the curve fit generated by
                 <em>DiffCircaPipeline</em>. At each time point (ZT), the
                 average normalized expression for each photoperiod is
                 displayed by a point and error bars indicate the standard error
                 of the mean (SEM). Dashed lines indicate the mean expression
                 across all time points for each photoperiod. Importantly,
                 the cosinor fit performed by <em>DiffCircaPipeline</em> is not
                 always a perfect fit, which can result in error estimating the
                 fit parameters.</p>
                 <p>The right plot displays the distribution of expression for
                 the selected gene as a box-and-whisker plot. The solid line
                 within each box plot indicates the median for the photoperiod.
                 The points represent individual samples, shaped in a violin to
                 represent the distribution of values for each photoperiod.
                 Note: The horizontal position of points is generated
                 stochastically and will change each time the plot is generated,
                 but the y-values are accurate indications of expression per
                 sample.</p>
                 <hr><b>Methodology:</b>
                 <p>Raw count matrices were normalized using <em>DESeq2</em>.
                 For
                 <em>DiffCircaPipeline</em>, known technical batch effects from
                 sample preparation were modeled using <em>Limma</em>.
                 Rhythmicity testing, Cosinor curve fitting, and rhythmic
                 parameter estimation was performed using
                 <em>DiffCircaPipeline</em>. Differential expression analysis
                 was performed using <em>DESeq2</em> and significant genes were
                 defined by a absolute <var>log<sub>2</sub>(Fold Change) > 1</var> and
                 <var>adjusted P-value < 0.05</var>.</p>
                 <hr><b>Citation Information:</b>
                 <ul>
                 <li>DiffCircaPipeline: <small>Xiangning Xue et al,
                 DiffCircaPipeline: a framework for multifaceted
                 characterization of differential rhythmicity\n Bioinformatics,
                 Volume 39, Issue 1, January 2023,
                 btad039\nhttps://doi.org/10.1093/bioinformatics/btad039</small>
                 </li>
                 <li>DESeq2: <small>Love MI, Huber W, Anders S (2014).
                 “Moderated estimation of fold change and dispersion for
                 RNA-seq data with DESeq2.”\nGenome Biology, 15, 550.
                 doi:10.1186/s13059-014-0550-8. </small></li>
                 <li>Limma: <small>Ritchie ME, Phipson B, Wu D, Hu Y, Law CW,
                 Shi W, Smyth GK (2015). “limma powers differential expression
                 analyses for RNA-sequencing and microarray studies.”
                 Nucleic Acids Research, 43(7), e47.\ndoi:10.1093/nar/gkv007.
                 </small></li></ul>"),
      # Action button to generate plot
      actionButton("plotButton", "Update", class = "btn-primary"),
      checkboxInput("show_advanced", "Show advanced options", FALSE),
      conditionalPanel(
        condition = "input.show_advanced == true",

        #advanced inputs
        hr(),
        # Select point size slider
        sliderInput("point.size", "Point size:", min = 0, max = 10, value = 4),
        # Select point size slider
        sliderInput("line.size", "Line size:", min = 0, max = 10, value = 2),
        # Color selector for LONG photoperiod
        colourpicker::colourInput("colLong",
                                  "Select colour",
                                  value = "#e69f00"),
        # Color selector for SHORT photoperiod
        colourpicker::colourInput("colShort",
                                  "Select colour",
                                  value = "#56b4e9")

      ),
      hr(),
      tags$small("This app was developed by Matt Cottam and JP Cartailler at ", a("Creative Data Solutions", href = "https://cds.vanderbilt.edu"))
    ),
    mainPanel(
      conditionalPanel(
        #condition = "input.gene !== 'Aaas'",
        condition = "input.plotButton == 0",
        #advanced inputs
        wellPanel(
          HTML("
            <h3>Welcome circadian and non-circadian scientists!</h3>
            <p>This website profiles Differential Rhythmicity and Differential Expression in the transcriptome of the mouse central circadian clock (suprachiasmatic nuclei; SCN) under summer-like (LONG day) vs. winter-like (SHORT day) seasonal photoperiods. The experimental conditions are LD (h Light: h Dark) 16:8 and LD 8:16, respectively. The direction of comparison is displayed as “the change from SHORT to LONG” (i.e. LONG – SHORT). Further details can be found in our paper, linked below. <em>*Please note that this dataset is restricted to protein-coding genes.</em></p>
            <p>Send questions or comments to <a href='mailto:'>douglas.g.mcmahon@vanderbilt.edu</a>; <a href='mailto:'>olivia.h.cox@vanderbilt.edu</a>. Please cite our work, should you make use of this resource: Cox et al., Journal of Biological Rhythms, 2024, Volume:pages.</p>
            
            <p>To proceed, please select a gene of interest in the drop-down menu and click on Update. The page will load results for gene expression in SHORT vs. LONG photoperiods using the following analyses:</p>
            
            <ul>
            <li>Differential Rhythmicity, via DiffCircaPipeline (Xue et al., 2023 doi: 10.1093/bioinformatics/btad039)</li>
            <li>Differential Expression, via DESeq2 (Love et al., 2014 doi: 10.1186/s13059-014-0550-8)</li>
            </ul>

            <p>Figure 1 illustrates the output of each of these types of analyses.</p>
            
            <img src='Fig1.png' alt='Figure 1' style='width:100%;'>
            
            <p><small><strong>Figure 1. Types of gene expression changes observed in rhythmic genes.</strong> A) Using DESeq2  we tested differential expression, which represents an overall change in the mean expression (timepoints collapsed) between conditions B) Using DiffCircaPipeline (DCP) we examined the following gene expression changes 1)Differential MESOR (Midline Estimating Statistic Of Rhythm) represents a change in the mean of the cosinor fit. 2) Differential peak phase represents a change in the peak expression time of a gene, as shown by a shift along the time axis. 3) Differential rhythmic amplitude represents a change in mean peak-to-trough amplitude.</small></p>
            
            <hr>

            <p><em>I’m not a circadian scientist, why is it important to take circadian rhythms and photoperiod into account when I plan my experiments?</em></p>
            
            <p><strong>Rhythmic gene expression</strong></p>

            <p>Virtually all tissues contain circadian clocks and these “peripheral clocks” allow for expression of genes in a manner specific and relevant to the particular tissue. In fact, between 40-80% of protein-coding genes in mammals are expressed rhythmically in at least one tissue (Mure LS, et. al, 2018 doi: 10.1126/science.aao0318; Zhang R, et. al, 2014 doi: 10.1073/pnas.1408886111). If your gene of interest is rhythmic, it may be important to sample at multiple timepoints. <strong>Fig. 2</strong> shows how single time point experiments can be confounded by shifts in the timing of circadian rhythms across treatments.</p>
            
           <img src='Fig2.png' alt='Figure 2' style='width:100%;'>

            
            <p><small><strong>Figure 2. Sampling across time to accurately measure gene expression in rhythmic genes.</strong> A.) Gene A represents a non-rhythmic gene, in which expression over time remains relatively stable in each individual condition (control = blue, treatment = orange). In this case, sampling at any time of day yields the same result. B.) Gene B represents a rhythmic gene, in which expression over time follows a circadian rhythm in each individual condition (control = blue, treatment = orange). In this case, sampling at different times of day yields different results. Therefore, it is important to collect samples at set intervals throughout the day for each condition. Sampling expression levels across 24-hour timepoints provides information on the temporal patterns of gene expression under different conditions or treatments.</small></p>
            <p><strong>Photoperiod affects rhythmic gene expression.</strong> This website shows that rhythmic gene expression in the mouse brain circadian clock is greatly affected by the light cycle that animals are housed under. Please see our paper to learn more about photoperiod and its influence on the transcriptome.</p>
            

            

            ")
        )
      ),  # condition end

      conditionalPanel(
        #condition = "input.gene !== 'Aaas'",
        condition = "input.plotButton > 0",

        # fluidrow start
        fluidRow(
          wellPanel(
            id = "well_light",
            h3("Characterization of differential rhythmicity (DiffCircaPipeline)"),
            hr(),
            # Show a plot of the gene rythmicity
            column(width = 12,
                   shinycssloaders::withSpinner(plotOutput("rhythmicPlot", height = "80%"))
            ),
            fluidRow(
              column(width = 4),
              column(width = 4,
                     uiOutput("saveUI1")),
              column(width = 4)
            ),
            br(), # spacer
            # Table of parameters
            fluidRow(
              column(
                width = 12,
                tableOutput("DCP_results_table_raw")
              )
            ),
            fluidRow(
              column(
                width = 12,
                tableOutput("DCP_results_table_diff")
              )
            ),

          )

        ), # fluidrow end

        # fluidrow start
        fluidRow(
          wellPanel(
            id = "well_light",
            h3("Differential gene expression (DESeq2)"),
            hr(),
            fluidRow(
              column(6, shinycssloaders::withSpinner(plotOutput("volcanoPlot", height = "100%"))),
              column(6, tableOutput("DE_results_table"))
              
            ),
            fluidRow(
              column(width = 2),
              column(width = 4,
                     uiOutput("saveUI_volcanoPlot")),
              column(width = 6)
            )

          )
        ) # fluidrow end

      ) # condition end

    )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output, session) {

  #shinyjs::logjs(reactive(input$gene))
  shinyjs::logjs(reactive(input$plotButton))


  # Gene selector update

  # Color input 1 (LONG) selector update
  colourpicker::updateColourInput(session, "colLong",
                                  label = "Long (LD 16:8) Color:",
                                  value = "#e69f00",
                                  showColour = "both",
                                  allowTransparent = TRUE)

  # Color input 2 (SHORT) selector update
  colourpicker::updateColourInput(session, "colShort",
                                  label = "Short (LD 8:16) Color:",
                                  value = "#56b4e9",
                                  showColour = "both",
                                  allowTransparent = TRUE)



  # Create a plotting function for the curve fit
  curvefit <- function(gene, point.size, line.size, colLong, colShort) {

    ## Gene info
    gene_ensemblid <- genemap[genemap$mgi_symbol %in% gene,]$ensembl_gene_id_version
    gene_symbol <- genemap[genemap$mgi_symbol %in% gene,]$mgi_symbol

    # Short data
    x1 <- DCP_rhythm$x1
    tod1 <- x1$time
    expr1 <- as.numeric(x1$data[match(gene_ensemblid, x1$gname), ])
    apar1 <- x1$rhythm[match(gene_ensemblid, x1$rhythm$gname), ]
    apar1$Rhythmic <- DCP_rhythm$rhythm.joint$TOJR[match(gene_ensemblid, DCP_rhythm$rhythm.joint$gname)]  
    apar1 <- apar1 %>% mutate(Rhythmic =case_when(Rhythmic =="rhyII"~ 0,
                                                  Rhythmic =="arrhy"~ 0,
                                                  Rhythmic =="rhyI"~ 1,
                                                  Rhythmic =="both"~ 1))
    specInfo1 <- "SHORT"

    # Long data
    x2 <- DCP_rhythm$x2
    tod2 <- x2$time
    expr2 <- as.numeric(x2$data[match(gene_ensemblid, x2$gname), ])
    apar2 <- x2$rhythm[match(gene_ensemblid, x2$rhythm$gname), ]
    apar2$Rhythmic <- DCP_rhythm$rhythm.joint$TOJR[match(gene_ensemblid, DCP_rhythm$rhythm.joint$gname) ]
    apar2<- apar2 %>% mutate(Rhythmic = case_when(Rhythmic =="rhyI"~ 0,
                                                 Rhythmic =="arrhy"~ 0,
                                                 Rhythmic =="rhyII"~ 1,
                                                 Rhythmic =="both"~ 1))
    specInfo2 <- "LONG"

    ## Wrangle the expression data for average line
    res <- cbind(DCP_rhythm$x1$data[gene_ensemblid, ],
                 DCP_rhythm$x2$data[gene_ensemblid, ])
    res <- t(res)
    colnames(res) <- "expr"
    res <- as.data.frame(res)
    res$sample_id <- rownames(res)
    res <- dplyr::left_join(res, metadata, by = c("sample_id"))
    res <- res %>%
      group_by(photoperiod) %>%
      summarize(avgexpr = mean(expr))

    # General settings
    period <- 24


    # Sinusoidal wave functions
    fun.cosinor1 <- function(x) {
      apar1$M + apar1$A * cos(2 * pi / period * x + apar1$phase)
    }

    fun.cosinor2 <- function(x) {
      apar2$M + apar2$A * cos(2 * pi / period * x + apar2$phase)
    }

    # Bind dataframes together
    df1 <- data.frame(Time = tod1, Expression = expr1, Photoperiod = specInfo1)
    df2 <- data.frame(Time = tod2, Expression = expr2, Photoperiod = specInfo2)
    df <- rbind(df1, df2)

    # Summarize to get SEM for points
    df <- df %>%
      group_by(Time, Photoperiod) %>%
      summarize(sd = sd(Expression),
                Expression = mean(Expression),
                sem = sd / sqrt(n()))
    #### New limits to add a buffer to current limits
    myylims <- c(min(df$Expression-df$sem-.1),
                 max(df$Expression+df$sem+.1))
    # Plotting start
    plot <- ggplot2::ggplot(df, aes(x = Time, y = Expression, color = Photoperiod)) +
      ggplot2::geom_point(size = point.size, stat = "summary") +
      geom_errorbar(aes(ymin = Expression - sem, ymax = Expression + sem), width = .8,
                    position = position_dodge(0.05)) +
      scale_color_manual(values = c("SHORT" = colShort, "LONG" = colLong),
                         labels = c("SHORT" = "SHORT (LD 8:16)", "LONG" = "LONG (LD 16:8)"),
                         guide = guide_legend(order = 1)) +
      {if(apar1$Rhythmic == 0){ggplot2::geom_function(fun = fun.cosinor1, size = line.size, color=colShort,linetype='dotdash')}
        else {ggplot2::geom_function(fun = fun.cosinor1, size = line.size, color=colShort)}}+
      {if(apar2$Rhythmic == 0){ggplot2::geom_function(fun = fun.cosinor2,size = line.size, color=colLong,linetype='dotdash')}
        else {ggplot2::geom_function(fun = fun.cosinor2, size = line.size, color=colLong)}}+
      ggplot2::scale_x_continuous(breaks = seq(0, 24, 4), limits = c(-0.4, 24)) +
      ggplot2::scale_y_continuous(n.breaks=6, limits = c(myylims[1], myylims[2]))+
      ggplot2::ggtitle(label = "Cosinor Curve Fit Plot", subtitle = gene_symbol) +
      ggplot2::labs(y = expression(paste(log[2], "(Normalized Expression)", sep = ""))) +
      ggplot2::xlab("Projected Zeitgeber Time") +
      ggpubr::theme_pubr() +
      ggplot2::theme(legend.position = "right",
                     plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
                     plot.subtitle = element_text(size = 18, hjust = 0.5),
                     axis.text = element_text(size = 20, hjust = 0.5),
                     axis.title = element_text(size = 20, hjust = 0.5),
                     legend.text = element_text(size = 14),
                     legend.title = element_text(size = 16),
                     #legend.title = element_blank(),
                     aspect.ratio = 0.6) +
      scale_fill_discrete(name = "", labels = c("SHORT (LD 8:16)", "LONG (LD 16:8)"),
                          guide = guide_legend(order = 1)) +

      # 1. original
      #ggplot2::geom_hline(data = res, aes(yintercept = avgexpr, color = photoperiod), linewidth = 1, linetype = 2) +
      # 2. split
      #ggplot2::geom_hline(data = res, aes(yintercept = avgexpr[1], color = photoperiod[1]), linewidth = 1, linetype = 2) +
      #ggplot2::geom_hline(data = res, aes(yintercept = avgexpr[2], color = photoperiod[2]), linewidth = 1, linetype = 2) +
      # 3. diff't
      geom_hline(aes(yintercept = res$avgexpr[1], linetype = "LONG (LD 16:8)"), size = line.size/2, colour= colLong) +
      geom_hline(aes(yintercept = res$avgexpr[2], linetype = "SHORT (LD 8:16)"), size = line.size/2, colour= colShort) +
      scale_linetype_manual(name = "", values = c('dashed', 'dashed'), guide = guide_legend(order = 2)) +
      guides(
        linetype = guide_legend(
          title="Mean Expression",
          override.aes = list(color = c(colLong, colShort))
          )
      )
    #plot

    return(plot)
  }

  # Create plotting function for the gene expression spread
  exprBoxPlot <- function(gene, point.size, colLong, colShort) {

    # gene="Rtca"
    # point.size = 1
    # colLong = "red"
    # colShort = "blue"

    ## Gene info
    gene_ensemblid = genemap[genemap$mgi_symbol %in% gene,]$ensembl_gene_id_version
    gene_symbol = genemap[genemap$mgi_symbol %in% gene,]$mgi_symbol

    ## Wrangle the expression data
    res <- cbind(DCP_rhythm$x1$data[gene_ensemblid,],
                 DCP_rhythm$x2$data[gene_ensemblid,])
    res <- t(res)
    colnames(res) <- "expr"
    res <- as.data.frame(res)
    res$sample_id <- rownames(res)
    res <- dplyr::left_join(res, metadata, by = c("sample_id"))

    ## Plotting start
    plot <- ggplot2::ggplot(res, aes(x = photoperiod, y = expr, fill = photoperiod)) +
      ggplot2::geom_boxplot(color = "black") +
      ggforce::geom_sina(color = "black", size = point.size - .99, alpha = 0.75) +
      scale_fill_manual(values = c("SHORT" = colShort, "LONG" = colLong), name = "Photoperiod", labels = c("SHORT" = "SHORT (LD 8:16)", "LONG" = "LONG (LD 16:8)")) +
      ggplot2::ggtitle(label = "Differential Box Plot ", subtitle = gene_symbol) +
      ggplot2::labs(y = expression(paste(log[2], "(Normalized Expression)", sep = ""))) +
      ggplot2::xlab("Photoperiod") +
      ggpubr::theme_pubr() +
      ggplot2::theme(legend.position = "bottom",
                     plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
                     plot.subtitle = element_text(size = 18, hjust = 0.5),
                     axis.text = element_text(size = 20, hjust = 0.5),
                     axis.title = element_text(size = 20, hjust = 0.5),
                     legend.text = element_text(size = 14, hjust = 0.5),
                     legend.title = element_text(size = 16, hjust = 0.5),
                     aspect.ratio = 1)

    return(plot)
  }

  # Generate box plot from DESeq2 data
  DEG_plotCounts <- function(gene, point.size, colLong, colShort) {

    ## Gene info
    gene_ensemblid <- genemap[genemap$mgi_symbol %in% gene,]$ensembl_gene_id_version
    gene_symbol <- genemap[genemap$mgi_symbol %in% gene,]$mgi_symbol

    res <- cbind(DEGs_normalized_counts_SHORT[gene_ensemblid,], counts_LONG <- DEGs_normalized_counts_LONG[gene_ensemblid,])
    res <- t(res)
    colnames(res) <- "expr"
    res <- as.data.frame(res)
    res$sample_id <- rownames(res)
    res <- dplyr::left_join(res, metadata, by = c("sample_id"))

    ## Plotting start
    plot <- ggplot2::ggplot(res, aes(x = photoperiod, y = expr, fill = photoperiod)) +
      ggplot2::geom_boxplot(color = "black") +
      ggforce::geom_sina(color = "black", size = point.size - .99, alpha = 0.75) +
      scale_fill_manual(values = c("SHORT" = colShort, "LONG" = colLong), name = "Photoperiod", labels = c("SHORT" = "SHORT (LD 8:16)", "LONG" = "LONG (LD 16:8)")) +
      ggplot2::ggtitle(label = "Counts Box Plot", subtitle = gene_symbol) +
      ggplot2::ylab("Normalized counts") +
      ggplot2::xlab("Photoperiod") +
      ggpubr::theme_pubr() +
      ggplot2::theme(legend.position = "bottom",
                     plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
                     plot.subtitle = element_text(size = 18, hjust = 0.5),
                     axis.text = element_text(size = 20, hjust = 0.5),
                     axis.title = element_text(size = 20, hjust = 0.5),
                     legend.text = element_text(size = 14, hjust = 0.5),
                     legend.title = element_text(size = 16, hjust = 0.5),
                     aspect.ratio = 1)

    return(plot)
  }


  volcanoPlot <- function(gene, point.size, colLong, colShort) {


    library(ggtext)
    # Gene info
    gene_ensemblid = genemap[genemap$mgi_symbol %in% gene,]$ensembl_gene_id_version
    gene_symbol = genemap[genemap$mgi_symbol %in% gene,]$mgi_symbol

    gene_label <- DEGS_all[DEGS_all$gene_symbol %in% gene_symbol,]

    # categorize for coloring
    DEGS_all$de <- "Non-sig."
    DEGS_all$de[DEGS_all$log2FoldChange >= 0.1 & DEGS_all$padj < 0.05] <- "Sig. URG"
    DEGS_all$de[DEGS_all$log2FoldChange < -0.1 & DEGS_all$padj < 0.05] <- "Sig. DRG"
    DEGS_all <- DEGS_all %>% arrange(de)

    DEGS_all$point.size.mag = (point.size * abs(DEGS_all$log2FoldChange)) * 0.35

    plot <- ggplot(DEGS_all, aes(log2FoldChange, -log10(padj))) +
      geom_point(aes(color = de), alpha = 0.5, size=DEGS_all$point.size.mag) +
      scale_color_manual(values = c("Sig. DRG" = color_blue, "Non-sig." = "grey", "Sig. URG" = color_red)) +
      ylab("-log<sub>10</sub>(Adjusted p-value)") +
      xlab("log<sub>2</sub>(Fold change)") +
      labs(color='DE genes') +
      geom_vline(xintercept = 0.1, linetype="dashed", color = "grey") +
      geom_vline(xintercept = -0.1, linetype="dashed", color = "grey") +
      geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "grey") +
    ggplot2::ggtitle(label = "Volcano Plot", subtitle = gene_symbol)

    plot <- plot +
      geom_point(data = gene_label,
                 aes(log2FoldChange, -log10(padj)),
                 col = "black", size = point.size + 1) +
      ggrepel::geom_label_repel(
        size = 8,
        data = gene_label,
        aes(label = gene_symbol),
        segment.color = "black",
        xlim=c(-6,-4),
        ylim=c(20, 30),
        color = "black") +
      ggpubr::theme_pubr() +
      ggplot2::theme(legend.position = "bottom",
                     plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
                     plot.subtitle = element_text(size = 18, hjust = 0.5),
                     axis.text = element_text(size = 20, hjust = 0.5),
                     axis.title = element_text(size = 20, hjust = 0.5),
                     legend.text = element_text(size = 14, hjust = 0.5),
                     legend.title = element_text(size = 16, hjust = 0.5),
                     axis.title.x = element_markdown(),
                     axis.title.y = element_markdown(),
                     aspect.ratio = 1)

    return(plot)
  }


  # Wrangling for the parameters results table
  DCP_params_subset <- DCP_params[, c(3, 1, 23:25, 28:30, 18:20)]
  # DCP_params_subset<- data.frame(lapply( DCP_params_subset, function(x) if(is.numeric(x)) round(x, 5) else x))
  DEGS_all <- DEGs # fo volcano plot
  DEGs_subset <- DEGs[, c(8, 3, 6:7)]
  # DEGs_subset<- data.frame(lapply( DEGs_subset, function(x) if(is.numeric(x)) round(x, 5) else x))
  parametersTable <- dplyr::left_join(DCP_params_subset, DEGs_subset, by = c("ensembl_gene_id" = "ensemblid")) # removed 2023-11-26
  parametersTable_for_DCP_tables <- DCP_params_subset
  parametersTable_for_DCP_tables$ensembl_gene_id <- NULL
  parametersTable <- parametersTable[!parametersTable %>% duplicated(),]

  ## more prep work
  dcp_raw_SHORT_table <-  dplyr::left_join(dcp_raw_SHORT, genemap, by = c("gname" = "ensembl_gene_id_version"))
  dcp_raw_SHORT_table$type <- "SHORT"
  dcp_raw_SHORT_table <- dcp_raw_SHORT_table[, c("type", "mgi_symbol", "M", "A", "peak", "pvalue", "qvalue", "sigma", "R2")]
  dcp_raw_SHORT_table <- data.frame(lapply(dcp_raw_SHORT_table, function(x) if(is.numeric(x)) round(x, 5) else x))

  dcp_raw_LONG_table <-  dplyr::left_join(dcp_raw_LONG, genemap, by = c("gname" = "ensembl_gene_id_version"))
  dcp_raw_LONG_table$type <- "LONG"
  dcp_raw_LONG_table <- dcp_raw_LONG_table[, c("type","mgi_symbol", "M", "A", "peak", "pvalue", "qvalue", "sigma", "R2")]
  dcp_raw_LONG_table <- data.frame(lapply(dcp_raw_LONG_table, function(x) if(is.numeric(x)) round(x, 5) else x))
  dcp_raw_COMBINED_table <- rbind(dcp_raw_SHORT_table, dcp_raw_LONG_table)


  rv <- reactiveValues(plot1 = NULL,
                       plot2 = NULL,
                       volcanoPlot = NULL,
                       DCP_results_table_raw = NULL,
                       dcp_raw_SHORT_table = NULL,
                       dcp_raw_LONG_table = NULL,
                       dcp_raw_COMBINED_table = NULL,
                       DE_results_table = NULL,
                       DCP_results_table_diff = NULL,
                       save = NULL) # Prime reactive values




 add_asterisk <- function(x) {
   # Check if x is numeric and determine formatting based on its range
   formatted_x <- if (is.numeric(x)) {
     scientific_range <- x > 0 & x < 1e-2
     ifelse(scientific_range,
            format(signif(x, 4), scientific = TRUE),  # Scientific notation for specified range
            sprintf("%.4f", round(x, digits = 4)))    # Round to 4 decimal places otherwise
   } else {
     format(x, nsmall = 4)  # Default formatting for non-numeric
   }

   # Append asterisks for significance levels after formatting
   ifelse(x <= 0.05, paste0(x, "*"), x)
 }



  # Observe plot button click
  observeEvent(input$plotButton, {

    # DCP cosinor curve fit ----
    rv$plot1 <- curvefit(input$gene, input$point.size, input$line.size, input$colLong, input$colShort)

    # DCP box plot ----
    rv$plot2 <- exprBoxPlot(input$gene, input$point.size, input$colLong, input$colShort)

    # DCP raw data table ----
    #
    # TOJR column in p-value sig for SHORT photoperiod if it says “rhyI”, for LONG
    # photoperiod if it says “rhyII”, in BOTH conditions if it says “Both” and not
    # significant in any condition if it says “arrhy”.
    # For example, if a gene says “rhyI”, then in the table for that gene, the
    # SHORT photoperiod will be “yes” and the information will be populated in
    # the rest of the columns. However, for this gene, because it only says rhyI,
    # the LONG photoperiod will be “no” and the rest of the columns for long
    # photoperiod will say “N/A”.

    dcp_raw_COMBINED_table <- dcp_raw_COMBINED_table %>% mutate(across('type', str_replace, 'SHORT', 'SHORT (8:16)'))
    dcp_raw_COMBINED_table <- dcp_raw_COMBINED_table %>% mutate(across('type', str_replace, 'LONG', 'LONG (16:8)'))

    # now join in DCP_rhythm$rhythm_joint
    d <- merge(dcp_raw_COMBINED_table, genemap, by="mgi_symbol")
    d2 <- merge(d, DCP_rhythm$rhythm.joint, by.x="ensembl_gene_id_version", by.y="gname")
    #head(d2)

    d2 = d2 %>% mutate(
      rhythmic = case_when(
        type=="SHORT (8:16)" & TOJR=="rhyI" ~ "Yes",
        type=="SHORT (8:16)" & TOJR=="both" ~ "Yes",
        type=="LONG (16:8)" & TOJR=="rhyII" ~ "Yes",
        type=="LONG (16:8)" & TOJR=="both" ~ "Yes",
        TRUE ~ "No"
      )
    )

    # if not rhytmic, set columns to NA (don't do string, we need to round)
    d2[d2$rhythmic == "No", c( "M", "A", "peak")] <- NA

    # if not rhytmic, set columns to NA (don't do string, we need to round)
    d2[d2$rhythmic == "Yes", c("rhythmic")] <- paste0('<span style="font-weight:bold; color:', color_green, '">Yes</span')
    d2[d2$rhythmic == "No", c("rhythmic")] <- paste0('<span style="color:', color_red, '">No</span')

    dcp_raw_COMBINED_table <- d2[, c("mgi_symbol", "type", "rhythmic", "R2", "M", "A", "peak")]

    # First check if there is any available data
    if(nrow(dcp_raw_COMBINED_table %>% filter(mgi_symbol == input$gene)) > 0){


      rv$DCP_results_table_raw <- dcp_raw_COMBINED_table %>%
        filter(mgi_symbol == input$gene) %>%
        #mutate(across(4:7, round, 3)) %>%
        knitr::kable(
          format = "html",
          digits = 20,
          col.names = c("Gene Symbol", "Photoperiod", "Rhythmic?", "Cosinor rhythm fitness (R<sup>2</sup>)", "Mesor (log<sub>2</sub>(TPM))", "Amplitude (log<sub>2</sub>(TPM)", "Peak phase (h)"),
          align = "c",
          booktabs = T,
          escape = F) %>%
        kable_styling(bootstrap_options = "bordered", full_width = T) %>%
        row_spec(0, background = "#f0f0f0") %>%
        #column_spec(4, background = "#eeeeee") %>%
        #column_spec(5, background = "#eeeeee") %>%
        #column_spec(6, background = "#eeeeee") %>%
        #column_spec(7, background = "#eeeeee") %>%
        add_header_above(c("Rhythmic Parameters" = 7),
                         bold = TRUE)

    }else{

      rv$DCP_results_table_raw <- "No data available for this gene"

    }

    # DCP raw data ----
    # First check if there is any available data
    if(nrow(parametersTable_for_DCP_tables %>%  filter(mgi_symbol == input$gene)) > 0){

      # reorder table
      parametersTable_for_DCP_tables_reordered <- parametersTable_for_DCP_tables[, c("mgi_symbol", "delta.M", "p.delta.M", "q.delta.M", "delta.A", "p.delta.A", "q.delta.A", "delta.peak", "p.delta.peak", "q.delta.peak")]

      rv$DCP_results_table_diff <- parametersTable_for_DCP_tables_reordered %>%
        filter(mgi_symbol == input$gene) %>%
         mutate(across(4, add_asterisk)) %>%
         mutate(across(6, add_asterisk)) %>%
         mutate(across(7, add_asterisk)) %>%
         mutate(across(9, add_asterisk)) %>%
         mutate(across(10, add_asterisk)) %>%
        knitr::kable(
          format = "html",
          digits = 100,
          col.names = c("Gene Symbol", "Δ Mesor (log<sub>2</sub>(TPM))", "p-value", "q-value", "Δ Amplitude (log<sub>2</sub>(TPM))", "p-value", "q-value", "Δ Peak phase (h)", "p-value", "q-value"),
          align = "c",
          booktabs = T,
          escape = FALSE) %>%
        # Cell styling
        kable_styling(bootstrap_options = "bordered", full_width = T) %>%
        row_spec(0, background = "#f0f0f0") %>%
        #column_spec(6, background = "red") %>%
        add_header_above(c("",
                           "Mesor" = 3,
                           "Amplitude" = 3,
                           "Peak phase" = 3), background = "#f0f0f0") %>%
        add_header_above(c("Differential Rhythmic Parameters (LONG vs SHORT)" = 10))


    }else{

      rv$DCP_results_table_diff <- "Either no data is available for this gene or only one photoperiod was rhythmic."

    }

    # DE volcano plot ----
    rv$volcanoPlot <- volcanoPlot(input$gene, input$point.size, input$colLong, input$colShort)

    # DE data table ----
    #
    # Data wrangling (we could do this offline...)

    DEGs_subset <- parametersTable[, c(1, 2, 12:14)]
    DEGs_subset <- DEGS_all[, c("gene_symbol", "ensemblid", "log2FoldChange", "pvalue", "padj")]
    colnames(DEGs_subset) <- c("mgi_symbol", "ensembl_gene_id", "log2FoldChange", "pvalue","padj")

    rownames(DEGs_normalized_counts_SHORT) <- sub("\\.[0-9]*$", "", rownames(DEGs_normalized_counts_SHORT))
    rownames(DEGs_normalized_counts_LONG) <- sub("\\.[0-9]*$", "", rownames(DEGs_normalized_counts_LONG))

    mean_short <- as.data.frame(rowMeans(DEGs_normalized_counts_SHORT))
    colnames(mean_short) <- c("mean_short")
    mean_short$ensemblid <- rownames(mean_short)

    mean_long <- as.data.frame(rowMeans(DEGs_normalized_counts_LONG))
    colnames(mean_long) <- c("mean_long")
    mean_long$ensemblid <- rownames(mean_long)

    DEGs_subset_expanded <- dplyr::left_join(DEGs_subset, mean_long, by = c("ensembl_gene_id" = "ensemblid"))
    DEGs_subset_expanded <- dplyr::left_join(DEGs_subset_expanded, mean_short, by = c("ensembl_gene_id" = "ensemblid"))
    
    DEGs_subset_expanded$ensembl_gene_id <- NULL

    DEGs_subset_expanded$significant[DEGs_subset_expanded$padj < 0.05] = "*"

    format_scientific <- function(x) {

      if (is.numeric(x)) {
        scientific_range <- abs(x) > 0 & abs(x) < 1e-04
        result <- ifelse(scientific_range, format(signif(x, 4), scientific = TRUE), format(round(x, digits = 4), nsmall = 4))
        return(result)
      } else {
        result <- format(x, nsmall = 4)
        return(result)
      }
    }

    DEGs_subset_expanded <- DEGs_subset_expanded %>% mutate(pvalue = format_scientific(pvalue))
    DEGs_subset_expanded <- DEGs_subset_expanded %>% mutate(padj = format_scientific(padj))

    DEGs_subset_expanded$padj <- ifelse(is.na(DEGs_subset_expanded$significant), DEGs_subset_expanded$padj, paste(DEGs_subset_expanded$padj, DEGs_subset_expanded$significant, sep = ""))
    DEGs_subset_expanded$significant <- NULL

    DEGs_subset_expanded_no_counts <- DEGs_subset_expanded %>%  select(-c("mean_long", "mean_short"))
    
    rv$DE_results_table <- DEGs_subset_expanded_no_counts %>%
      filter(mgi_symbol == input$gene) %>%
      `colnames<-`(c("Gene Symbol", "log<sub>2</sub>(Fold change)", "p-value", "Adj. p-value")) %>%
      knitr::kable(
        format = "html",
        align = "c",
        booktabs = T,
        escape = FALSE) %>%
      kable_styling(bootstrap_options = "bordered", full_width = T) %>%
      row_spec(0, background = "#f0f0f0") %>%
      add_header_above(c("Differential Expression of LONG vs SHORT photoperiod" = 4))


    rv$save <- TRUE

  })

  # Generate the curve fit plot
  output$rhythmicPlot <- renderPlot(
    rv$plot1,
    height = 400
  )

  # Generate the expression box plot
  output$boxPlot <- renderPlot(
    rv$plot2,
    height = 600
  )

  # Generate the results table
  output$DCP_results_table_raw <- renderText(rv$DCP_results_table_raw)
  output$DCP_results_table_diff <- renderText(rv$DCP_results_table_diff)


  # Plot saving functions
  output$savePlot1 <- downloadHandler(
    file = paste0("DCP_CurveFit__", input$gene, ".png"),
    content = function(file) {
      psave1 <- rv$plot1 +
        labs(caption = paste0("If using this image, please cite our work:\n", citation))
      ggsave(plot = psave1,
             filename = file,
             dpi = 300,
             width = 5,
             height = 6,
             units = "in",
             scale = 2)
    }
  )

  output$saveUI1 <- renderUI({
    if (!is.null(rv$save)) {
      downloadButton("savePlot1", "Download curve fit plot")
    }
  })


  output$savePlot2 <- downloadHandler(
    file = paste0("DCP_MeanExpression__", input$gene, ".png"),
    content = function(file) {
      psave2 <- rv$plot2 +
        labs(caption = paste0("If using this image, please cite our work:\n", citation))
      ggsave(plot = psave2,
             filename = file,
             dpi = 300,
             width = 5,
             height = 6,
             units = "in",
             scale = 2)
    })


  ## DESeq2 section
 
  # Generate the expression box plot
  output$volcanoPlot <- renderPlot(
    rv$volcanoPlot,
    height = 400
  )

  output$savePlot_volcanoPlot <- downloadHandler(
    file = paste0("DEG_volcanoplot__", input$gene, ".png"),
    content = function(file) {
      psave2 <- rv$volcanoPlot +
        labs(caption = paste0("If using this image, please cite our work:\n", citation))
      ggsave(plot = psave2,
             filename = file,
             dpi = 300,
             width = 5,
             height = 6,
             units = "in",
             scale = 1.5)
    })

  output$saveUI_volcanoPlot <- renderUI({
    if (!is.null(rv$save)) {
      downloadButton("savePlot_volcanoPlot", "Download DE volcano plot")
    }
  })

  output$DE_results_table <- renderText(rv$DE_results_table)

}

# Run the application
if(password_protect == TRUE){
  shinyApp(polished::secure_ui(ui), polished::secure_server(server))
}else{
  shinyApp(ui = ui, server = server)
}


