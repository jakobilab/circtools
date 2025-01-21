#!/usr/bin/env Rscript

# suppress loading messages
suppressMessages(library(formattable))
suppressMessages(library(kableExtra))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
# suppressMessages(library(colortools))

# switch to red warning color if more blast hits are found
high_count_number = 0

args <- commandArgs(trailingOnly = TRUE)

experiment_name <- args[2]

no_svg_flag <- args[3]
#svg_dir <- args[4]

# set output to HTML
options(knitr.table.format = 'html')

# generic HTML header with bootstrap JS libraries
#############################################################################################################
html_header="
<html>
<head>

<meta charset=\"utf-8\">
<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">

  <link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css\">
  <script src=\"https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js\"></script>
  <script src=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js\"></script>

<script>
$(document).ready(function(){
    $('[data-toggle=\"popover\"]').popover();
});
</script>

<script>
$(document).ready(function(){
    $('[data-toggle=\"tooltip\"]').tooltip();
});
</script>

   <style type=\"text/css\">
        /* The max width is dependant on the container (more info below) */
        .popover{
            max-width: 50%; /* Max Width of the popover (depending on the container!) */
        }
    </style>

</head>
<body>"

html_header <- paste(html_header,"<h1>Circtools padlock probe design results for experiment ", experiment_name , "</h1>" , sep="")

html_header <- paste(html_header, "
<div class='my-legend'>
  <div class='legend-title'>Colour legends for TM and GC columns </div>
  <div class='legend-scale'>
    <ul class='legend-labels'>
      <li><span style='background:#4575B4;'></span>Low</li>
      <li><span style='background:#ABD9E9;'></span></li>
      <li><span style='background:#FEE090;'></span>Medium</li>
      <li><span style='background:#F46D43;'></span></li>
      <li><span style='background:#D73027;'></span>High</li>
    </ul>
  </div>
  </div>
  
    <style type='text/css'>
    .my-legend .legend-title {
      text-align: right;
      margin-bottom: 8px;
      font-weight: bold;
      font-size: 90%;
      }
    .my-legend .legend-scale ul {
      margin: 0;
      padding: 0;
      float: right;
      list-style: none;
      }
    .my-legend .legend-scale ul li {
      display: block;
      float: left;
      width: 50px;
      margin-bottom: 6px;
      text-align: center;
      font-size: 80%;
      list-style: none;
      }
    .my-legend ul.legend-labels li span {
      display: block;
      float: right;
      height: 15px;
      width: 50px;
      }
    .my-legend .legend-source {
      font-size: 70%;
      color: #999;
      clear: both;
      }
    .my-legend a {
      color: #777;
      }
  </style>", sep="\n")

#############################################################################################################

# generate a divergent color scale with 11 shades
color_palette <- rev(brewer.pal(n = 8, name = 'RdYlBu'))

#default TM value
default_tm_value <- 50
default_gc_value <- 40
default_product_value <- 1
#default_product_value <- "preferred"

construct_color_column <- function(column, default_value, palette, top_val, bottom_val)
{
    ## give different min and max values for Tm and GC content per column
    #top_val <- (max(column, na.rm=T) - default_value)
    #bottom_val <- (default_value - min(column, na.rm=T))

    if (top_val > bottom_val){
        from <- default_value - top_val
        to <- default_value + top_val
    } else if (top_val == bottom_val) {
        from <- default_value - 1
        to <- default_value + 1
    } else {
        from <- default_value - bottom_val
        to <- default_value + bottom_val
    }

    return(as.character(cut(column,seq( from, to, length.out= length(palette)+1 ), labels=palette, include.lowest = T)))
}

# read data file name from args
data_file_name <- args[1]

# read whole file into data table
data_table <- read.csv(data_file_name, header = FALSE, sep = "\t")
data_table$circid <- paste(data_table$V1,data_table$V2,data_table$V3,data_table$V4,data_table$V5,data_table$V6,sep="_")
#print(paste(data_table$V1,data_table$V2,data_table$V3, data_table$V4, data_table$V5,data_table$V6))
data_table$circid <- paste(sep="","<img src=",data_table$circid,".svg>")
#print(head(data_table$V6))

## remove unused columns
# commented this because probe design script does not have column 6 -> junction flag from CircCoordinate and 9 and 10 because probe junctions do not have primer start/end regions
#data_table <- data_table[-c(6,9,10)]

# correctly name the columns
colnames(data_table) <- c(  "Annotation",
                            "Chr",
                            "Start",
                            "Stop",
                            "Strand",
                            "Left_",
                            "Right_",
                            "TM_left",
                            "TM_right",
                            "TM_Full",
                            "GC_left",
                            "GC_right",
                            "Product_size",
                            "BLAST_left",
                            "BLAST_right",
                            "ID"
                            )

## for padlock probe designing, product_size is not really needed column. So replaced this with ligation junction preference column. 
## Ligation junction 1 if preferred or 0 if neutral
#data_table$Product_size <- ifelse(data_table$Product_size == "preferred", 1, 0)

data_table$right_tm_color  = construct_color_column(data_table$TM_right,default_tm_value,color_palette, 70, 50)
data_table$left_tm_color   = construct_color_column(data_table$TM_left,default_tm_value,color_palette, 70, 50)
data_table$full_tm_color   = construct_color_column(data_table$TM_Full,default_tm_value,color_palette, 82, 68)

data_table$left_gc_color   = construct_color_column(data_table$GC_left,default_gc_value,color_palette, 60, 35)
data_table$right_gc_color  = construct_color_column(data_table$GC_right,default_gc_value,color_palette, 60, 35)

#data_table$product_color  = construct_color_column(data_table$Product_size,default_product_value,color_palette)
data_table$product_color <- ifelse(data_table$Product_size == "preferred", "#99d594", "#ffffbf")  # manually adding colors for ligation junction as this is a character value

colnames_final <- c(        "Annotation",
                            "Chr",
                            "Start",
                            "Stop",
                            "Strand",
                            "TM RBD5",
                            "TM RBD3",
                            "TM Full",
                            "GC% RBD5",
                            "GC% RBD3",
                            "Ligation Junction",
                            "RBD5",
                            "BLAST",
                            "RBD3",
                            "BLAST"
                    )

# run_primer_design a column with BLAST hit counts
data_table$BLAST_left_count <- lengths(regmatches(data_table$BLAST_left, gregexpr(";", data_table$BLAST_left))) + 1
data_table$BLAST_right_count <- lengths(regmatches(data_table$BLAST_right, gregexpr(";", data_table$BLAST_right))) + 1

data_table$BLAST_left_count[data_table$BLAST_left_count == 1] = 0
data_table$BLAST_right_count[data_table$BLAST_right_count == 1] = 0

# replace ; with HTML linebreaks for hover popover text
data_table$BLAST_left <- gsub(";", "<br/><br/>", data_table$BLAST_left)
data_table$BLAST_right <- gsub(";", "<br/><br/>", data_table$BLAST_right)

# remove 0 entries from location columns for provided circRNA FASTA files
data_table$Chr <- gsub("\\b0\\b", "", data_table$Chr )
data_table$Start <- gsub("\\b0\\b", "", data_table$Start )

data_table$Stop  <- ifelse( data_table$Start == "" , "", data_table$Stop )

data_table$Strand<- gsub("\\b0\\b", "", data_table$Strand )


# clean up rownames to hide them lateron
rownames(data_table) <- c()

# main output table generation
if (no_svg_flag == "FALSE"){
  output_table <- data_table %>%
    mutate(
      #Product_size = color_bar(product_color)(Product_size),
      
      Forward = cell_spec(escape = F, Left_, popover = spec_popover( title = "Graphical represensation of designed probes\"data-html=\"True\"", position = "left", content =ID ), 
                          background = ifelse(BLAST_left_count > high_count_number, "red", "darkgreen"),
                          color = ifelse(BLAST_left_count > high_count_number, "white", "white")),
      
      L = cell_spec(paste(BLAST_left_count),
                    popover = spec_popover(content = BLAST_left, title = "Blast Hits\"data-html=\"True\"", position = "right"),
                    background = ifelse(BLAST_left_count > high_count_number, "red", "darkgreen"),
                    color = ifelse(BLAST_left_count > high_count_number, "white", "white"), bold = "true"),
      
      Reverse = cell_spec(escape = F, Right_, popover = spec_popover( title = "Graphical represensation of designed probes\"data-html=\"True\"", position = "left", content =ID ), 
                background = ifelse(BLAST_right_count > high_count_number, "red", "darkgreen"),
                color = ifelse(BLAST_right_count > high_count_number, "white", "white")),
      
      R = cell_spec(paste(BLAST_right_count),
                    popover = spec_popover(content = BLAST_right, title = "Blast Hits\"data-html=\"True\"", position = "left"),
                    background = ifelse(BLAST_right_count > high_count_number, "red", "darkgreen"),
                    color = ifelse(BLAST_right_count > high_count_number, "white", "white"), bold = "true"),
      TM_left = color_bar(left_tm_color)(TM_left),
      TM_right = color_bar(right_tm_color)(TM_right),
      TM_Full = color_bar(full_tm_color)(TM_Full),
      
      
      GC_left = color_bar(left_gc_color)(GC_left),
      GC_right = color_bar(right_gc_color)(GC_right),
      
      
      Strand = formatter('span', style = style(font.weight = "bold"))(Strand)
    )
} else {
  output_table <- data_table %>%
    mutate(
      Forward = cell_spec(escape = F, Left_, background = ifelse(BLAST_left_count > high_count_number, "red", "darkgreen"), 
                                                           color = ifelse(BLAST_left_count > high_count_number, "white", "white")),
      
      L = cell_spec(paste(BLAST_left_count),
                    popover = spec_popover(content = BLAST_left, title = "Blast Hits\"data-html=\"True\"", position = "right"),
                    background = ifelse(BLAST_left_count > high_count_number, "red", "darkgreen"),
                    color = ifelse(BLAST_left_count > high_count_number, "white", "white"), bold = "true"),
      
      Reverse = cell_spec(escape = F, Right_, background = ifelse(BLAST_right_count > high_count_number, "red", "darkgreen"), 
                                                            color = ifelse(BLAST_left_count > high_count_number, "white", "white")),
      
      R = cell_spec(paste(BLAST_right_count),
                    popover = spec_popover(content = BLAST_right, title = "Blast Hits\"data-html=\"True\"", position = "left"),
                    background = ifelse(BLAST_right_count > high_count_number, "red", "darkgreen"),
                    color = ifelse(BLAST_right_count > high_count_number, "white", "white"), bold = "true"),
      TM_left = color_bar(left_tm_color)(TM_left),
      TM_right = color_bar(right_tm_color)(TM_right),
      TM_Full = color_bar(full_tm_color)(TM_Full),
      
      
      GC_left = color_bar(left_gc_color)(GC_left),
      GC_right = color_bar(right_gc_color)(GC_right),
      
      
      Strand = formatter('span', style = style(font.weight = "bold"))(Strand)
      )
}
 
 output_table <- output_table %>%
    select(- Left_) %>%
    select(- Right_) %>%
    select(- BLAST_left) %>%
    select(- BLAST_right) %>%
    select(- BLAST_left_count) %>%
    select(- BLAST_right_count) %>%
    select(- ID) %>%
    select(- right_tm_color) %>%
    select(- left_tm_color) %>%
    select(- full_tm_color) %>%
    select(- right_gc_color) %>%
    select(- left_gc_color) %>%
    select(- product_color) %>%
    select(Annotation, everything()) %>%
    kable("html", escape = F, col.names=colnames_final) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "responsive"), full_width = T) %>%
    # column_spec(5, width = "3cm")
    add_header_above(c("Input IDs" = 5, "Designed Probes" = 10))
    # group_rows("Group 1", 4, 7) %>%
    # group_rows("Group 1", 8, 10)
    # collapse_rows(columns = 1)
write(paste(html_header, output_table, sep=""), file = "")
