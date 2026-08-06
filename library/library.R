
############################################################
#Functions for multiple subsampling (rarefaction) of species counts tables
# and analysis of rarefied diversity metrics

#############################################################################
#perform rarefaction
multiple_subsamples = function(x = NULL, depth = NULL, iterations = NULL){
    # Performs multiple subsamples with vegan::rrarefy and returns 
    # a list of length iterations containing the resulting rarefied tables.
    # Attempts to convert non matrix objects (e.g., data.frame) to 
    # matrix and filters out columns with less than depth.
    # Input: x is a table of samples (rows) x species (columns)
    # depth is the desired sample depth per sample
    # iterations is the number of samples
    
    if(require(vegan) != T){
        install.packages(vegan)
    }
    library(vegan)
    
    if(is.matrix(x) == F){
        x = as.matrix(x)
    }
    
    x.min = x[rowSums(x) >= depth,]
    
    x_subsamples = list()
    
    for(i in 1:iterations){
        #rrarefy apparently throws a warnings if there are no counts of 1 in the data
        # annoying... but if you get that warning ignore
        x_subsamples[[i]] = rrarefy(
            x = x.min,
            sample = depth
        )
    }
    return(x_subsamples)
}

#perform rarefaction, no auto discard of samples below min
multiple_subsamples_no_min = function(x = NULL, depth = NULL, iterations = NULL){
    # Performs multiple subsamples with vegan::rrarefy and returns
    # a list of length iterations containing the resulting rarefied tables.
    # Attempts to convert non matrix objects (e.g., data.frame) to matrix
    # Input: x is a table of samples (rows) x species (columns)
    # depth is the desired sample depth per sample
    # iterations is the number of samples
    
    if(require(vegan) != T){
        install.packages(vegan)
    }
    library(vegan)
    
    if(is.matrix(x) == F){
        x = as.matrix(x)
    }
    
    x_subsamples = list()
    
    for(i in 1:iterations){
        #rrarefy apparently throws a warnings if there are no counts of 1 in the data
        # annoying... but if you get that warning ignore
        x_subsamples[[i]] = rrarefy(
            x = x,
            sample = depth
        )
    }
    return(x_subsamples)
}

#############################################################################
#calculate richness of each sample (i.e., number of non zero entries per row)
richness_calc = function(x){
    x[x>0] = 1
    return(rowSums(x))
}
#############################################################################
#function to log+1 transform counts before dist
log_dist = function(x, method = "bray"){
    x_log = log(x+1)
    vegdist(x_log, method = method, binary = F, diag = T, upper = T)
}

#############################################################################
#Take avgs over calculated metrics and/or the subsampled count tables
#i.e., average a list of matrices
avg_matrix_list = function(x){
    #function to average a list of matrices element-wise
    # silently converts list elements to matrix (from e.g., dist or df)
    # returns matrix with original row and col names if any
    list_len = length(x)
    
    #convert to matrix
    x = lapply(x, as.matrix)
    
    #get row and col names.
    names_row = rownames(x[[list_len]])
    names_col = colnames(x[[list_len]])
    
    sum_mat = Reduce('+', x)
    avg_mat = sum_mat/list_len
    
    #add names
    rownames(avg_mat) = names_row
    colnames(avg_mat) = names_col
    
    return(avg_mat)
}
#############################################################################


#############################################################################
############################
#ggplot2 themes and palettes

library(ggplot2)
my_gg_theme = theme(
    panel.background = element_rect(fill='white', colour='black'),
    panel.grid.major=element_blank(),
    panel.grid.minor= element_blank(),
    text=element_text(family="sans"),
    axis.text=element_text(size=15, color="black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust=0, size=20),
    axis.title = element_text(size=17),
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    strip.text = element_text(size = 15),
    axis.title.x = element_text(margin = margin(t= 10)),
    axis.title.y = element_text(margin = margin(r=10)),
    legend.key = element_blank()
)

my_gg_theme.def_size = theme(
    panel.background = element_rect(fill='white', colour='black'),
    panel.grid.major=element_blank(),
    panel.grid.minor= element_blank(),
    text=element_text(family="sans"),
    axis.text=element_text(color="black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust=0),
    #axis.title = element_text(size=17),
    #legend.title = element_blank(),
    #legend.text = element_text(size = 19),
    #strip.text = element_text(size = 15),
    axis.title.x = element_text(margin = margin(t= 10)),
    axis.title.y = element_text(margin = margin(r=10)),
legend.key = element_blank()
)


#function to convert integer dates back to breaks and date labels
#use with scale_x_gradient2 if it is desired to define a midpoint
date_breaks <- function(x){
    #breaks <- c(min(x),median(x),max(x))
    breaks <- quantile(x, probs = seq(0, 1, 0.25))
    attr(breaks,"labels") <- month.name[as.numeric(format(as.Date(breaks, origin="1970-01-01"),"%m"))]
    names(breaks) <- attr(breaks,"labels")
    return(breaks)
}
#April    May August
#19460  19502  19586
#mean of min-max is 19523

#April   April     May    June  August
#19460.0 19477.5 19502.0 19516.0 19586.0

#palette for date breaks
five_cols_gradient_palette = c('#b2182b','#d6604d','#f4a582','white','#2166ac')


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

twelvePaired <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

############################
#############################################################################
