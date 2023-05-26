#' @import methods

setClass('TCRgrapher',
         slots = c(
           clonoset  = 'data.table',
           metadata = 'data.table'
         ))

setClass('TCRgrapherCounts',
         contains = 'TCRgrapher',
         slots = c(
           count_table = 'data.frame',
           feature_info = 'data.table'
         ))

# setClass('TCRgrALICE',
#          contains = 'TCRgrapher',
#          slots = c(
#            Q_val = 'numeric',
#            thres_counts = 'numeric',
#            N_neighbors_thres = 'numeric',
#            p_adjust_method = "character",
#            chain = 'character',
#            stats = 'character',
#            model = 'character',
#            OLGAVJ = 'list'
#          ))
