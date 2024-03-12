#     TCRgrapher: R package for identifying condition associated T cell clonotypes
#     Copyright (C) 2024 Kseniia Lupyr
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @import methods

setClass('TCRgrapher',
         slots = c(
           clonoset  = 'data.table',
           metadata = 'data.table',
           edges = 'matrix'
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
