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

#' VJ combinations probabilities for mouse TRB model
#'
#' A data frame with probabilities of VJ combinations
#' All probabilities based on OLGA mouse T beta model
#'
#' @source \url{https://github.com/statbiophys/OLGA/tree/master/olga/default_models/mouse_T_beta}
#'
#' @author Pavel V. Shelyakin
"OLGAVJ_MOUSE_TRB"

#' VJ combinations probabilities for mouse TRB model obtained from C57BL/6
#' unproductive sequences
#'
#' A data frame with probabilities of VJ combinations
#' All probabilities based on unproductive sequences from C57BL/6 mice line
#'
#' @source \url{https://github.com/KseniaMIPT/tcrgrapher/tree/master/models/C57BL.6_mouseTRB}
#'
#' @author Ksenia R. Lupyr
"C57BL6_MOUSE_TRB"

#' VJ combinations probabilities for human TRB model
#'
#' A data frame with probabilities of VJ combinations
#' All probabilities based on OLGA human T beta model
#'
#' @source \url{https://github.com/pogorely/ALICE/blob/master/OLGA_V_J_hum_beta.rda}
"OLGAVJ_HUMAN_TRB"

#' VJ combinations probabilities for human TRA model
#'
#' A data frame with probabilities of VJ combinations
#' All probabilities based on OLGA human T beta model
#'
#' @source \url{https://github.com/pogorely/ALICE/blob/master/OLGA_A_hum_alpha.rda}
"OLGAVJ_HUMAN_TRA"
