#----------------------------------------------------------------------------#
####   Function to check if transition probability array/matrix  is valid ####
#----------------------------------------------------------------------------#
#' Check if transition array is valid
#'
#' \code{check_transition_probability} checks if transition probabilities are in \[0, 1\].
#'
#' @param a_P A transition probability array.
#' @param err_stop Logical variable to stop model run if set up as TRUE. Default = FALSE.
#' @param verbose Logical variable to indicate print out of messages. 
#' Default = FALSE
#'
#' @return
#' This function stops if transition probability array is not valid and shows 
#' what are the entries that are not valid
#' @import utils
#' @export
check_transition_probability <- function(a_P,
                                         err_stop = FALSE, 
                                         verbose = FALSE) {
  a_P <- as.array(a_P)
  m_indices_notvalid <- arrayInd(which(a_P < 0 | a_P > 1), 
                                 dim(a_P))
  if(dim(m_indices_notvalid)[1] != 0){
    v_rows_notval   <- rownames(a_P)[m_indices_notvalid[, 1]]
    v_cols_notval   <- colnames(a_P)[m_indices_notvalid[, 2]]
    v_cycles_notval <- dimnames(a_P)[[3]][m_indices_notvalid[, 3]]
    df_notvalid <- data.frame(`Transition probabilities not valid:` = 
                                matrix(paste0(paste(v_rows_notval, v_cols_notval, sep = "->"),
                                              "; at cycle ",
                                              v_cycles_notval), ncol = 1), 
                              check.names = FALSE)
    if(err_stop) {
      stop("Not valid transition probabilities\n",
           paste(capture.output(df_notvalid), collapse = "\n"))
    }
    if(verbose){
      warning("Not valid transition probabilities\n",
              paste(capture.output(df_notvalid), collapse = "\n"))
    } 
  }
}
#----------------------------------------------------------------------------#
####   Function to check if sum of transition probabilities equal to one  ####
#----------------------------------------------------------------------------#
#' Check if the sum of transition probabilities equal to one. 
#'
#' \code{check_sum_of_transition_array} checks if each of the rows of the 
#' transition matrices sum to one. 
#' 
#' @param a_P A transition probability array.
#' @param n_states Number of health states.
#' @param n_t Number of cycles.
#' @param err_stop Logical variable to stop model run if set up as TRUE. Default = FALSE.
#' @param verbose Logical variable to indicate print out of messages. 
#' Default = FALSE
#' @return 
#' The transition probability array and the cohort trace matrix.
#' @import dplyr
#' @export
check_sum_of_transition_array <- function(a_P,
                                          n_states,
                                          n_t,  
                                          err_stop = FALSE, 
                                          verbose  = FALSE) {
  a_P <- as.array(a_P)
  valid <- (apply(a_P, 3, function(x) sum(rowSums(x))) == n_states)
  if (!isTRUE(all.equal(as.numeric(sum(valid)), as.numeric(n_t)))) {
    if(err_stop) {
      stop("This is not a valid transition matrix")
    }
    if(verbose){
      warning("This is not a valid transition matrix")
    } 
  }
}











