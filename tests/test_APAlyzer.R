test_judge_rep <- function(){
    checkEquals(.judge_rep("CON1", "TRT1"), 'single')
}

test_absMIN <- function(){
    checkEquals(.absMIN(c(-1, 2, 3, 4, -10)), 1)
}
