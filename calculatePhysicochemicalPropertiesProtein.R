#' @name calculateGravyScore
#' 
#' @title Calculate Gravy scores of an amino acid/protein
#' 
#' @description 
#' The function `calculateGravyScore` calculates the GRAVY (grand average of 
#' hydropathy) value of a given amino acid/protein sequence. 
#' The GRAVY value is calculated by adding the hydropathy value for each 
#' residue and dividing by the length of the sequence.
#'
#' @details
#' The hydropathy values are taken from Table 2 in Kyte and Doolittle (1982).
#' 
#' @references
#' Kyte and Doolittle (1982):
#' 'A Simple Method for Displaying the  Hydropathic Character of a Protein',
#' Journal of Molecular Biology, 175, 105-132.
#' doi: 10.1016/0022-2836(82)90515-0
calculateGravyScore <- function(aa) {

    if (!is.character(aa)) stop("'aa' has to be a character vector")
    if (length(aa) != 1) stop("'aa' has to be of length 1")

    mappings <- rbind(
        c("I", "Ile", "isoleucine", 4.5),
        c("V", "Val", "valine", 4.2),
        c("L", "Leu", "leucine", 3.8),
        c("F", "Phe", "phenylalanine", 2.8),
        c("C", "Cys", "cysteine", 2.5),
        c("M", "Met", "methionine", 1.9),
        c("A", "Ala", "alanine", 1.8),
        c("G", "Gly", "glycine", -0.4),
        c("T", "Thr", "threonine", -0.7),
        c("W", "Trp", "tryptophan", -0.9),
        c("S", "Ser", "serine", -0.8),
        c("Y", "Tyr", "tyrosine", -1.3),
        c("P", "Pro", "proline", -1.6),
        c("H", "His", "histidine", -3.2),
        c("E", "Glu", "glutamic acid", -3.5),
        c("Q", "Gln", "glutamine", -3.5),
        c("D", "Asp", "aspartic acid", -3.5),
        c("N", "Asn", "asparagine", -3.5),
        c("K", "Lys", "lysine", -3.9),
        c("R", "Arg", "arginine", -4.5))

    mappings <- as.data.frame(mappings)
    colnames(mappings) <- c("SINGLE_LETTER", "THREE_LETTER", "NAME", "SCORE")
    rownames(mappings) <- mappings$SINGLE_LETTER
    mappings$SCORE <- as.numeric(mappings$SCORE)
    score <- mappings[, "SCORE", drop = FALSE]

    ## get all the characters of the sequence
    aa_split <- strsplit(aa, split = "")[[1]]

    ## get the corresponding scores for the sequence characters and calculate
    ## the mean of the sequence
    aa_score <- score[aa_split, ]
    mean(aa_score)
}


## unit tests
library(testthat)

## sample sequences taken from
## https://www.bioinformatics.org/sms2/protein_gravy.html
test_that("calculateGravyScore", {

    expect_error(calculateGravyScore(1),
        "'aa' has to be a character vector")
    expect_error(calculateGravyScore(c("MQ", "AELS")),
        "'aa' has to be of length 1")
    expect_equal(calculateGravyScore("xycAME"), NaN)
    seq <- "MQKSPLEKASFISKLFFSWTTPILRKGYRHHLELSDIYQAPSADSADHLSEKLEREWDREQASKKNPQLIHALRRCFFWRFLFYGILLYLGEVTKAVQPVLLGRIIASYDPENKVERSIAIYLGIGLCLLFIVRTLLLHPAIFGLHRIGMQMRTAMFSLIYKKTLKLSSRVLDKISIGQLVSLLSNNLNKFDEGLALAHFIWIAPLQVTLLMGLLWDLLQFSAFCGLGLLIILVIFQAILGKMMVKYRDQRAAKINERLVIT"
    expect_equal(calculateGravyScore(seq), 0.308, tolerance = 1e-05)
    seq <- "SEIIDNIYSVKAYCWESAMEKMIENLREVELKMTRKAAYMRFFTSSAFFFSGFFVVFLSVLPYTVINGIVLRKIFTTISFCIVLRMSVTRQFPTAVQIWYDSFGMIRKIQDFLQKQEYKVLEYNLMTTGI"
    expect_equal(calculateGravyScore(seq), 0.321, tolerance = 1e-03)
    seq <- "IMENVTAFWEEGFGELLQKAQQSNGDRKHSSDENNVSFSHLCLVGNPVLKNINLNIEKGEMLAITGSTGLGKTSLLMLILGELEASEGIIKHSGRVSFCSQFSWIMPGTIKENIIFGVSYDEYRYKSV"
    expect_equal(calculateGravyScore(seq), -0.127, tolerance = 1e-03)
})


#' @name calculateIsoelectricPoint
#' 
#' @title Calculate isoelectric point of an amino acid/protein
#' 
#' @description 
#' The function `calculateIsoelectricPoint` calculates the isoelectric point 
#' of a given amino acid/protein sequence. 
#' Several methods are implemented that use different pKa values for the 
#' ionizable groups of proteins. 
#'
#' @details
#' The pKa values are taken from Table 4 in Kozlowski (2016).
#' 
#' For polypeptides, the isoelectric point depends primarily on the dissociation 
#' constants (pKa) for the ionizable groups of seven charged amino acids: 
#' glutamate (\delta-carboxyl group), aspartate (\beta-carboxyl group), cysteine 
#' (thiol group), tyrosine (phenol group), histidine (imidazole side chains), 
#' lysine (\epsilon-ammonium group) and arginine (guanidinium group). Moreover, 
#' the charge of the terminal groups (NH2 and COOH) can greatly affect the pI 
#' of short peptides. Generally, the Glu, Asp, Cys, and Tyr ionizable groups 
#' are uncharged below their pKa and negatively charged above their pKa. 
#' Similarly, the His, Lys, and Arg ionizable groups are positively charged 
#' below their pKa and uncharged above their pKa.
#' 
#' @references
#' Kozlowski (2016):
#' 'IPC - Isoelectric Point Calculator',
#' Biology Direct, 11, 55.
#' doi: doi.org/10.1186/s13062-016-0159-9
calculateIsoelectricPoint <- function(aa, 
    method = c("EMBOSS", "DTASelect", "Solomon", "Sillero", "Rodwell", 
        "Lehninger", "Toseland", "Thurlkill", "Nozaki", 
        "IPC_protein", "IPC_peptide")) {
    
    method <- match.arg(method)
    
    if (method == "EMBOSS")
        pKa <- c(NH2 = 8.6, COOH = 3.6, C = 8.5, D = 3.9, E = 4.1, 
            H = 6.5, K = 10.8, R = 12.5, Y = 10.1)
    
    if (method == "DTASelect") 
        pKa <- c(NH2 = 8, COOH = 3.1, C = 8.5, D = 4.4, E = 4.4, 
            H = 6.5, K = 10, R = 12, Y = 10)
    
    if (method == "Solomon")
        pKa <- c(NH2 = 9.6, COOH = 2.4, C = 8.3, D = 3.9, E = 4.3, 
            H = 6, K = 10.5, R = 12.5, Y = 10.1)

    if (method == "Sillero")
        pKa <- c(NH2 = 8.2, COOH = 3.2, C = 9, D = 4, E = 4.5, 
            H = 6.4, K = 10.4, R = 12, Y = 10)
    
    if (method == "Rodwell")
        pKa <- c(NH2 = 8, COOH = 3.1, C = 8.33, D = 3.68, E = 4.25, 
            H = 6, K = 11.5, R = 11.5, Y = 10.07)
    
    if (method == "Lehninger")
        pKa <- c(NH2 = 9.69, COOH = 2.34, C = 8.33, D = 3.86, E = 4.25, 
            H = 6, K = 10.5, R = 12.4, Y = 10)
    
    if (method == "Toseland")
        pKa <- c(NH2 = 8.71, COOH = 3.19, C = 6.87, D = 3.6, E = 4.29, 
            H = 6.33, K = 10.45, R = 12, Y = 9.61)
    
    if (method == "Thurlkill")
        pKa <- c(NH2 = 8, COOH = 3.67, C = 8.55, D = 3.67, E = 4.25, 
            H = 6.54, K = 10.4, R = 12, Y = 9.84)
    
    if (method == "Nozaki")
        pKa <- c(NH2 = 7.5, COOH = 3.8, C = 9.5, D = 4, E = 4.4, 
            H = 6.3, K = 10.4, R = 12, Y = 9.6)
    
    if (method == "IPC_protein")
        pKa <- c(NH2 = 9.094, COOH = 2.869, C = 7.555, D = 3.872, E = 4.412, 
            H = 5.637, K = 9.052, R = 11.84, Y = 10.85)
    
    if (method == "IPC_peptide")
        pKa <- c(NH2 = 9.564, COOH = 2.383, C = 8.297, D = 3.887, E = 4.317, 
            H = 6.018, K = 10.517, R = 12.503, Y = 10.071)
    
    ## get all the characters of the sequence
    aa_split <- strsplit(aa, split = "")[[1]]
    residues <- c("C", "D", "E", "H", "K", "R", "Y")
    num <- table(aa_split)[residues]
    num <- num[!is.na(names(num))]
    num_missing <- residues[!residues %in% names(num)]
    if (length(num_missing) > 0) {
        num_missing_vec <- numeric(length(num_missing))
        names(num_missing_vec) <- num_missing 
        num <- c(num, num_missing_vec)
    }
    num <- c("NH2" = 1, "COOH" = 1, num)
    
    ## the net charge of the peptide/protein is related to the solution pH,
    ## use the Henderson-Hasselbalch equation to calculate the charge at a 
    ## certain pH
    pH <- 7
    pHprev <- 0
    pHnext <- 14
    pI <- FALSE
    
    ## implement the bisection algorithm, which in each iteration halves
    ## the search space and then moves higher or lower by 3.5 (half of 7)
    ## depending on the charge. In the next ieration, the pH is changed
    ## by 1.75 (half of 3.5), ... The process is repeated until the algorithm
    ## reaches the desired precision
    while (!pI) {
        ## for negatively charged residues
        ## the Glu (E), Asp (D), Cys (C), and Tyr (Y) ionizable groups are 
        ## uncharged below their pKa and negatively charged above their pKa,
        ## charge of COOH group: negative
        
        sum_COOH <- sum_E <- sum_D <- sum_C <- sum_Y <- 0
        
        ## calculate the sum of the Henderson-Hasselbalch equation along the
        ## residues
        sum_COOH <- -num[["COOH"]] / (1 + 10^(pKa[["COOH"]] - pH))
        sum_E <- -num[["E"]] / (1 + 10^(pKa[["E"]] - pH))
        sum_D <- -num[["D"]] / (1 + 10^(pKa[["D"]] - pH))
        sum_C <- -num[["C"]] / (1 + 10^(pKa[["C"]] - pH))
        sum_Y <- -num[["Y"]] / (1 + 10^(pKa[["Y"]] - pH))
        sum_neg <- sum_COOH + sum_E + sum_D + sum_C + sum_Y

        ## for positively charged residues
        ## His (H), Lys (K), and Arg (R) are positively charged below their pKa  
        ## and uncharged above their pKa. 
        ## charge of NH2 group: positive
        sum_NH2 <- sum_H <- sum_K <- sum_R <- 0
        
        ## calculate the sum of the Henderson-Hasselbalch equation along the 
        ## residues
        sum_NH2 <- num[["NH2"]] / (1 + 10^(pH - pKa[["NH2"]]))
        sum_H <- num[["H"]] / (1 + 10^(pH - pKa[["H"]]))
        sum_K <- num[["K"]] / (1 + 10^(pH - pKa[["K"]]))
        sum_R <- num[["R"]] / (1 + 10^(pH - pKa[["R"]]))
        sum_pos <- sum_NH2 + sum_H + sum_K + sum_R
        
        ## charge of a macromolecule at a given pH is the sum of the positive
        ## and negative charges of the individual amino acids 
        if (sum_pos + sum_neg < 0) {
            ## the new pH value must be smaller
            pH_tmp <- pH
            pH <- pH - (pH - pHprev) / 2
            pHnext <- pH_tmp
        } else {
            ## the new pH value must be higher
            pH_tmp <- pH
            pH <- pH + (pHnext - pH) / 2
            pHprev <- pH_tmp
        }
        
        if (pH - pHprev < 0.001 & pHnext - pH < 0.001) 
            pI <- TRUE
        
        if (pH > 14) stop("'pH' is greater than 14")
        if (pH < 0) stop("'pH' is smaller than 14")
    }
    
    return(pH)
}

test_that("calculateIsoelectricPoint", {
    
    expect_error(calculateIsoelectricPoint(seq, "foo"), 
        "'arg' should be one of")
    
    ## values taken from http://isoelectric.org/calculate.php
    ## P99029-1, experimental isoelectric points (different sources): 
    ## 7.84, 7.65, 7.54
    seq <- "MLQLGLRVLGCKASSVLRASTCLAGRAGRKEAGWECGGARSFSSSAVTMAPIKVGDAIPSVEVFEGEPGKKVNLAELFKGKKGVLFGVPGAFTPGCSKTHLPGFVEQAGALKAKGAQVVACLSVNDVFVIEEWGRAHQAEGKVRLLADPTGAFGKATDLL
LDDSLVSLFGNRRLKRFSMVIDNGIVKALNVEPDGTGLTCSLAPNILSQL"
    expect_equal(calculateIsoelectricPoint(seq, "EMBOSS"), 9.147, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "DTASelect"), 8.848, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Solomon"), 9.101,
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Sillero"), 9.291, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Rodwell"), 9.03,
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Lehninger"), 9.127,
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Toseland"), 8.18, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Thurlkill"), 9.02, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Nozaki"), 9.542, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "IPC_protein"), 8.054,
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "IPC_peptide"), 9.101, 
        tolerance = 5e-02)
    
    ## P40185-1 experimental isoelectric point 6.89
    seq <- "MFLRNSVLRTAPVLRRGITTLTPVSTKLAPPAAASYSQAMKANNFVYVSGQIPYTPDNKPVQGSISEKAEQVFQNVKNILAESNSSLDNIVKVNVFLADMKNFAEFNSVYAKHFHTHKPARSCVGVASLPLNVDLEMEVIAVEKN"
    expect_equal(calculateIsoelectricPoint(seq, "EMBOSS"), 9.764, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "DTASelect"), 9.269, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Solomon"), 9.694, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Sillero"), 9.542, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Rodwell"), 9.918, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Lehninger"), 9.667, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Toseland"), 9.357, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Thurlkill"), 9.443, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Nozaki"), 9.421, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "IPC_protein"), 8.64, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "IPC_peptide"), 9.687, 
        tolerance = 5e-02)
    
    ## Q9XFT3-1 experimental isoelectric point 8.07
    seq <- "MASMGGLHGASPAVLEGSLKINGSSRLNGSGRVAVAQRSRLVVRAQQSEETSRRSVIGLVAAGLAGGSFVQAVLADAISIKVGPPPAPSGGLPAGTDNSDQARDFALALKDRFYLQPLPPTEAAARAKESAKDIINVKPLIDRKAWPYVQNDLRSKASYLRYDLNTIISSKPKDEKKSLKDLTTKLFDTIDNLDYAAKKKSPSQAEKYYAETVSALNEVLAKLG"
    expect_equal(calculateIsoelectricPoint(seq, "EMBOSS"), 10.222, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "DTASelect"), 9.645, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Solomon"), 10.046, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Sillero"), 9.921, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Rodwell"), 10.498, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Lehninger"), 10.017, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Toseland"), 9.824,
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Thurlkill"), 9.866, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Nozaki"), 9.784,
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "IPC_protein"), 8.958, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "IPC_peptide"), 10.046, 
        tolerance = 5e-02)
 
    ## Q9BVP2-1 experimental isoelectric point 6.47
    seq <- "MKRPKLKKASKRMTCHKRYKIQKKVREHHRKLRKEAKKRGHKKPRKDPGVPNSAPFKEALLREAELRKQRLEELKQQQKLDRQKELEKKRKLETNPDIKPSNVEPMEKEFGLCKTENKAKSGKQNSKKLYCQELKKVIEASDVVLEVLDARDPLGCRCPQVEEAIVQSGQKKLVLILNKSDLVPKENLESWLNYLKKELPTVVFRASTKPKDKGKITKRVKAKKNAAPFRSEVCFGKEGLWKLLGGFQETCSKAIRVGVIGFPNVGKSSIINSLKQEQMCNVGVSMGLTRSMQVVPLDKQITIIDSPSFIVSPLNSSSALALRSPASIEVVKPMEAASAILSQADARQVVLKYTVPGYRNSLEFFTVLAQRRGMHQKGGIPNVEGAAKLLWSEWTGASLAYYCHPPTSWTPPPYFNESIVVDMKSGFNLEELEKNNAQSIRAIKGPHLANSILFQSSGLTNGIIEEKDIHEELPKRKERKQEEREDDKDSDQETVDEEVDENSSGMFAAEETGEALSEETTAGEQSTRSFILDKIIEEDDAYDFSTDYV"
    expect_equal(calculateIsoelectricPoint(seq, "EMBOSS"), 9.41,
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "DTASelect"), 8.933, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Solomon"), 9.249, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Sillero"), 9.304,
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Rodwell"), 9.579, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Lehninger"), 9.237, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Toseland"), 8.903, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Thurlkill"), 9.149,
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "Nozaki"), 9.363, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "IPC_protein"), 8.081, 
        tolerance = 5e-02)
    expect_equal(calculateIsoelectricPoint(seq, "IPC_peptide"), 9.251, 
        tolerance = 5e-02)
})
