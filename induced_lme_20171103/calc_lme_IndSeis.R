


library(nlme)
source("RNsubset.R")

#files <- system("ls IndSeisResid2/*.csv", intern = TRUE)
files <- system("ls resid_*.csv", intern = TRUE)
short <- c("CEUS", "NGAE_p01", "NGAE_p03", "NGAW2")

for(i in 1:length(files)){
    df1 <- read.csv(files[i])
    df1$evid <- df1$originTime
    df1[df1 == -999] <- NA
    # add evid
    df1$evid <- as.character(df1$originTime)
    imcols <- names(df1)[grep('resid_', names(df1))]
    for(j in 1:length(imcols)){
        imcol <- imcols[j]
        df1 <- df1[!is.na(df1[, imcol]), ]
        eqidcol <- which(names(df1) == 'evid')
        df <- RNsubset(df1, 4, eqidcol, 200, 5)
        m <- lme(formula(paste(imcol, "~ 1")), random = ~1 | evid, data = df,
                 control=lmeControl(opt = "optim"))
        df[,paste(imcol, "_c_ci95l",sep="")] <- as.numeric(intervals(m)$fixed)[1]
        df[,paste(imcol, "_c",sep="")] <- as.numeric(intervals(m)$fixed)[2]
        df[,paste(imcol, "_c_ci95u",sep="")] <- as.numeric(intervals(m)$fixed)[3]
        df[,paste(imcol, "_phi_ci95l")] <- as.numeric(intervals(m)$sigma)[1]
        df[,paste(imcol, "_phi")] <- as.numeric(intervals(m)$sigma)[2]
        df[,paste(imcol, "_phi_ci95u")] <- as.numeric(intervals(m)$sigma)[3]
        df[,paste(imcol, "_tau_ci95l")] <-
            as.numeric(intervals(m)$reStruct$evid)[1]
        df[,paste(imcol, "_tau")] <- as.numeric(intervals(m)$reStruct$evid)[2]
        df[,paste(imcol, "_tau_ci95u")] <-
            as.numeric(intervals(m)$reStruct$evid)[3]
        df[,paste(imcol, "_intra")]  <- as.numeric(residuals(m, level = 1))
        inter <- as.numeric(m$coefficients$random$evid)
        inter_evid <- rownames(m$coefficients$random$evid)
        mm <- match(df$evid, inter_evid)
        df[,paste(imcol, "_inter")]  <- inter[mm]
#        write.csv(df, file = paste("forMorgan/", short[i], "_", imcols[j], "_", imcol, '_lme.csv', sep=""), row.names=FALSE)
        write.csv(df, file = paste(short[i], "_", imcols[j], "_", imcol, '_lme.csv', sep=""), row.names=FALSE)
    }
}

# Redo NGA W2 for M>3.95
df1 <- read.csv("IndSeisResid2/resid_NGAW2_evDateTime_allM.csv")
df1$evid <- df1$originTime
df1[df1 == -999] <- NA
df1 <- df1[df1$mag >= 3.95, ]
# add evid
df1$evid <- as.character(df1$originTime)
imcols <- names(df1)[grep('resid_', names(df1))]
for(j in 1:length(imcols)){
    imcol <- imcols[j]
    df1 <- df1[!is.na(df1[, imcol]), ]
    eqidcol <- which(names(df1) == 'evid')
    df <- RNsubset(df1, 4, eqidcol, 200, 5)
    m <- lme(formula(paste(imcol, "~ 1")), random = ~1 | evid, data = df,
             control=lmeControl(opt = "optim"))
    df[,paste(imcol, "_c_ci95l",sep="")] <- as.numeric(intervals(m)$fixed)[1]
    df[,paste(imcol, "_c",sep="")] <- as.numeric(intervals(m)$fixed)[2]
    df[,paste(imcol, "_c_ci95u",sep="")] <- as.numeric(intervals(m)$fixed)[3]
    df[,paste(imcol, "_phi_ci95l")] <- as.numeric(intervals(m)$sigma)[1]
    df[,paste(imcol, "_phi")] <- as.numeric(intervals(m)$sigma)[2]
    df[,paste(imcol, "_phi_ci95u")] <- as.numeric(intervals(m)$sigma)[3]
    df[,paste(imcol, "_tau_ci95l")] <-
        as.numeric(intervals(m)$reStruct$evid)[1]
    df[,paste(imcol, "_tau")] <- as.numeric(intervals(m)$reStruct$evid)[2]
    df[,paste(imcol, "_tau_ci95u")] <-
        as.numeric(intervals(m)$reStruct$evid)[3]
    df[,paste(imcol, "_intra")]  <- as.numeric(residuals(m, level = 1))
    inter <- as.numeric(m$coefficients$random$evid)
    inter_evid <- rownames(m$coefficients$random$evid)
    mm <- match(df$evid, inter_evid)
    df[,paste(imcol, "_inter")]  <- inter[mm]
    write.csv(df, file = paste("forMorgan/", short[i], "_",
                               imcols[j], "_", imcol, '_lme_Mmin3.95.csv', sep=""),
              row.names=FALSE)
}





