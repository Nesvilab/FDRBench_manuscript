options(repr.matrix.max.cols=150, repr.matrix.max.rows=200)

library(readr)
library(dplyr)
library(ggplot2)
library(DBI)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(readr)
library(dplyr)
library(stringr)

# report_file: report.tsv file, pep_file: peptide level pair file
run_diann_fdp_analysis=function(report_file="",level="protein",pep_file=NULL,prefix="test",k_fold=1,pick_one_protein_method="first",out_dir=NULL,r=NULL) {
    a <- read_tsv(report_file)
    n_run <- a %>% select(Run) %>% distinct() %>% nrow()
    if(level=="protein"){
        if(n_run>=2){
            ## multiple runs
            cat("Multiple runs in the report file:",n_run,"\n")
            b <- a %>% select(`Protein.Group`,`Lib.PG.Q.Value`) %>% distinct()
            b$q_value <- b$Lib.PG.Q.Value
        }else{
            cat("Single run in the report file\n")
            b <- a %>% select(`Protein.Group`,`PG.Q.Value`) %>% distinct()
            b$q_value <- b$PG.Q.Value
        }
        #b$protein <- sapply(b$Protein.Group,function(x){ y<-  str_split(x,pattern = ";") %>% unlist; y[ length(y)]})
        b$protein <- b$Protein.Group
        set.seed(2024)
        b$score <- sample(x = 1:nrow(b),size = nrow(b),replace = FALSE)
        b <- b %>% arrange(q_value,score) %>% mutate(score=row_number())

        cat("The number of proteins:",nrow(b),"\n")

        if(is.null(out_dir)){
            out_dir <- dirname(report_file)
        }else{
            # create the directory if it does not exist
            if(!dir.exists(out_dir)){
                dir.create(out_dir)
            }
        }
        out_file <- paste(out_dir,"/",prefix,"-fdp_protein_input.tsv",sep="")
        write_tsv(b,out_file)
        fdp_file <- paste(out_dir,"/",prefix,"-diann_fdp_protein.csv",sep="")

        if(!is.null(r)){
            cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -level protein -o ",fdp_file, " -score score:0"," -r ",r," -pick ",pick_one_protein_method,sep="")
        }else{
            cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -level protein -o ",fdp_file, " -score score:0"," -fold ",k_fold," -pick ",pick_one_protein_method,sep="")
        }
        cat("Running ",cmd,"\n")
        out <- system(cmd,intern = TRUE)
        cat(paste(out,collapse = "\n"),"\n")
        return(fdp_file)
    }else if(level=="peptide" || level=="precursor"){ 
        if(n_run>=2){
            ## multiple runs
            cat("Multiple runs in the report file:",n_run,"\n")
            b <- a %>% select(`Run`,`Stripped.Sequence`,`Modified.Sequence`,`Precursor.Charge`,`Lib.Q.Value`,`PEP`,`Protein.Group`) %>% 
                rename(q_value=`Lib.Q.Value`,run=`Run`,peptide=`Stripped.Sequence`,mod_peptide=`Modified.Sequence`,charge=`Precursor.Charge`,protein=`Protein.Group`)

            ## keep the top precursor
            b <- b %>% group_by(peptide,mod_peptide,charge) %>% arrange(q_value,PEP) %>% filter(row_number()==1) %>% ungroup()
        }else{
            cat("Single run in the report file\n")
            # b <- a %>% select(`Run`,`Stripped.Sequence`,`Modified.Sequence`,`Precursor.Charge`,`Global.Q.Value`,`PEP`,`Protein.Group`) %>% distinct()
            # b <- a %>% select(`Run`,`Stripped.Sequence`,`Modified.Sequence`,`Precursor.Charge`,`Global.Q.Value`,`PEP`,`Protein.Group`) %>% 
            #   rename(q_value=`Global.Q.Value`,run=`Run`,peptide=`Stripped.Sequence`,mod_peptide=`Modified.Sequence`,charge=`Precursor.Charge`,protein=`Protein.Group`)
            b <- a %>% select(`Run`,`Stripped.Sequence`,`Modified.Sequence`,`Precursor.Charge`,`Q.Value`,`PEP`,`Protein.Group`) %>% distinct()
            b <- a %>% select(`Run`,`Stripped.Sequence`,`Modified.Sequence`,`Precursor.Charge`,`Q.Value`,`PEP`,`Protein.Group`) %>% 
                rename(q_value=`Q.Value`,run=`Run`,peptide=`Stripped.Sequence`,mod_peptide=`Modified.Sequence`,charge=`Precursor.Charge`,protein=`Protein.Group`)
        }
        set.seed(2024)
        b <- b %>% arrange(q_value,PEP) %>% mutate(score=row_number())

        cat("The number of peptides:",nrow(b),"\n")

        if(is.null(out_dir)){
            out_dir <- dirname(report_file)
        }else{
            # create the directory if it does not exist
            if(!dir.exists(out_dir)){
                dir.create(out_dir)
            }
        }
        out_file <- paste(out_dir,"/",prefix,"-fdp_precursor_input.tsv",sep="")
        write_tsv(b,out_file)
        fdp_file <- paste(out_dir,"/",prefix,"-diann_fdp_precursor.csv",sep="")

        if(!is.null(pep_file)){
            if(!is.null(r)){
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -pep ", pep_file, " -level precursor -o ",fdp_file, " -r ",r," -score score:0",sep="")
            }else{
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file," -fold ",k_fold, " -pep ", pep_file, " -level precursor -o ",fdp_file, " -score score:0",sep="")
            }
            cat("Running ",cmd,"\n")
            out <- system(cmd,intern = TRUE)
            cat(paste(out,collapse = "\n"),"\n")
            return(fdp_file)
        }else{
            cat("No paired peptide file\n")
        }

    }

}

## protein.tsv
run_fragpipe_fdp_analysis=function(report_file="",level="protein",pep_file=NULL,prefix="test",k_fold=1,pick_one_protein_method="first",out_dir=NULL,r=NULL,require_unique_peptide=FALSE) {
    a <- read_tsv(report_file)
    # n_run <- a %>% select(Run) %>% distinct() %>% nrow()
    if(require_unique_peptide){
        n_total <- nrow(a)
        a <- a %>% filter(`Unique Peptides`>=1)
        cat("Proteins without unique peptides:",n_total-nrow(a),"\n")
    }
    if(level=="protein"){
        a$protein <- ""
        for(i in 1:nrow(a)){
            if(!is.na(a$`Indistinguishable Proteins`[i])){
                a$protein[i] <- paste(a$Protein[i],a$`Indistinguishable Proteins`[i],sep=";")
            }else{
                a$protein[i] <- a$Protein[i]
            }
        }
        b <- a %>% rename(raw_protein=Protein) %>% select(protein,raw_protein,`Indistinguishable Proteins`, `Top Peptide Probability`) %>% distinct()
        b$q_value <- 0.01
        set.seed(2024)
        b <- b %>% arrange(q_value,desc(`Top Peptide Probability`)) %>% mutate(score=row_number())

        cat("The number of proteins:",nrow(b),"\n")

        if(is.null(out_dir)){
            out_dir <- dirname(report_file)
        }else{
            # create the directory if it does not exist
            if(!dir.exists(out_dir)){
                dir.create(out_dir)
            }
        }
        out_file <- paste(out_dir,"/",prefix,"-fdp_protein_input.tsv",sep="")
        write_tsv(b,out_file)
        fdp_file <- paste(out_dir,"/",prefix,"-fragpipe_fdp_protein.csv",sep="")

        if(!is.null(r)){
            cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -level protein -o ",fdp_file, " -score score:0"," -r ",r," -pick ",pick_one_protein_method,sep="")
        }else{
            cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -level protein -o ",fdp_file, " -score score:0"," -fold ",k_fold," -pick ",pick_one_protein_method,sep="")
        }
        cat("Running ",cmd,"\n")
        out <- system(cmd,intern = TRUE)
        cat(paste(out,collapse = "\n"),"\n")
        return(fdp_file)
    }else{
        stop("Error: not supported yet!\n")
    }

}


# report_file: elib file, pep_file: peptide level pair file
run_encyclopedia_fdp_analysis=function(report_file="",level="protein",pep_file=NULL,prefix="test",k_fold=1,pick_one_protein_method="first",out_dir=NULL,r=NULL) {
    mydb <- dbConnect(RSQLite::SQLite(), report_file)
    # dbListTables(mydb)
    if(level=="protein"){
        a <- dbReadTable(mydb,"proteinscores")
        ## remove decoy matches
        n_decoys <- a %>% select(ProteinAccession,QValue,MinimumPeptidePEP,IsDecoy) %>% filter(IsDecoy==1) %>% nrow 
        n_targets <- a %>% select(ProteinAccession,QValue,MinimumPeptidePEP,IsDecoy) %>% filter(IsDecoy==0) %>% nrow
        b <- a %>% select(ProteinAccession,QValue,MinimumPeptidePEP,IsDecoy) %>% filter(IsDecoy==0) %>%
            mutate(protein=ProteinAccession) %>% rename(q_value=QValue,PEP=MinimumPeptidePEP) %>%
            arrange(q_value,PEP) %>%
            mutate(score=row_number())

        if(str_detect(b$protein[1],pattern = fixed("sp|"))){
            cat("Uniprot protein accession\n")
            b$original_protein <- b$protein
            b$protein <- sapply(b$original_protein,function(x){
                y <- str_split(x,pattern = "\\|") %>% unlist
                if(length(y)>3){
                    cat(x,"\n")
                    stop("Error!")
                }
                if(length(y)>1){
                    return(y[2])
                }else{
                    return(x)
                }
            })
        }
        
        cat("The number of proteins:",nrow(b),"\n")

        if(is.null(out_dir)){
            out_dir <- dirname(report_file)
        }else{
            # create the directory if it does not exist
            if(!dir.exists(out_dir)){
                dir.create(out_dir)
            }
        }
        # out_dir <- dirname(report_file)
        out_file <- paste(out_dir,"/",prefix,"-fdp_protein_input.tsv",sep="")
        write_tsv(b,out_file)
        fdp_file <- paste(out_dir,"/",prefix,"-encyclopedia_fdp_protein.csv",sep="")

        #   cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -level protein -o ",fdp_file, " -score score:0",sep="")
        if(!is.null(r)){
            cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -level protein -o ",fdp_file, " -score score:0"," -r ",r," -pick ",pick_one_protein_method,sep="")
        }else{
            cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -level protein -o ",fdp_file, " -score score:0"," -fold ",k_fold," -pick ",pick_one_protein_method,sep="")
        }

        cat("Running ",cmd,"\n")
        out <- system(cmd,intern = TRUE)
        cat(paste(out,collapse = "\n"),"\n")
        return(fdp_file)
    }else if(level=="peptide"){ 
        a <- dbReadTable(mydb,"peptidescores")
        # pep2pro <- dbReadTable(mydb,"peptidetoprotein")
        n_targets <- a %>% select(PeptideSeq,QValue,PosteriorErrorProbability,IsDecoy) %>% filter(IsDecoy==0)
        n_decoys <- a %>% select(PeptideSeq,QValue,PosteriorErrorProbability,IsDecoy) %>% filter(IsDecoy==1)
        b <- a %>% select(PeptideSeq,QValue,PosteriorErrorProbability,IsDecoy) %>% filter(IsDecoy==0) %>%
            mutate(peptide=PeptideSeq) %>% rename(q_value=QValue,PEP=PosteriorErrorProbability) %>%
            arrange(q_value,PEP) %>%
            mutate(score=row_number())

        cat("The number of peptides:",nrow(b),"\n")

        if(is.null(out_dir)){
            out_dir <- dirname(report_file)
        }else{
            # create the directory if it does not exist
            if(!dir.exists(out_dir)){
                dir.create(out_dir)
            }
        }
        # out_dir <- dirname(report_file)
        out_file <- paste(out_dir,"/",prefix,"-fdp_peptide_input.tsv",sep="")
        write_tsv(b,out_file)
        fdp_file <- paste(out_dir,"/",prefix,"-encyclopedia_fdp_peptide.csv",sep="")

        if(!is.null(pep_file)){
            #   cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -pep ", pep_file, " -level peptide -o ",fdp_file, " -score score:0",sep="")
            if(!is.null(r)){
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -pep ", pep_file, " -level peptide -o ",fdp_file, " -r ",r," -score score:0",sep="")
            }else{
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -pep ", pep_file, " -level peptide -o ",fdp_file, " -score score:0",sep="")
            }
            cat("Running ",cmd,"\n")
            out <- system(cmd,intern = TRUE)
            cat(paste(out,collapse = "\n"),"\n")
            return(fdp_file)
        }else{
            cat("No paired peptide file\n")
        }

    }

}



run_spectronaut_fdp_analysis=function(report_file="",level="protein",pep_file=NULL,prefix="test",k_fold=1,pick_one_protein_method="first",out_dir=NULL,r=NULL) {
    a <- read_tsv(report_file)
    if(level=="protein"){
        ## for protein level report
        ## [1] "R.Condition"            "R.FileName"             "R.Replicate"            "PG.ProteinAccessions"   "PG.ProteinDescriptions" "PG.ProteinNames"       
        ## [7] "PG.Coverage"            "PG.IsSingleHit"         "PG.Qvalue"              "PG.Quantity"  
        if(str_detect(report_file,fixed("_Report_BGS Factory Report (Normal).tsv"))){
            b <- a %>% filter(!is.na(`PG.Qvalue`)) %>% select(`PG.ProteinGroups`,`PG.Qvalue`,`PG.PEP`) %>%
                    distinct() %>%
                    rename(q_value=PG.Qvalue,PEP=PG.PEP,protein=`PG.ProteinGroups`) %>%
                    arrange(q_value,PEP) %>%
                    mutate(score=row_number())
        }else{
            if("PG.PEP" %in% names(a)){
            b <- a %>% filter(!is.na(`PG.Qvalue`),!is.na(`PG.Quantity`)) %>% select(`PG.ProteinAccessions`,`PG.Qvalue`,`PG.PEP`) %>%
                rename(q_value=PG.Qvalue,PEP=PG.PEP,protein=`PG.ProteinAccessions`) %>%
                arrange(q_value,PEP) %>%
                mutate(score=row_number())
            }else{
                b <- a %>% filter(!is.na(`PG.Qvalue`)) %>% select(`PG.ProteinAccessions`,`PG.Qvalue`) %>%
                rename(q_value=PG.Qvalue,protein=`PG.ProteinAccessions`) %>%
                arrange(q_value) %>%
                mutate(score=row_number())
            }
        }
        
        cat("The number of proteins:",nrow(b),"\n")

        if(is.null(out_dir)){
            out_dir <- dirname(report_file)
        }else{
            # create the directory if it does not exist
            if(!dir.exists(out_dir)){
                dir.create(out_dir)
            }
        }

        #out_dir <- dirname(report_file)
        out_file <- paste(out_dir,"/",prefix,"-fdp_protein_input.tsv",sep="")
        write_tsv(b,out_file)
        fdp_file <- paste(out_dir,"/",prefix,"-spectronaut_fdp_protein.csv",sep="")

        if(!is.null(r)){
            cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -level protein -o ",fdp_file, " -score score:0"," -r ",r," -pick ",pick_one_protein_method,sep="")
        }else{
            cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -level protein -o ",fdp_file, " -score score:0"," -fold ",k_fold," -pick ",pick_one_protein_method,sep="")
        }

        cat("Running ",cmd,"\n")
        out <- system(cmd,intern = TRUE)
        cat(paste(out,collapse = "\n"),"\n")
        return(fdp_file)

    }else if(level=="peptide" || level=="precursor"){ 
        ## for protein level report
        ## [1] "R.Condition"            "R.FileName"             "R.Replicate"            "PG.ProteinAccessions"   "PG.ProteinDescriptions" "PG.ProteinNames"       
        ## [7] "PG.Coverage"            "PG.IsSingleHit"         "PG.Qvalue"              "PG.Quantity"    
        a <- read_tsv(report_file)
        if("EG.PEP" %in% names(a)){
            b <- a %>% filter(!is.na(`EG.Qvalue`)) %>% select(`EG.ModifiedSequence`,`PEP.StrippedSequence`,`EG.Qvalue`,`EG.PEP`,`FG.Charge`,`PG.ProteinAccessions`) %>%
                rename(q_value=EG.Qvalue,PEP=EG.PEP,mod_peptide=`EG.ModifiedSequence`,peptide=`PEP.StrippedSequence`,charge=`FG.Charge`,protein=`PG.ProteinAccessions`) %>%
                arrange(q_value,PEP) %>%
                mutate(score=row_number())
        }else{
            b <- a %>% filter(!is.na(`EG.Qvalue`)) %>% select(`EG.ModifiedSequence`,`PEP.StrippedSequence`,`EG.Qvalue`,`FG.Charge`,`PG.ProteinAccessions`) %>%
                rename(q_value=EG.Qvalue,mod_peptide=`EG.ModifiedSequence`,peptide=`PEP.StrippedSequence`,charge=`FG.Charge`,protein=`PG.ProteinAccessions`) %>%
                arrange(q_value) %>%
                mutate(score=row_number())
        }

        if(is.null(out_dir)){
            out_dir <- dirname(report_file)
        }else{
            # create the directory if it does not exist
            if(!dir.exists(out_dir)){
                dir.create(out_dir)
            }
        }

        # out_dir <- dirname(report_file)
        out_file <- paste(out_dir,"/",prefix,"-fdp_precursor_input.tsv",sep="")
        write_tsv(b,out_file)
        fdp_file <- paste(out_dir,"/",prefix,"-spectronaut_fdp_precursor.csv",sep="")

        if(!is.null(pep_file)){
            if(!is.null(r)){
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -pep ", pep_file, " -level precursor -o ",fdp_file, " -r ",r," -score score:0",sep="")
            }else{
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -pep ", pep_file, " -level precursor -o ",fdp_file, " -score score:0",sep="")
            }
            cat("Running ",cmd,"\n")
            out <- system(cmd,intern = TRUE)
            cat(paste(out,collapse = "\n"),"\n")
            return(fdp_file)
        }else{
            cat("No paired peptide file\n")
            if(!is.null(r)){
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -level precursor -o ",fdp_file, " -r ",r," -score score:0",sep="")
            }else{
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -level precursor -o ",fdp_file, " -score score:0",sep="")
            }
            cat("Running ",cmd,"\n")
            out <- system(cmd,intern = TRUE)
            cat(paste(out,collapse = "\n"),"\n")
            return(fdp_file)
        }

    }
}


run_maxquant_fdp_analysis=function(report_file="",level="protein",pep_file=NULL,prefix="test",k_fold=1,pick_one_protein_method="first",out_dir=NULL,r=NULL) {
    a <- read_tsv(report_file)
    if(level=="protein"){
        ## remove reverse proteins
        b <- a %>% filter(is.na(Reverse)) %>% mutate(protein=`Protein IDs`,q_value=`Q-value`) %>% arrange(q_value,desc(Score)) %>% mutate(score=row_number()) %>%
            select(protein,q_value,score)
        
        cat("The number of proteins:",nrow(b),"\n")

        if(is.null(out_dir)){
            out_dir <- dirname(report_file)
        }else{
            # create the directory if it does not exist
            if(!dir.exists(out_dir)){
                dir.create(out_dir)
            }
        }

        #out_dir <- dirname(report_file)
        out_file <- paste(out_dir,"/",prefix,"-fdp_protein_input.tsv",sep="")
        write_tsv(b,out_file)
        fdp_file <- paste(out_dir,"/",prefix,"-maxquant_fdp_protein.csv",sep="")

        if(!is.null(r)){
            cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -level protein -o ",fdp_file, " -score score:0"," -r ",r," -pick ",pick_one_protein_method,sep="")
        }else{
            cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -level protein -o ",fdp_file, " -score score:0"," -fold ",k_fold," -pick ",pick_one_protein_method,sep="")
        }

        cat("Running ",cmd,"\n")
        out <- system(cmd,intern = TRUE)
        cat(paste(out,collapse = "\n"),"\n")
        return(fdp_file)

    }
}


# report_file: mokapot.peptides.txt
run_mokapot_fdp_analysis=function(report_file="",level="peptide",pep_file=NULL,prefix="test",k_fold=1,pick_one_protein_method="first",out_dir=NULL,r=NULL) {
    a <- read_tsv(report_file)
    # n_run <- a %>% select(Run) %>% distinct() %>% nrow()
    if(level=="protein"){
        # TODO: need to update this part
        if(n_run>=2){
            ## multiple runs
            cat("Multiple runs in the report file:",n_run,"\n")
            b <- a %>% select(`Protein.Group`,`Lib.PG.Q.Value`) %>% distinct()
            b$q_value <- b$Lib.PG.Q.Value
        }else{
            cat("Single run in the report file\n")
            b <- a %>% select(`Protein.Group`,`PG.Q.Value`) %>% distinct()
            b$q_value <- b$PG.Q.Value
        }
        #b$protein <- sapply(b$Protein.Group,function(x){ y<-  str_split(x,pattern = ";") %>% unlist; y[ length(y)]})
        b$protein <- b$Protein.Group
        set.seed(2024)
        b$score <- sample(x = 1:nrow(b),size = nrow(b),replace = FALSE)
        b <- b %>% arrange(q_value,score) %>% mutate(score=row_number())

        cat("The number of proteins:",nrow(b),"\n")

        if(is.null(out_dir)){
            out_dir <- dirname(report_file)
        }else{
            # create the directory if it does not exist
            if(!dir.exists(out_dir)){
                dir.create(out_dir)
            }
        }
        out_file <- paste(out_dir,"/",prefix,"-fdp_protein_input.tsv",sep="")
        write_tsv(b,out_file)
        fdp_file <- paste(out_dir,"/",prefix,"-diann_fdp_protein.csv",sep="")

        if(!is.null(r)){
            cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -level protein -o ",fdp_file, " -score score:0"," -r ",r," -pick ",pick_one_protein_method,sep="")
        }else{
            cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -level protein -o ",fdp_file, " -score score:0"," -fold ",k_fold," -pick ",pick_one_protein_method,sep="")
        }
        cat("Running ",cmd,"\n")
        out <- system(cmd,intern = TRUE)
        cat(paste(out,collapse = "\n"),"\n")
        return(fdp_file)
    }else if(level=="peptide" || level=="precursor"){ 
        b <- a %>% select(Peptide,`mokapot q-value`,`mokapot PEP`,`Proteins`) %>% 
                mutate(mod_peptide=Peptide) %>%
                rename(q_value=`mokapot q-value`,peptide=Peptide,protein=`Proteins`,PEP=`mokapot PEP`) %>%
                mutate(peptide=str_replace_all(peptide,pattern = "[^A-Z]",replacement = ""))

        set.seed(2024)
        b <- b %>% arrange(q_value,PEP) %>% mutate(score=row_number())

        cat("The number of peptides:",nrow(b),"\n")
        cat("The number of peptides passed 1% FDR:",nrow(b %>% filter(q_value<=0.01)),"\n")


        if(is.null(out_dir)){
            out_dir <- dirname(report_file)
        }else{
            # create the directory if it does not exist
            if(!dir.exists(out_dir)){
                dir.create(out_dir)
            }
        }
        out_file <- paste(out_dir,"/",prefix,"-fdp_peptide_input.tsv",sep="")
        write_tsv(b,out_file)
        fdp_file <- paste(out_dir,"/",prefix,"-mokapot_fdp_peptide.csv",sep="")

        if(!is.null(pep_file)){
            if(!is.null(r)){
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -pep ", pep_file, " -level peptide -o ",fdp_file, " -r ",r," -score score:0",sep="")
            }else{
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file," -fold ",k_fold, " -pep ", pep_file, " -level peptide -o ",fdp_file, " -score score:0",sep="")
            }
            cat("Running ",cmd,"\n")
            out <- system(cmd,intern = TRUE)
            cat(paste(out,collapse = "\n"),"\n")
            return(fdp_file)
        }else{
            cat("No paired peptide file\n")
        }

    }

}

# report_file: mokapot.peptides.txt
run_percolator_fdp_analysis=function(report_file="",level="peptide",pep_file=NULL,prefix="test",k_fold=1,pick_one_protein_method="first",out_dir=NULL,r=NULL) {
    a <- read_tsv(report_file)
    # n_run <- a %>% select(Run) %>% distinct() %>% nrow()
    if(level=="protein"){
        # TODO: need to update this part
    }else if(level=="peptide" || level=="precursor"){
        if (level=="precursor") {
            b <- a %>% mutate(charge=str_sub(peptide, -3, -3)) %>% select(peptide,`q-value`,`posterior_error_prob`,`proteinIds`,charge) %>%
              mutate(mod_peptide=peptide) %>%
              rename(q_value=`q-value`,peptide=peptide,protein=`proteinIds`,PEP=`posterior_error_prob`) %>%
              # -.YNHDLALLELDEPLVLNSYVTPLC[57.0215]LADK3.-
              mutate(peptide=str_replace_all(peptide,pattern = "^..",replacement = "")) %>%
              mutate(peptide=str_replace_all(peptide,pattern = "...$",replacement = "")) %>%
              mutate(peptide=str_replace_all(peptide,pattern = "[^A-Z]",replacement = ""))
        } else {
            b <- a %>% select(peptide,`q-value`,`posterior_error_prob`,`proteinIds`) %>%
              mutate(mod_peptide=peptide) %>%
              rename(q_value=`q-value`,peptide=peptide,protein=`proteinIds`,PEP=`posterior_error_prob`) %>%
              # -.GMGGHGYGGAGDASSGFHGGHFVHMR.-
              mutate(peptide=str_replace_all(peptide,pattern = "^..",replacement = "")) %>%
              mutate(peptide=str_replace_all(peptide,pattern = "..$",replacement = "")) %>%
              mutate(peptide=str_replace_all(peptide,pattern = "[^A-Z]",replacement = ""))
        }

        set.seed(2024)
        b <- b %>% arrange(q_value,PEP) %>% mutate(score=row_number())

        cat("The number of peptides:",nrow(b),"\n")
        cat("The number of peptides passed 1% FDR:",nrow(b %>% filter(q_value<=0.01)),"\n")


        if(is.null(out_dir)){
            out_dir <- dirname(report_file)
        }else{
            # create the directory if it does not exist
            if(!dir.exists(out_dir)){
                dir.create(out_dir)
            }
        }
        out_file <- paste(out_dir,"/",prefix,"-fdp_peptide_input.tsv",sep="")
        write_tsv(b,out_file)
        fdp_file <- paste(out_dir,"/",prefix,"-percolator_fdp_peptide.csv",sep="")

        if(!is.null(pep_file)){
            if(!is.null(r)){
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -pep ", pep_file, " -level ", level, " -o ",fdp_file, " -r ",r," -score score:0",sep="")
            }else{
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file," -fold ",k_fold, " -pep ", pep_file, " -level ", level, " -o ",fdp_file, " -score score:0",sep="")
            }
            cat("Running ",cmd,"\n")
            out <- system(cmd,intern = TRUE)
            cat(paste(out,collapse = "\n"),"\n")
            return(fdp_file)
        }else{
            if(!is.null(r)){
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -level ", level, " -o ",fdp_file, " -r ",r," -score score:0",sep="")
            }else{
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file," -fold ",k_fold, " -level ", level, " -o ",fdp_file, " -score score:0",sep="")
            }
            cat("Running ",cmd,"\n")
            out <- system(cmd,intern = TRUE)
            cat(paste(out,collapse = "\n"),"\n")
            return(fdp_file)
        }

    }

}


run_percolator_reset_fdp_analysis=function(report_file="",level="peptide",pep_file=NULL,prefix="test",k_fold=1,pick_one_protein_method="first",out_dir=NULL,r=NULL) {
    a <- read_tsv(report_file) %>% filter(Label==1)
    # n_run <- a %>% select(Run) %>% distinct() %>% nrow()
    if(level=="protein"){
        # TODO: need to update this part
    }else if(level=="peptide" || level=="precursor"){ 
        b <- a %>% select(Peptide,q_val,SVM_score,Proteins) %>% 
                distinct() %>%
                mutate(mod_peptide=Peptide) %>%
                rename(q_value=q_val,peptide=Peptide,protein=`Proteins`,score=SVM_score) # %>%
                #mutate(peptide=str_replace_all(peptide,pattern = "^..",replacement = "")) %>%
                #mutate(peptide=str_replace_all(peptide,pattern = "..$",replacement = "")) %>%
                #mutate(peptide=str_replace_all(peptide,pattern = "[^A-Z]",replacement = ""))

        set.seed(2024)
        b <- b %>% arrange(q_value,desc(score)) %>% mutate(score=row_number())

        cat("The number of peptides:",nrow(b),"\n")
        cat("The number of peptides passed 1% FDR:",nrow(b %>% filter(q_value<=0.01)),"\n")


        if(is.null(out_dir)){
            out_dir <- dirname(report_file)
        }else{
            # create the directory if it does not exist
            if(!dir.exists(out_dir)){
                dir.create(out_dir)
            }
        }
        out_file <- paste(out_dir,"/",prefix,"-fdp_peptide_input.tsv",sep="")
        write_tsv(b,out_file)
        fdp_file <- paste(out_dir,"/",prefix,"-percolator_fdp_peptide.csv",sep="")

        if(!is.null(pep_file)){
            if(!is.null(r)){
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -pep ", pep_file, " -level peptide -o ",fdp_file, " -r ",r," -score score:0",sep="")
            }else{
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file," -fold ",k_fold, " -pep ", pep_file, " -level peptide -o ",fdp_file, " -score score:0",sep="")
            }
            cat("Running ",cmd,"\n")
            out <- system(cmd,intern = TRUE)
            cat(paste(out,collapse = "\n"),"\n")
            return(fdp_file)
        }else{
            cat("No paired peptide file\n")
        }

    }

}


# report_file: mokapot.peptides.txt
run_sage_fdp_analysis=function(report_file="",level="peptide",pep_file=NULL,prefix="test",k_fold=1,pick_one_protein_method="first",out_dir=NULL,r=NULL) {
    a <- read_tsv(report_file)
    
    if(level=="peptide"){ 
        b <- a %>% filter(!str_detect(proteins,"rev_")) %>% 
                select(peptide,peptide_q,posterior_error,proteins) %>% 
                arrange(peptide_q,posterior_error) %>%
                group_by(peptide) %>%
                filter(row_number()==1) %>% 
                mutate(mod_peptide=peptide) %>%
                mutate(peptide=str_replace_all(peptide,pattern = "[^A-Z]",replacement = "")) %>%
                rename(q_value=peptide_q,protein=proteins,PEP=posterior_error) %>%
                ungroup()
        set.seed(2024)
        b <- b %>% arrange(q_value,PEP) %>% mutate(score=row_number())

        cat("The number of peptides:",nrow(b),"\n")
        cat("The number of peptides passed 1% FDR:",nrow(b %>% filter(q_value<=0.01)),"\n")


        if(is.null(out_dir)){
            out_dir <- dirname(report_file)
        }else{
            # create the directory if it does not exist
            if(!dir.exists(out_dir)){
                dir.create(out_dir)
            }
        }
        out_file <- paste(out_dir,"/",prefix,"-fdp_peptide_input.tsv",sep="")
        write_tsv(b,out_file)
        fdp_file <- paste(out_dir,"/",prefix,"-sage_fdp_peptide.csv",sep="")

        if(!is.null(pep_file)){
            if(!is.null(r)){
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -pep ", pep_file, " -level peptide -o ",fdp_file, " -r ",r," -score score:0",sep="")
            }else{
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file," -fold ",k_fold, " -pep ", pep_file, " -level peptide -o ",fdp_file, " -score score:0",sep="")
            }
            cat("Running ",cmd,"\n")
            out <- system(cmd,intern = TRUE)
            cat(paste(out,collapse = "\n"),"\n")
            return(fdp_file)
        }else{
            cat("No paired peptide file\n")
        }

    }

}


run_msgf_fdp_analysis=function(report_file="",level="peptide",pep_file=NULL,prefix="test",k_fold=1,pick_one_protein_method="first",out_dir=NULL,r=NULL) {
    a <- read_tsv(report_file)
    
    if(level=="peptide"){ 
        b <- a %>% select(Peptide,Protein,SpecEValue,EValue,PepQValue) %>% 
                arrange(PepQValue,SpecEValue) %>%
                group_by(Peptide) %>%
                filter(row_number()==1) %>%
                ungroup() %>%
                mutate(mod_peptide=Peptide) %>%
                mutate(peptide=str_replace_all(Peptide,pattern = "^..",replacement = "")) %>%
                mutate(peptide=str_replace_all(peptide,pattern = "..$",replacement = "")) %>%
                mutate(peptide=str_replace_all(peptide,pattern = "[^A-Z]",replacement = "")) %>%
                rename(q_value=PepQValue,protein=Protein,PEP=SpecEValue) %>%
                select(peptide,mod_peptide,protein,q_value, PEP) %>%
                distinct()
                
        set.seed(2024)
        b <- b %>% arrange(q_value,PEP) %>% mutate(score=row_number())

        cat("The number of peptides:",nrow(b),"\n")
        cat("The number of peptides passed 1% FDR:",nrow(b %>% filter(q_value<=0.01)),"\n")


        if(is.null(out_dir)){
            out_dir <- dirname(report_file)
        }else{
            # create the directory if it does not exist
            if(!dir.exists(out_dir)){
                dir.create(out_dir)
            }
        }
        out_file <- paste(out_dir,"/",prefix,"-fdp_peptide_input.tsv",sep="")
        write_tsv(b,out_file)
        fdp_file <- paste(out_dir,"/",prefix,"-msgf_fdp_peptide.csv",sep="")

        if(!is.null(pep_file)){
            if(!is.null(r)){
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -pep ", pep_file, " -level peptide -o ",fdp_file, " -r ",r," -score score:0",sep="")
            }else{
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file," -fold ",k_fold, " -pep ", pep_file, " -level peptide -o ",fdp_file, " -score score:0",sep="")
            }
            cat("Running ",cmd,"\n")
            out <- system(cmd,intern = TRUE)
            cat(paste(out,collapse = "\n"),"\n")
            return(fdp_file)
        }else{
            cat("No paired peptide file\n")
        }

    }

}


## only peptide level is supported
run_fdp_analysis=function(report_file="",level="peptide",pep_file=NULL,prefix="test",k_fold=1,pick_one_protein_method="first",out_dir=NULL,r=NULL) {
    a <- read_tsv(report_file)
    if(level=="peptide"){         
        set.seed(2024)
        b <- a %>% arrange(q_value,score) %>% mutate(score=row_number())
        cat("The number of peptides:",nrow(b),"\n")
        cat("The number of peptides passed 1% FDR:",nrow(b %>% filter(q_value<=0.01)),"\n")
        if(is.null(out_dir)){
            out_dir <- dirname(report_file)
        }else{
            # create the directory if it does not exist
            if(!dir.exists(out_dir)){
                dir.create(out_dir)
            }
        }
        out_file <- paste(out_dir,"/",prefix,"-fdp_peptide_input.tsv",sep="")
        write_tsv(b,out_file)
        fdp_file <- paste(out_dir,"/",prefix,"-std_fdp_peptide.csv",sep="")

        if(!is.null(pep_file)){
            if(!is.null(r)){
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file, " -pep ", pep_file, " -level peptide -o ",fdp_file, " -r ",r," -score score:0",sep="")
            }else{
                cmd <- paste("java -jar ~/github/FDRBench/target/fdrbench-0.0.1/fdrbench-0.0.1.jar -i ", out_file," -fold ",k_fold, " -pep ", pep_file, " -level peptide -o ",fdp_file, " -score score:0",sep="")
            }
            cat("Running ",cmd,"\n")
            out <- system(cmd,intern = TRUE)
            cat(paste(out,collapse = "\n"),"\n")
            return(fdp_file)
        }else{
            cat("No paired peptide file\n")
        }

    }

}


plot_fdp_fdr=function(fdp_file="",fdr_max=NULL,fig_title=NULL,scale_xy=TRUE,add_numbers=FALSE) {
    x <- read_csv(fdp_file)
    if("FDP_1B" %in% names(x)){
        dat <- x %>% mutate(FDP_min=n_p/(n_p+n_t)) %>% select(q_value,FDP,FDP_1B,FDP_min) %>% distinct() %>% 
            rename(`Combined entrapment`=FDP,`Paired entrapment`=FDP_1B) %>%
            gather(key = "Method",value = "FDP",-`q_value`) %>% select(q_value,FDP,Method)
    }else{
        dat <- x %>% mutate(FDP_min=n_p/(n_p+n_t)) %>% select(q_value,FDP,FDP_min) %>% distinct() %>% 
            rename(`Combined entrapment`=FDP) %>%
            gather(key = "Method",value = "FDP",-`q_value`) %>% select(q_value,FDP,Method)
    }
    

    max_fdp <- max(c(dat$FDP,dat$q_value))
    if(!is.null(fdr_max)){
        max_fdp <- min(c(fdr_max,max_fdp))
    }
    gg1 <- ggplot(dat,aes(x=q_value,y=FDP,color=Method)) + 
            geom_abline(slope = 1,intercept = 0,color="gray")+
            #rasterise(geom_line(), dpi = 300) + 
            geom_line()+
            xlab("FDR")+
            ylab("FDP")+
            theme_bw()+
            #geom_segment(x=0.01,xend=0.01,y=0,yend=0.01,linewidth=0.3,color="blue",linetype=2)+
            theme_pubr(base_size = 12,border = TRUE)

    if(scale_xy){
        gg1 <- gg1 + geom_vline(xintercept = 0.01,linetype=2,color="blue")+
            xlim(0,max_fdp)+
            ylim(0,max_fdp)+
            scale_y_continuous(labels = scales::percent,limits =c(0,max_fdp))+
            #scale_y_pct()+
            #scale_x_pct()+
            scale_x_continuous(labels = scales::percent,limits =c(0,max_fdp))
    }
    #theme(legend.position = "top",plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"))+    
    ## legend on the botton right
    ## no background color for the legend
    gg1 <- gg1 + theme(legend.position = c(0.65, 0.16),plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"),legend.background = element_blank())
            
    if(!is.null(fig_title)){
        gg1 <- gg1 + ggtitle(fig_title)
    }

    if(add_numbers){
        # add numbers on the top left size of the figure using annotation, text align to left
        y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.4f%%",FDP001*100))
        #n_t <- 
        if(abs(max_fdp-0.01)<=0.02){
            # text right align
            gg1 <- gg1 + annotate("text", x = max_fdp*0.1, y = 0.9*max_fdp, label = paste("Total discoveries:",nrow(x %>% filter(q_value<=0.01)),"\n",paste(y$Method,y$ratio,sep=":",collapse = "\n"),sep=""), color = "black", size = 3,hjust = 0)
        }else{
            gg1 <- gg1 + annotate("text", x = 0.01, y = 0.9*max_fdp, label = paste("Total discoveries:",nrow(x %>% filter(q_value<=0.01)),"\n",paste(y$Method,y$ratio,sep=":",collapse = "\n"),sep=""), color = "black", size = 3,hjust = 0)
        }
    }
   
    # library(ggpubr)
    # options(repr.plot.width = 12, repr.plot.height = 6)
    # gg <- ggarrange(gg1,gg2,ncol = 2,  common.legend = TRUE)
    # options(jupyter.plot_mimetypes = "image/png")
    return(gg1)

    #pdf("fdp_fdr_lib_qvalue_protein.pdf",width = 8,height = 4)
    #print(gg)
    #dev.off()
}




plot_fdp_fdr_v2=function(fdp_file="",fdr_max=NULL,fig_title=NULL,scale_xy=TRUE,add_numbers=FALSE,r=1,fixed_fdr_max=FALSE,max_x=NA,max_y=NA,
                         color_mapping=NULL,
                         legend_position=NULL,
                         fdr_decimal_place=2,
                         return_data=FALSE,
                         text_position=NULL,
                         add_max_qvalue=FALSE) {
    x <- read_csv(fdp_file)
    if("FDP_1B" %in% names(x)){
        if(r>=2){
            dat <- x %>% mutate(FDP_min=n_p/(n_p+n_t)) %>% select(q_value,FDP,FDP_1B,FDP_min) %>% distinct() %>% 
                rename(`Combined method`=FDP,`Matched method`=FDP_1B,`Lower bound`=FDP_min) %>%
                gather(key = "Method",value = "FDP",-`q_value`) %>% select(q_value,FDP,Method)
        }else{
            dat <- x %>% mutate(FDP_min=n_p/(n_p+n_t)) %>% select(q_value,FDP,FDP_1B,FDP_min) %>% distinct() %>% 
                rename(`Combined method`=FDP,`Paired method`=FDP_1B,`Lower bound`=FDP_min) %>%
                gather(key = "Method",value = "FDP",-`q_value`) %>% select(q_value,FDP,Method)
        }
    }else{
        dat <- x %>% mutate(FDP_min=n_p/(n_p+n_t)) %>% select(q_value,FDP,FDP_min) %>% distinct() %>% 
            rename(`Combined method`=FDP,`Lower bound`=FDP_min) %>%
            gather(key = "Method",value = "FDP",-`q_value`) %>% select(q_value,FDP,Method)
    }
    

    max_fdp <- max(c(dat$FDP,dat$q_value))
    if(!is.null(fdr_max)){
        if(fixed_fdr_max){
            max_fdp <- fdr_max
        }else{
            max_fdp <- min(c(fdr_max,max_fdp))
        }
    }
    
    gg1 <- ggplot(dat,aes(x=q_value,y=FDP,color=Method)) + 
            geom_abline(slope = 1,intercept = 0,color="gray")+
            #rasterise(geom_line(), dpi = 300) + 
            geom_line()+
            xlab("FDR threshold")+
            ylab("Estimated FDP")+
            theme_bw()+
            #geom_segment(x=0.01,xend=0.01,y=0,yend=0.01,linewidth=0.3,color="blue",linetype=2)+
            theme_pubr(base_size = 12,border = TRUE)

    if(!is.null(color_mapping)){
        gg1 <- gg1 + scale_color_manual(values = color_mapping)
    }

    if(scale_xy){

        if(!is.na(max_x) || !is.na(max_y)){
            gg1 <- gg1 + geom_vline(xintercept = 0.01,linetype=2,color="blue")
            if(!is.na(max_x)){
                gg1 <- gg1 + xlim(0,max_x) + scale_x_continuous(labels = scales::percent,limits =c(0,max_x))
            }else{
                gg1 <- gg1 + xlim(0,max_fdp)+ scale_x_continuous(labels = scales::percent,limits =c(0,max_fdp))
            }
            if(!is.na(max_y)){
                gg1 <- gg1 + ylim(0,max_y) + scale_y_continuous(labels = scales::percent,limits =c(0,max_y))
            }else{
                gg1 <- gg1 + ylim(0,max_fdp) + scale_y_continuous(labels = scales::percent,limits =c(0,max_fdp))
            }
        }else{
            gg1 <- gg1 + geom_vline(xintercept = 0.01,linetype=2,color="blue")+
                xlim(0,max_fdp)+
                ylim(0,max_fdp)+
                scale_y_continuous(labels = scales::percent,limits =c(0,max_fdp))+
                #scale_y_pct()+
                #scale_x_pct()+
                scale_x_continuous(labels = scales::percent,limits =c(0,max_fdp))
        }
    }
    #theme(legend.position = "top",plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"))+    
    ## legend on the botton right
    ## no background color for the legend
    if(is.null(legend_position)){
        gg1 <- gg1 + theme(legend.position = c(0.65, 0.16),plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"),legend.background = element_blank(),legend.text=element_text(size=12),legend.title=element_text(size=12))    
    }else{
        gg1 <- gg1 + theme(legend.position = legend_position,plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"),legend.background = element_blank(),legend.text=element_text(size=12),legend.title=element_text(size=12))    
    }
    
            
    if(!is.null(fig_title)){
        gg1 <- gg1 + ggtitle(fig_title)
    }

    added_numbers <- NULL
    if(add_numbers){
        # add numbers on the top left size of the figure using annotation, text align to left
        # y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.4f%%",FDP001*100))
        if(fdr_decimal_place==1){
            y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% arrange(desc(q_value)) %>% filter(row_number()==1) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.1f%%",FDP001*100))
        }else if(fdr_decimal_place==2){
            y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% arrange(desc(q_value)) %>% filter(row_number()==1) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.2f%%",FDP001*100))
        }else{
            y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% arrange(desc(q_value)) %>% filter(row_number()==1) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.4f%%",FDP001*100))
        }
        #n_t <- 
        if(abs(max_fdp-0.01)<=0.02){
            # text right align
            added_numbers <- paste("Total discoveries:",nrow(x %>% filter(q_value<=0.01)),"\n",paste(y$Method,y$ratio,sep=":",collapse = "\n"),sep="")
            if(add_max_qvalue){
                added_numbers <- paste(added_numbers,"\n","Max q-value:",sprintf("%.2e",max(x$q_value)),sep="")
            }
            if(is.null(text_position)){
                gg1 <- gg1 + annotate("text", x = max_fdp*0.1, y = 0.9*max_fdp, label = added_numbers, color = "black", size = 3,hjust = 0)
            }else{
                gg1 <- gg1 + annotate("text", x = text_position[1], y = text_position[2], label = added_numbers, color = "black", size = 3,hjust = 0)
            }
            
        }else{
            added_numbers <- paste("Total discoveries:",nrow(x %>% filter(q_value<=0.01)),"\n",paste(y$Method,y$ratio,sep=":",collapse = "\n"),sep="")
            if(add_max_qvalue){
                added_numbers <- paste(added_numbers,"\n","Max q-value:",sprintf("%.2e",max(x$q_value)),sep="")
            }
            if(is.null(text_position)){
                gg1 <- gg1 + annotate("text", x = 0.01*1.05, y = 0.9*max_fdp, label = added_numbers, color = "black", size = 3,hjust = 0)
            }else{
                gg1 <- gg1 + annotate("text", x = text_position[1], y = text_position[2], label = added_numbers, color = "black", size = 3,hjust = 0)
            }
        }
    }
   
    # library(ggpubr)
    # options(repr.plot.width = 12, repr.plot.height = 6)
    # gg <- ggarrange(gg1,gg2,ncol = 2,  common.legend = TRUE)
    # options(jupyter.plot_mimetypes = "image/png")
    if(return_data){
        return(list(gg=gg1,data=dat,added_numbers=added_numbers))
    }else{
        return(gg1)
    }

    #pdf("fdp_fdr_lib_qvalue_protein.pdf",width = 8,height = 4)
    #print(gg)
    #dev.off()
}



plot_fdp_fdr_multiple=function(dat,added_numbers=NULL,n_row_plots=2,n_col_plots=2,fdr_max=NULL,fig_title=NULL,scale_xy=TRUE,add_numbers=FALSE,r=1,fixed_fdr_max=FALSE,max_x=NA,max_y=NA,
                         color_mapping=NULL,
                         legend_position=NULL,
                         fdr_decimal_place=2,
                         return_data=FALSE) {

    max_fdp <- max(c(dat$FDP,dat$q_value))
    if(!is.null(fdr_max)){
        if(fixed_fdr_max){
            max_fdp <- fdr_max
        }else{
            max_fdp <- min(c(fdr_max,max_fdp))
        }
    }
    
    gg1 <- ggplot(dat,aes(x=q_value,y=FDP,color=Method)) + 
            geom_abline(slope = 1,intercept = 0,color="gray")+
            #rasterise(geom_line(), dpi = 300) + 
            geom_line()+
            facet_wrap(.~tool,nrow=n_row_plots,ncol=n_col_plots)+
            xlab("FDR threshold")+
            ylab("Estimated FDP")+
            theme_bw()+
            #geom_segment(x=0.01,xend=0.01,y=0,yend=0.01,linewidth=0.3,color="blue",linetype=2)+
            theme_pubr(base_size = 12,border = TRUE)

    if(!is.null(color_mapping)){
        gg1 <- gg1 + scale_color_manual(values = color_mapping)
    }

    if(scale_xy){

        if(!is.na(max_x) || !is.na(max_y)){
            gg1 <- gg1 + geom_vline(xintercept = 0.01,linetype=2,color="blue")
            if(!is.na(max_x)){
                gg1 <- gg1 + xlim(0,max_x) + scale_x_continuous(labels = scales::percent,limits =c(0,max_x))
            }else{
                gg1 <- gg1 + xlim(0,max_fdp)+ scale_x_continuous(labels = scales::percent,limits =c(0,max_fdp))
            }
            if(!is.na(max_y)){
                gg1 <- gg1 + ylim(0,max_y) + scale_y_continuous(labels = scales::percent,limits =c(0,max_y))
            }else{
                gg1 <- gg1 + ylim(0,max_fdp) + scale_y_continuous(labels = scales::percent,limits =c(0,max_fdp))
            }
        }else{
            gg1 <- gg1 + geom_vline(xintercept = 0.01,linetype=2,color="blue")+
                xlim(0,max_fdp)+
                ylim(0,max_fdp)+
                scale_y_continuous(labels = scales::percent,limits =c(0,max_fdp))+
                #scale_y_pct()+
                #scale_x_pct()+
                scale_x_continuous(labels = scales::percent,limits =c(0,max_fdp))
        }
    }
    #theme(legend.position = "top",plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"))+    
    ## legend on the botton right
    ## no background color for the legend
    if(is.null(legend_position)){
        gg1 <- gg1 + theme(legend.position = c(0.65, 0.16),plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"),legend.background = element_blank())    
    }else{
        gg1 <- gg1 + theme(legend.position = legend_position,plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"),legend.background = element_blank())    
    }
    
            
    if(!is.null(fig_title)){
        gg1 <- gg1 + ggtitle(fig_title)
    }

    added_numbers <- NULL
    if(add_numbers){
        
    }
   
    return(gg1)
}
