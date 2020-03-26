


#----------------------------------------------------------------
# Setup working directories and output folder labels 
#-----------------------------------------------------------------
# Set Working directory file, where assessments are stored 
File = "C:/Work/Research/GitHub/JABBAbeta/JABBAv2.0development"
# Set working directory for JABBA R source code
JABBA.file = File
# JABBA version
# NA
# Set Assessment file: assement folder within File that includes .csv input files
assessment = "bet.iccat"
# add specifier for assessment (File names of outputs)
setwd(file.path(File))

source("jabba.functions.R")

load(file.path(getwd(),assessment,"bet.iccat.rdata"),verbose=T)


# Load jabba base pkgs
jabba_libs()
# Build JABBA model
jbinput = build_jabba(catch=bet$catch,cpue=bet$cpue,se=bet$se,assessment=assessment,model.type = "Fox")
# fit JABBA model
fit = fit_jabba(jbinput)
# Run all plots
jbplots(fit)
