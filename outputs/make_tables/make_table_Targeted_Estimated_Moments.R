# ===============================================================
# Prepare tex table with target and estimated parameters
# ======================================================
# Define format of figures:
nb.dec <- 2 # number of decimal numbers
format.nb <- paste("%.",nb.dec,"f",sep="")


latex.table <- rbind(
  paste("$\\mathbb{E}_0$($T_{AT,2040}$)&",
        make.entry(target_vector[15],format.nb),"\\degree C&",
        make.entry(mim[15],format.nb),"\\degree C&",
        "RCP4.5+RCP6.0",
        "\\\\",sep=""),
  paste("$\\mathbb{E}_0$($T_{AT,2100}$)&",
        make.entry(target_vector[1],format.nb),"\\degree C&",
        make.entry(mim[1],format.nb),"\\degree C&",
        "RCP4.5+RCP6.0",
        "\\\\",sep=""),
  paste("$\\mathbb{V}\\mathrm{ar}_0$($T_{AT,2100}$)&",
        make.entry(target_vector[2],format.nb),"$^2$\\degree C&",
        make.entry(mim[2],format.nb),"$^2$\\degree C&",
        "RCP4.5+RCP6.0",
        "\\\\",sep=""),
  # paste("Expected contribution of FL to $T_{AT,2100}$&",
  #       make.entry(target_vector[3],format.nb),"\\degree C&",
  #       make.entry(mim[3],format.nb),"\\degree C&",
  #       "\\citet*{Burke_Hartley_Jones_2012}",
  #       "\\\\",sep=""),
  paste("$\\mathbb{E}_0$($Cum_{\\mathcal{E},2100,FL}$|$\\mu=0$) &",
        make.entry(target_vector[4],paste("%.",nb.dec=0,"f",sep="")),"\\,GtCO$_2$&", 
        make.entry(mim[4],paste("%.",nb.dec=0,"f",sep="")),"\\,GtCO$_2$&",
        "\\citet{Burke_etal_2013}", #Burke_Hartley_Jones_2012/ Schaefer et al, 2014
        "\\\\",sep=""),
  paste("$\\mathbb{C}\\mathrm{ov}_0$($T_{AT,2100},Cum_{D,2100}$)/                 
        $\\mathbb{V}\\mathrm{ar}_0$($T_{AT,2100}$)&",
        make.entry(target_vector[5],format.nb),"&",
        make.entry(mim[5],format.nb),"&",
        "\\citet{Burke_Hsiang_Miguel_2015}",
        "\\\\",sep=""),                                                         #linear regression slope
  paste("Long-term rate in 2100 &",
        make.entry(target_vector[6],format.nb),"\\%","&",
        make.entry(mim[6],format.nb),"\\%","&",
        "US Treasury",
        "\\\\",sep=""),
  paste("$\\mathbb{E}_0(H_{2100})$ &",
        make.entry(target_vector[7],format.nb),"m&",
        make.entry(mim[7],format.nb),"m&",
        "RCP4.5+RCP6.0",
        "\\\\",sep=""),
  paste("$\\mathbb{V}\\mathrm{ar}_0$($H_{2100}$) &",
        make.entry(target_vector[8],format.nb),"$^2$m","&",
        make.entry(mim[8],format.nb),"$^2$m","&",
        "\\citet{Mengel_et_al_2016}",
        "\\\\",sep="")
)


latex.file <- paste("outputs/Tables/table_moments.txt",sep="")

write(latex.table, file = latex.file)

#* Explain other targets in the note.
