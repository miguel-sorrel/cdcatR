.onAttach <- function(libname, pkgname){
  packageStartupMessage(StartWelcomeMessage())
}

StartWelcomeMessage <- function(){
  paste(c("==============================\n",
          "cdcatR Package",
        " [Version ", utils::packageDescription("cdcatR")$Version,
        "; ",utils::packageDescription("cdcatR")$Date, "]\n",
        "More information: https://github.com/miguel-sorrel/cdcatR\n",
        "==============================\n"),
        sep="")
}
