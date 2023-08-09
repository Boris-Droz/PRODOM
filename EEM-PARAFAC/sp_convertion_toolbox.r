###### R-Toolbox to open spectroflourometer (sp files) data obtained with FLWInlab software ######

# authors: Stefan Seeger & Elena Fernandez Pascual
# Created: 2015-11-18
# Last revision: 2023-07-17

#----------------------------------------------------

# contained functions and their purpose:

# loadSpBinary()				   load data from one binary sp-file
# loadSpASCII()					   load data from one ASCII sp-file
# loadSp()						     wrapper function for loadSpBinary() and loadSpASCII()
# loadSpList()					   load data from a list of sp-files and merge them to a excitation matrix

#----------------------------------------------------

# function to load data from binary FLWinlab sp-files
loadSpBinary <- function(filepath){

	getPos <- function(filepath, searchValue, minSearchPos=1, maxSearchPos=2000){
		for(pos in minSearchPos:maxSearchPos){
			con <- file(filepath, 'rb')
			buff <- readBin(con, "raw", n=pos)
			x <- readBin(con, "raw", n=length(searchValue))
			close(con)
			if(all(x == searchValue)) return(pos+length(searchValue))
		}
		return(NA)
	}
	pos1 <- getPos(filepath,c(0x2b, 0x8c, 0x0e, 0x00, 0x00, 0x00, 0x1c, 0x75))
	pos2 <- getPos(filepath,c(0x72, 0x8b, 0x12, 0x00, 0x00, 0x00, 0x1d, 0x75),pos1)
	filesize <- file.info(filepath)$size
	
	metaInfo<-list()
	con<-file(filepath, 'rb')
	buff <- readBin(con, "raw", n=pos1)
	metaInfo$exLength <- readBin(con, "numeric", n=1)
	buff <- readBin(con, "raw", n= pos2-pos1-8)
	metaInfo$emStart = readBin(con, "numeric", n=1)
	metaInfo$emEnd = readBin(con, "numeric", n=1)
	buff <- readBin(con, "raw", n=8)
	metaInfo$valMin <- readBin(con, "numeric", n=1)
	metaInfo$valMax <- readBin(con, "numeric", n=1)
	buff <- readBin(con, "raw", n=8)
	metaInfo$emStepWidth <- readBin(con, "numeric", n=1)
	metaInfo$emStepCount <- (metaInfo$emEnd - metaInfo$emStart)/metaInfo$emStepWidth +1
	dataStartPos = filesize - metaInfo$emStepCount*8
	buff <- readBin(con, "raw", n=dataStartPos - (pos2+56))
	rawData <- readBin(con, "numeric", n = metaInfo$emStepCount*8)
	close(con)
	
	emissionWavelengths <- seq(metaInfo$emStart, metaInfo$emEnd, by= metaInfo$emStepWidth)
	if(length(emissionWavelengths) != length(rawData)) warning(sprintf("inconsistent data read from %s!", filepath))
	return(list(excitationLevels = rawData, emissionWavelengths = emissionWavelengths, metaInfo = metaInfo))
}

# function to load data from ASCII FLWinlab sp-files
loadSpASCII <- function(filepath){
	metaInfoRaw <- readLines(filepath, n=54)
	metaInfo = list(emEnd = as.numeric(metaInfoRaw[10]),
					emStart = as.numeric(metaInfoRaw[48]),
					emStepWidth = as.numeric(metaInfoRaw[49]),
					emStepCount = as.numeric(metaInfoRaw[50]),
					exLength = as.numeric(metaInfoRaw[22]),
					valMax = as.numeric(metaInfoRaw[52]),
					valMin = as.numeric(metaInfoRaw[53]))
	scanData <- read.table(filepath, skip=54, sep='\t')
	
	if(nrow(scanData) != metaInfo$emStepCount) warning(sprintf("inconsistent data read from %s!", filepath))
	return(list(excitationLevels = scanData[,2], emissionWavelengths = scanData[,1], metaInfo = metaInfo))
}

# metaFunction to load data from binary and ASCII FLWinlab sp-files
loadSp <-function(filepath){
	firstLine <- readLines(filepath, n=1, warn=FALSE)
	if(firstLine == "PEPE2D constant interval DataSet file") return(loadSpBinary(filepath))
	return(loadSpASCII(filepath))
}

# function to load a set of FLWinlab sp files and merge them into a matrix
loadSpList <- function(filelist){
  for(i in 1:length(filelist)){
    cat(filelist[i],'\n')
    dat <- loadSp(filelist[i])
    if(i==1){
      excitationLengths <- dat$metaInfo$exLength
      EEM <- dat$excitationLevels
      refMetaInfo <- dat$metaInfo[c("emStart","emEnd","emStepWidth")]
    }
    else{
      excitationLengths <- c(excitationLengths,dat$metaInfo$exLength)
      EEM <- cbind(EEM, dat$excitationLevels)
    }
  }
  return(list(dataMatrix = EEM, excitationLengths = excitationLengths, emissionLengths = dat$emissionWavelengths))
}


# function to load all FLWinlab sp files and combine them into one matrix
loadSpDir <- function(spDir){
	filenameList <- dir(spDir)
	filenameList <- filenameList[grep("\\.sp$", filenameList)]
	return(loadSpList(sprintf("%s/%s", spDir, filenameList)))
}

