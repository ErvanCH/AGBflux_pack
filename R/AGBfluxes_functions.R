#' CTFS-formated data preparation
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Main routine to format, detect major obvious errors, and gap-fill those errors in CTFS-formated data
#' @param site the full name of your site (in lower case); e.g., 'barro colorado island'
#' @param stem TRUE or FALSE, using the stem data (stem=TRUE) rather than the tree data (i.e. called 'full', stem=FALSE)
#' @param taper.correction TRUE or FALSE, apply Cushman et al (2014) taper correction
#' @param fill.missing TRUE or FALSE, interpolate missing DBH values
#' @param use.palm.allometry TRUE or FALSE, if TRUE, compute biomass of palm trees using a palm-specific allometric model from Goodman et al. (2013)
#' @param flag.stranglers TRUE or FALSE, individuals of known strangler fig species greater than 'dbh.stranglers' are flagged 
#' @param dbh.stranglers (optional) minimal diameter (in mm) of strangler figs, default = 500 
#' @param maxrel a numeric value: the threshold for flagging major errors in productivity, applied as absval(individual tree productivity)>maxrel*(average productivity per hectare)
#' @param output.errors TRUE or FALSE, output all records for trees with major errors in productivity to a csv file
#' @param DATA_path the pathname where the data are located 
#' @param exclude.interval NULL by default. If needed a vector (e.g. c(1,2)) indicating which census interval(s) must be discarded from computation due, for instance, to a change in measurement protocol
#' @return a data.table (data.frame) with all relevant variables.
#' @import data.table  # load data.table package
#' @export

data.prep <- function(site,stem,taper.correction,fill.missing,use.palm.allometry,flag.strangler,dbh.stranglers,maxrel,graph.problem.trees,output.errors,DATA_path,exclude.interval=NULL) {
	site <- tolower(site)
	INDEX <- match(tolower(site),site.info$site)
	if (is.na(INDEX)) {			stop("Site name should be one of the following: \n",paste(levels(factor(site.info$site)),collapse=" - ")) }
	
	if(!exists("DATA_path")){
		DATA_path <- paste0(path_folder,"/data/")
	}
	
	files <-list.files(DATA_path)
	ifelse(stem,files <- files[grep("stem",files)], files <- files [grep("full",files)])
	ifelse(!dir.exists(file.path(paste0(path_folder,"/output"))), dir.create(file.path(paste0(path_folder,"/output"))), FALSE)
	
	
	# Create the receiving data.frame
	df <- data.frame("treeID"=NA, "stemID"=NA, "tag"=NA, "StemTag"=NA, "sp"=NA, "quadrat"=NA, "gx"=NA, "gy"=NA,"dbh"=NA,"hom"=NA, "ExactDate"=NA, "DFstatus"=NA, "codes"=NA, "date"=NA,"status"=NA,"agb"=NA,"CensusID"=NA,"year"=NA)
	
	for (i in 1:length(files)) {
		temp <- setDT(LOAD(paste(DATA_path,files[i],sep="/")))
		temp$CensusID <- i
		temp[,year:= round(mean(as.numeric(format(as.Date(date, origin='1960-1-1'),'%Y')),na.rm=T))]
		temp  <- temp[,match(names(df),names(temp)),with=FALSE]
		df <- rbind(df,temp)
	}
	rm(temp)
	table(df$year,df$CensusID)
	df <- df[-1,]
	
	df <- data.correction(df,taper.correction,fill.missing,stem)
	print("Step 1: data correction done.")
	
	df <- computeAGB(df,site,palm,DATA_path=NULL)
	print("Step 2: AGB calculation done.")
	
	# save(df,file=paste0(path_folder,"/output/",site,"_rawcensus_corrected.Nov17.Rdata"))
	# stem=T;taper.correction=T;fill.missing=T;palm=T;flag.stranglers=T;maxrel=0.2;graph.problem.trees=F;output.errors=F;exclude.interval=1;assign("dbh.stranglers",500)
	# df <- LOAD("C:/Users/Ervan/Documents/Work/R/CTFS/Code_compil/biomass_fluxes_BCI/output/barro colorado island_rawcensus_corrected.Nov17.Rdata")
	# 
	DF <- format.interval(df,flag.stranglers,dbh.stranglers)
	print("Step 3: data formating done.")
	
	DF <- flag.errors(DF,site,flag.stranglers=flag.stranglers,maxrel=maxrel,graph.problem.trees=graph.problem.trees,output.errors=output.errors,exclude.interval=exclude.interval)
	print("Step 4: errors flagged. Saving corrected data into output folder.")
	
	save(DF,file=paste0(path_folder,"/output/",site,"_census_intervals_20perc.Rdata"))
	# save(DF,file=paste0(path_folder,"/output/",site,"_census_intervals_20perc_nobrokenstems.Rdata"))
	rm(list=setdiff(ls(), c("DF","path_folder","SITE", lsf.str())))
	return(DF)
}

#' Trim and correct data
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Stack all censuses together and correct DBH, if required
#' @param taper.correction TRUE or FALSE, are you willing to apply Cushman et al (2014) taper correction?
#' @param fill.missing TRUE or FALSE, are you willing to extrapolate missing DBH from surrounding DBH?
#' @param stem is the function applied over stem (stem=TRUE) or tree (stem=FALSE) data?
#' @import data.table  # load data.table package
#' @return a data.table (data.frame) with all relevant variables.
#' @export
data.correction <- function(df,taper.correction,fill.missing,stem) {
	if (stem) {
		df[,"id" :=paste(treeID,stemID,sep="-")] # creat a unique tree-stem ID
		df <- df[order(id,CensusID)]
	} else { df[,"id" := treeID]}
	
	df[,status1:=normal.stat(.SD),by=id] # check that the status is consistent over all censuses (e.g. can't be dead and alive at next census)
	df <- df[!status1%in%c("P","Dr")] # discard all priors & replicated dead trees
	
	# Add average date of census when missing for alive trees (mandatory for interpolating DBHs)
	df[,"year":=round(mean(as.numeric(format(as.Date(ExactDate, origin='1960-1-1'),'%Y')),na.rm=T)),by=CensusID] # Assign 1 year per census
	DATE <- df[,.(date=mean(date,na.rm=T)),by=year]
	df <- within(df,date[is.na(date)] <- DATE$date[match(df[is.na(date),"year"]$year,DATE$year)])
	
	# remove stems without any measurement in any census
	NO.MEASURE <- df[,all(is.na(dbh)),by=id]
	df <- df[!id%in%NO.MEASURE$treeID[NO.MEASURE$V1]] 
	
	# Taper correction or missing values: -> fill gaps for missing values
	system.time(df[, c("dbh2","hom2") := corDBH(.SD,taper.correction=taper.correction,fill.missing=fill.missing), by=id]) # might be time consuming (25 minutes for BCI)
	table(is.na(df$dbh2),df$status1)
	NO.MEASURE <- df[,all(is.na(dbh2)),by=treeID]
	table(NO.MEASURE$V1)
	df <- df[!treeID%in%NO.MEASURE$treeID[NO.MEASURE$V1]] # remove trees without any measurement
	return(df)
}

#' Data correction
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Perform two mains tasks: (a) apply a taper correction when POM is > 130 cm, and (b) linear interpolation values when missing DBHs. Interpolation of missing values is done by averaging surrounding available DBH values.
#' @param DF a data.table
#' @param taper.correction TRUE or FALSE, are you willing to apply Cushman et al (2014) taper correction?
#' @param fill.missing TRUE or FALSE, are you willing to extrapolate missing DBH from surrounding DBH?
#' @import data.table  # load data.table package
#' @return a data.table (data.frame) with all relevant variables.
#' @export

# Debug
# for (i in levels(factor(df$id))) {
# 	corDBH(df[id==i],taper.correction==T,fill.missing==T)
# }
#   DT <- df[id==i]

corDBH <- function(DT,taper.correction,fill.missing) {
	dbh2 <- DT[status1=="A",'dbh'][[1]]
	hom2 <- DT[status1=="A",round(hom*100)/100]   # round hom to be at 1.3 meter (avoiding odd rounding)
	loc <- which(is.na(dbh2))
	DATE <- DT[status1=="A","date"][[1]]
	hom2[is.na(hom2) & !is.na(dbh2)] <- 1.3
	
	# 1. Apply Cushman's correction to trees with POM changed
	if (taper.correction & any(hom2 > 1.3,na.rm=T)) {
		hom2[is.na(hom2)] <- 1.3
		if (taper.correction & any(hom2 > 1.3,na.rm=T)) {
			loc1 <- hom2!=1.3 & !is.na(dbh2)
			if(!any(is.na(loc1))) {
				dbh2[loc1] <- round(dbh2[loc1]*exp(0.0247*(hom2[loc1]-1.3)),1)
				hom2[loc1] <- 1.3
			}}
	}
	
	# 2. Fill gaps
	if (fill.missing & any(is.na(dbh2))) {
		
		# dead stem
		if (any(grep("D",DT$status1))) {
			ifelse(length(dbh2)%in%loc,RULE<-2,RULE<-1)  # if last value is NA, last valid DBH measure is replicated (RULE=2)
			tryCatch(approx(DATE,dbh2,rule=RULE,xout=DT$date[loc])$y,error=function(e) print(unique(paste("Stem id",DT$id,"generates a problem."))))
			if(any(is.na(approx(DATE,dbh2,rule=RULE,xout=DATE[loc])$y))) {	print(DT$id)}
			dbh2[loc] <- approx(DATE,dbh2,rule=RULE,xout=DATE[loc])$y
			hom2[loc] <- approx(DATE,hom2,rule=RULE,xout=DATE[loc])$y
			
			dbh2 <- c(dbh2,NA)  # add NA at year of death
			hom2 <- c(hom2,NA)
		} # end of dead stems
		else {
			# alive stems
			ifelse(length(dbh2)%in%loc,RULE<-2,RULE<-1)  # if last value is NA, last valid BDBH measure is replicated (RULE=2)
			tryCatch(approx(DT$date,dbh2,rule=RULE,xout=DT$date[loc])$y,error=function(e) print(paste("Stem id",DT$id,"generates a problem.")))
			dbh2[loc] <- approx(DT$date,dbh2,rule=RULE,xout=DT$date[loc])$y
			hom2[loc] <- approx(DT$date,hom2,rule=RULE,xout=DT$date[loc])$y
		} # end of alive stems
		
		# Resprouts
		RESP <- grep("\\bR\\b",DT$codes[is.na(DT$dbh)])
		if(length(RESP)>0) {
			dbh2[RESP] <- NA
			hom2[RESP] <- 1.30
		}
		
	} else {# end of fill.missing
		
		# Replicate dbh & hom if no NA values or fill.missing = FALSE
		dbh2 <- DT[,'dbh'][[1]]
		hom2 <- DT[,round(hom*100)/100]
		hom2[is.na(hom2)] <- 1.3
		if (taper.correction & any(hom2 > 1.3,na.rm=T)) {
			loc1 <- hom2!=1.3 & !is.na(dbh2)
			if(!any(is.na(loc1))) {
				dbh2[loc1] <- round(dbh2[loc1]*exp(0.0247*(hom2[loc1]-1.3)),1)
				hom2[loc1] <- 1.3
			}}
	} # end if fill.missing = F or no missing dbh
	return(list(dbh2,hom2))
}

#' Biomass computation
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Allocate wood density and compute above-ground biomass using the updated model of Chave et al. (2014), given in Rejou-Mechain et al. (2017). Palm trees (palm=T) are computed using a different allometric model (Goodman et al. 2013).
#' @param site provide the full name of your site (in lower case) i.e. 'barro colorado island'
#' @param palm TRUE or FALSE, if TRUE, biomass of palm trees is computed through a specific allometric model (Goodman et al. 2013)
#' @param DATA_path allows to provide a different path where the data are located
#' @import data.table  # load data.table package
#' @return a data.table (data.frame) with all relevant variables.
#' @export

computeAGB <- function(df,site,use.palm.allometry=T) {
	if(!exists("DATA_path")){
		DATA_path <<- paste0(path_folder,"/data/")
	}
	## Allocate wood density
	df <- density.ind(df,site,wsg=WSG)
	
	# Compute biomass
	df$agb <- AGB.comp(site,df$dbh2, df$wsg,H = NULL)
	
	# Compute biomass for palms
	if (use.palm.allometry) {
		SP <-  LOAD(paste0(DATA_path,list.files(DATA_path)[grep("spptable",list.files(DATA_path))]))
		if(is.na(match("genus",tolower(names(SP))))) {
			trim <- function (x) gsub("^\\s+|\\s+$", "", x)
			SP$genus <-  trim(substr(SP$Latin,1,regexpr(" ",SP$Latin)))
			SP$species <-  trim(substr(SP$Latin,regexpr(" ",SP$Latin),50))
			SP <- SP[,c("sp","genus","species","Family")]
			SP$name <- paste(SP$genus,SP$species,sep=" ")
			names(SP) <- c("sp","genus","species","Family","name")
			df <- merge(df,SP,by="sp",all.x=T)
			agbPalm <- function(D) { exp(-3.3488 + 2.7483*log(D/10) + ((0.588)^2)/2)/1000 }
			df['Family'=="Arecaceae","agb":= agbPalm(dbh2)]
		} else {
			SP <- SP[,c("sp","Genus","Species","Family")]
			SP$name <- paste(SP$Genus,SP$Species,sep=" ")
			names(SP) <- c("sp","genus","species","Family","name")
			df <- merge(df,SP,by="sp",all.x=T)
			agbPalm <- function(D) { exp(-3.3488 + 2.7483*log(D/10) + ((0.588)^2)/2)/1000 }
			df['Family'=="Arecaceae","agb":= agbPalm(dbh2)]
		}
	}
	return(df)
}

#' Format census intervals
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Create census intervals (i.e. put consecutive census side by side), assign status by tree (i.e. alive (A),dead (D), recruited (R) or resprout (Rsp))
#' @param df a data.table
#' @param flag.stranglers TRUE or FALSE, individuals of known strangler fig species greater than 'dbh.stranglers' are flagged 
#' @param dbh.stranglers (optional) minimal diameter (in mm) of strangler figs, default = 500 
#' @import data.table  # load data.table package
#' @return a formated data.table.
#' @export

format.interval <- function(df,flag.stranglers,dbh.stranglers,code.broken=NULL) {
	YEAR <- unique(df$year)
	
	# Receiveing data set
	DF2 <- data.table("treeID"=NA,"dbh1"=NA,"dbhc1"=NA,"status1"=NA,"code1"=NA,"hom1"=NA,"sp"=NA,"agb1"=NA,"date1"=NA,"dbh2"=NA,"dbhc2"=NA,"status2"=NA,"code2"=NA,"hom2"=NA,"agb2"=NA,"agbl"=NA,"date2"=NA,"agb1.surv"=NA,"interval"=NA,"year"=NA)
	
	for (j in 1:(length(YEAR)-1)) {  # 4 minutes to run
		i1 <- df[year==YEAR[j] & status1 != "D", .I[which.max(dbh2)], by = treeID] # keep only information for the biggest alive stem per treeID
		A1 <- df[i1$V1,c("treeID","dbh","dbh2","status1","codes","hom","agb","id","sp")]
		names(A1) <- c("treeID","dbh1","dbhc1","status1","code1","hom1","agb","id1","sp")
		B1 <- df[year==YEAR[j] & status1 != "D",list("agb1"=sum(agb,na.rm=T),"date1"=mean(date,na.rm=T)),by=treeID]
		BB <- merge(B1,A1,by="treeID",all.x=T)
		cens1 <- BB[,c("treeID","dbh1","dbhc1","status1","code1","hom1","agb","id1","sp","agb1","date1")]
		
		i2 <- df[year==YEAR[j+1] & status1 != "D" , .I[which.max(dbh2)], by = treeID]
		A2 <- df[i2$V1,c("treeID","dbh","dbh2","hom","status1","codes","id")]
		names(A2) <- c("treeID","dbh2","dbhc2","hom2","status2","code2","id2")
		B2 <- df[year==YEAR[j+1],list("agb2"=sum(agb[status1!="D"],na.rm=T),"date2"=mean(date,na.rm=T)),by=treeID]
		BB2 <- merge(B2,A2,by="treeID",all.x=T)
		cens2 <- BB2[,c("treeID","dbh2","dbhc2","status2","code2","hom2","agb2","date2","id2")]
		cens2 <- within(cens2,status2[is.na(status2)] <- "D")
		
		ID <- data.table(treeID=unique(c(cens1$treeID,cens2$treeID)))
		ID <- merge(ID,cens1,by='treeID',all.x=T)
		ID <- merge(ID,cens2,by='treeID',all.x=T)
		
		ID[,"broken":=ifelse(id1!=id2 & dbhc2<dbhc1-4*(0.0038*dbhc1 + 0.927) & dbhc1>100 & sp!="oenoma" & hom2<=hom1,1,0)] # flag large broken main stems to be added to losses
		ID[broken==1,"agbl":=agb]
		ID[broken==1,"agb1.surv":=agb1-agbl]
		# ID$broken <- 0  ## when not flaggigng broken stems
		# ID$agbl <- as.numeric(NA)  ## when not flaggigng broken stems
		# ID$agb1.surv <- as.numeric(NA)  ## when not flaggigng broken stems
		if(exists("code.broken")) {
			idx <- ID[broken==1,.I[!grepl(paste0("\\b",code.broken,"\\b"),code2) ]]	
		} else { 
			idx <- ID[broken==1,.I[!grepl("\\bR\\b",code2) ]]
		}
		ID[broken==1,] <- within(ID[broken==1,],code2[idx] <- "broken")
		ID[,c("agb","broken","id1","id2"):=NULL]
		ID$interval=j
		ID$year=YEAR[j+1]
		DF2 <- rbind(DF2,ID)
	}
	DF2 <- DF2[-1,]
	
	# Add coordinates & other mandatory information
	COORD <- df[,.(gx=round(10*gx[!is.na(gx)])/10,gy=round(10*gy[!is.na(gx)])/10,quadrat=quadrat,name=name),by=treeID]
	COORD <- unique(COORD[!duplicated(COORD$treeID)])
	DF2 <- merge(DF2,COORD,by="treeID",all.x=T)
	
	# Add average date of census when missing
	DATE <- df[,.(date=mean(date,na.rm=T)),by=year]
	DATE[,"year2" := c(year[2:length(year)],NA)]
	DF2 <- within(DF2,date1[is.na(date1)] <- DATE$date[match(DF2[is.na(date1),"year"]$year,DATE$year2)])
	DF2 <- within(DF2,date2[is.na(date2)] <- DATE$date[match(DF2[is.na(date2),"year"]$year,DATE$year)])
	DF2$int <- (DF2$date2 - DF2$date1)/365.5  # census interval in days
	# 
	# Update status for recruited trees
	DF2[, nrow := seq_len(.N), by = treeID]
	DF2 <- within(DF2,status1[is.na(status1) & !is.na(dbh2) & nrow==1] <- "P")  # recruited trees or Z code (=?)
	DF2[,nrow:= NULL]
	
	# Assign status
	DF2[, c("code","dHOM") := assign.status(.SD), by=treeID]
	
	# Remove unnecessary rows
	DF2[code=="D", nrow := seq_len(.N), by = treeID]
	DF2 <- within(DF2,nrow[is.na(nrow)] <- 0)
	DF2 <- 	DF2[nrow<2,]  # keeps only 1 line for dead trees
	DF2[,"nrow":=NULL]
	
	# Compute annualized fluxes
	DF2[code%in%c("A","AC","B"),prod.g := (agb2-agb1)/int,by=treeID] # annual prod for alive trees & multi-stems
	DF2[code=="B",prod.g := (agb2-agb1.surv)/int,by=treeID] # annual prod for alive trees & multi-stems
	DF2[code%in%c("R","Rsp"),prod.r:=agb2/int,by=treeID] # annual prod for resprouts and recruits
	DF2[code=="D",agbl:= agb1,by=treeID]
	DF2[,loss:=agbl/int,by=treeID] # annualized loss for dead trees
	
	# Flag large strangler figs
	if(flag.stranglers) {
		DF2$ficus <- 0
		ficus$name <- paste(ficus$Genus,ficus$Species,sep=" ")
		FIC <- match(DF2$name,ficus$name[ficus$Strangler=="Yes"])
		if(!exists("dbh.stranglers")) {assign("dbh.stranglers",500) }
		DF2 <- within(DF2,ficus[!is.na(FIC) & dbhc1>dbh.stranglers]<-1)
	}
	return(DF2)
}


#' Flag major errors
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Identify trees with major errors in DBH measurments. A major error correspond to a relative individal productivity (or growth) is above a given percentage (set by 'maxrel') of the mean productivity computed at a site. Additionnaly, flagged trees that died at next census interval are also flagged. Option to see DBH measurement (=draw.graph) of flagged trees or print a csv (output.errors) are given.
#' @param DF a data.table
#' @param site provide the full name of your site (in lower case) i.e. 'barro colorado island'
#' @param strangler TRUE or FALSE, if TRUE, strangler figs tree are flagged (upon a list to published soon)
#' @param maxrel a numeric value setting the threshold over which relative productivity is assumed to be too high (usually set at 20 percents)
#' @param output.errors TRUE or FALSE, output all records for trees with major errors in productivity to a csv file
#' @param exclude.interval a vector (i.e. c(1,2)) indicating if a set of census intervals must be discarded from computation due for instance to a change in  protocol of measurment
#' @import data.table  # load data.table package
#' @return a data.table (data.frame) with all relevant variables.
#' @export
#'
flag.errors <- function(DF,site,flag.stranglers,maxrel,graph.problem.trees,output.errors,exclude.interval) {
	mean.prod <- determine.mean.prod(DF,site,flag.stranglers,exclude.interval)
	DF[,prod.rel:=as.numeric(NA),]
	DF[,prod.rel:=prod.g/mean.prod] # relative contribution to average total productivity
	DF[,error:=0]
	DF <- within(DF,error[prod.rel>maxrel & dHOM==0 & code!="D"] <- 1)
	DF <- within(DF,error[prod.rel<(-maxrel) & dHOM==0 & code!="D"] <- -1)
	
	# Flag census after a major error
	POSI <- DF[,.I[error!=0]+1,]
	POSI2 <- DF[POSI,.I[error==0 & dHOM==0 & code=="D"],]
	DF[,error.loss:=0]
	DF <- within(DF,error.loss[POSI[POSI2]] <- 1) # flag consecutive census
	if (flag.stranglers) {
		DF <- within(DF,error.loss[ficus==1 & code=="D"] <- 1) # flag dead strangler figs	
	}
	ID <- DF[error!=0 & !code%in%c("D","R"),nrow(.SD)>1,by=treeID]
	ID <- ID$treeID[ID$V1]
	if (length(ID)/12 > 10) {
		A <- 	menu(c("Y", "N"), title="There are more than 144 trees (10 pages) to be printed. Do you want to print them?") 
		ifelse(A==1,graph.problem.trees<-T,graph.problem.trees<-F)
	} 
	if(graph.problem.trees) { # Plot trees with large major error
		YEAR <- levels(factor(DF$year))
		CX=2
		a=0
		i = 0
		GRAPH = list()
		
		for (n in 1:length(ID)){
			i = i + 1
			DF[,year:=as.numeric(year)]
			DF$dbh1 <- round(DF$dbh1,1)
			X <- DF[treeID==ID[n] & !code%in%c("R")][order(year)]
			X$point <- 0
			X$point[X$error!=0] <- 1
			Y <- DF[treeID==ID[n] & !code%in%c("D","R")][order(year)]
			YY <- Y[,.(year=max(year),name=unique(name),d2=round(dbhc2[year==max(year)],1),d02=round(dbh2[year==max(year)],1),hom2=round(hom2[year==max(year)],2)),by=treeID]
			Y$line <- 0
			Y$line[Y$dHOM==0] <- 1
			Y$point <- 0
			Y$point[Y$error!=0] <- 1
			
			GRAPH[[i]] = ggplot(X,aes(x=year,y=dbhc1)) + geom_point(size=2) + geom_segment(data=Y,aes(x=year,y=dbhc1,xend=year+5,yend=dbhc2,linetype=as.factor(line))) + geom_point(data=X[point==1],aes(x=year+5,y=dbhc2),col=2)	+ labs(title=paste0(unique(X$name)," (",ID[n],")"), x=" ",y="dbh (mm)") + geom_text(data=Y,aes(x=year,y=dbh1-(0.05*dbh1)),label=round(Y$hom1,2),cex=CX)	+ geom_text(data=YY,aes(x=year+5,y=d02-(0.05*d02)),label=round(YY$hom2,2),cex=CX) + geom_text(data=Y,aes(x=year,y=0.3*max(dbhc2)),label=Y$dbh1,cex=CX,angle=90,vjust=1) + geom_text(data=YY,aes(x=year+5,y=0.3*max(d2)),label=YY$d02,cex=CX,angle=90,vjust=1) + theme(plot.title = element_text(size=5*CX,face="bold"),axis.title.y= element_text(size=5*CX,,face="bold"),axis.text.y = element_text(size=4*CX),axis.text.x = element_text(size=4*CX,vjust=0,angle=30),panel.background = element_blank(),strip.text = element_text(size = 4*CX,face="bold"),strip.background = element_rect("lightgrey"),panel.spacing = unit(0.1, "lines")) + scale_linetype_manual(values=c('0'="dashed",'1'="solid")) + guides(linetype=F,colour=F) + scale_x_continuous(limits=c(min(as.numeric(YEAR))-3,max(as.numeric(YEAR))+3),breaks = as.numeric(YEAR)) + scale_y_continuous(limits=c(0.2*max(YY$d2),max(X$dbhc2,X$dbhc1)))
			
			
			if (i %% 15 == 0) { ## print 8 plots on a page
				a=a+1	
				plot(do.call(grid.arrange,  GRAPH))
				ggsave(do.call(grid.arrange,  GRAPH),file=paste0(path_folder,"/output/trees_with_major_errors_",a,".pdf"),width = 29.7, height = 20.1, units = "cm")
				GRAPH = list() # reset plot
				i = 0 # reset index
			}}
	} # end of graph 
	
	if (length(ID)==0){ 
		print(paste0("No tree productivity above",maxrel, "% or below",-maxrel,"% of mean productivity at your plot. You may eventually want to try a lower threshold.")) 
	}
	if (output.errors & length(ID)>0){
		write.csv(DF[treeID%in%ID],file=paste0(path_folder,"/output/trees_with_major_errors.csv"))
	}
	return(DF)
} # end of major.error

#' Set mean productivity
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Set mean productivity (Mg/ha/yr) across all census intervals for a given site
#' @param DF a data.table
#' @param site provide the full name of your site (in lower case) i.e. 'barro colorado island'
#' @param strangler TRUE or FALSE, if TRUE, strangler figs tree are flagged (upon a list to published soon)
#' @param exclude.interval a vector (i.e. c(1,2)) indicating if a set of census intervals must be discarded from computation due for instance to a change in  protocol of measurment
#' @return mean productivity in Mg/ha/yr.
#' @export

determine.mean.prod <- function(DF,site,flag.stranglers,exclude.interval) {
	AREA <- site.info$size[match(site,site.info$site)]
	if (missing(exclude.interval)){
		ifelse(flag.stranglers,
			PRODA <- data.table(DF)[ficus!=1,.(prod=rec.flux(sum(agb1,na.rm=T),sum(agb2[code!="D"],na.rm=T),sum(agb1[code%in%c("A","AC")],na.rm=T),AREA,mean(int,na.rm=T))), by=interval],
			PRODA <- data.table(DF)[, .(prod=rec.flux(sum(agb1,na.rm=T),sum(agb2[code!="D"],na.rm=T),sum(agb1[code%in%c("A","AC")],na.rm=T),AREA,mean(int,na.rm=T))), by=interval])
	}	else {
		ifelse(flag.stranglers,
			PRODA <- data.table(DF)[ficus!=1 & !interval%in%c(exclude.interval),.(prod=rec.flux(sum(agb1,na.rm=T),sum(agb2[code!="D"],na.rm=T),sum(agb1[code%in%c("A","AC")],na.rm=T),AREA,mean(int,na.rm=T))), by=interval],
			PRODA <- data.table(DF)[!interval%in%c(exclude.interval),.(prod=rec.flux(sum(agb1,na.rm=T),sum(agb2[code!="D"],na.rm=T),sum(agb1[code%in%c("A","AC")],na.rm=T),AREA,mean(int,na.rm=T))), by=interval])
	}
	mPROD <- mean((PRODA$prod))
	return(mPROD)
}

#' Loess
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description a wrapper to get smoothed predictions of AGB fluxes using a loess function (library 'locfit')
#' @param x a data.table
#' @param var the name of the variable to be smoothed again intial AGB
#' @param range the range of initial AGB to be used for prediction (i.e. 5th and 95th percentiles of the whole distribution)
#' @return a smoothed prediction of the variable of interest
#' @export
loess.fun <- function(x,var,range,alpha=NULL)  {
	fit <- locfit(var ~ lAGB, data=x)
	pred <- predict(fit,newdata=list(lAGB=range))
	return(data.frame(lAGB=range,y=as.numeric(pred)))
}

# loess.fun <- function(x,var)  {
# 	fit <- locfit(var ~ lAGB, data=x)
# 	Xrange = as.numeric(quantile(x,seq(0.01,1,.01),na.rm=T)) # 100 values to estimate the smooth spline
# 	pred <- predict(fit,newdata=list(lAGB=Xrange))
# 	return(data.frame(lAGB=Xrange,y=as.numeric(pred)))
# }

#' Normalized tree status
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Check the consistency of stem/tree status over time (i.e. a tree that is 'alive' at last census can not be 'dead' inbetween)
#' @param x a data.table
#' @return a data.table (data.frame) with a new colum "status1" where values can be "P"(prior),"A"(alive),"D"(dead) and "Dr"(dead replicated).
#' @export
#' 

normal.stat <-function(DT) {
	if(!"status"%in%names(DT)) {
		DT$status <- NA
		DT$status[DT$DFstatus=="alive"] <- "A"
		DT$status[DT$DFstatus=="dead"] <- "D"
	}
	STAT <- rep("A",nrow(DT))
	if (all(is.na(DT$dbh))) {
		STAT <- rep("Dr",nrow(DT))
	} else if (any(is.na(DT$dbh))|any(grep("\\bD\\b",DT$status))) { # look for dbh=NA or dead (D) status
		locA <- which(DT$status=="A" & !is.na(DT$dbh))
		if (length(locA)!=0) {
			if (min(locA) >1) {
				STAT[is.na(DT$dbh)][1:min(locA)-1] <- "P" }  # Prior (to be discarded)
			if (max(locA) < nrow(DT)) {
				STAT[(max(locA)+1)] <- "D" }
			if (max(locA)+2 <= nrow(DT)) {
				STAT[(max(locA)+2):nrow(DT)] <- "Dr"	} # Dead replicated (to be discarded)
		} # end of loca>0
	}
	return(STAT)
}

#' Assign wood density
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Assign wood density using CTFS wood density data base (WSG)
#' @param df a data.table
#' @param site provide the full name of your site (in lower case) i.e. 'barro colorado island'
#' @param wsgdata a list of tree species (mnenomic) and corresponding wood density
#' @param denscol the variable to be return ("wsg" by default). No other option is implemented.
#' @return a data.table (data.frame) with all relevant variables.
#' @export

density.ind <- function (DF, site, wsgdata, denscol = "wsg") {
	if(!exists("DATA_path")){
		DATA_path <<- paste0(path_folder,"/data/")
	}
	#Add genus & species to data
	SP <-  LOAD(paste(DATA_path,list.files(DATA_path)[grep("spptable",list.files(DATA_path))],sep="/"))
	SP <- subset(SP,,c("sp","Genus","Species","Family"))
	SP$name = paste(SP$Genus,SP$Species,sep=" ")
	
	# Import & format CTFS wood data base
	wsgdatamatch = which(WSG$site %in% site.info$wsg.site.name[site.info$site==tolower(site)])
	if (length(wsgdatamatch) == 0)	stop("Site name doesn't match!")
	wsg = unique(WSG[wsgdatamatch,c("genus","species","spwood") ])
	names(wsg) <- c("genus","species","wd")
	wsg <- wsg[!is.na(wsg$wd),]
	
	# Assign WD by taxon
	A <- getWoodDensity(SP$Genus, SP$Species, stand = rep(site,nrow(SP)), family = NULL, region = "World",addWoodDensityData = wsg)
	A$name <- paste(A$genus,A$species,sep=" ")
	
	SP <- unique(merge(SP,A[,c("name","meanWD")],by="name",all.x=T))
	names(SP) <- c("name", "sp", "genus", "species", "Family", "wsg")
	if (any(grep("name",names(DF)))) {
		DT <- merge(DF,subset(SP,,c("name","wsg")),by="name",all.x=T)	
	} else {
		DT <- merge(DF,subset(SP,,c("sp","wsg")),by="sp",all.x=T)
	}
	# Allocate mean WD to species not in the list
	if (any(is.na(DT$wsg))){
	print(paste0("There are ",nrow(DT[is.na(wsg)])," individuals without WD values. Plot-average value (", round(mean(DT$wsg, na.rm=T),2),") was assigned."))
	DT <- within(DT,wsg[is.na(wsg)] <- mean(wsg, na.rm=T))
	}
	return(DT)
}

#' AGB computation
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Compute above-ground biomass using the updated model of Chave et al. (2014), given in Rejou-Mechain et al. (2017)
#' @param site provide the full name of your site (in lower case) i.e. 'barro colorado island'
#' @param D a column with tree diameters (in mm)
#' @param WD a column with wood density values
#' @param H a column with tree heights (optional)
#' @return a vector with AGB values (in Mg)
#' @export

AGB.comp <- function (site,D,WD, H = NULL) {
	if (length(D) != length(WD))
		stop("D and WD have different lenghts")
	if (is.null(site)) {
		stop("You must specified your site") }
	if (!is.null(H)) {
		if (length(D) != length(H))
			stop("H and WD have different length")
		if (any(is.na(H)) & !any(is.na(D)))
			warning("NA values are generated for AGB values because of missing information on tree height, \nyou may construct a height-diameter model to overcome that issue (see ?HDFunction and ?retrieveH)")
		if (any(is.na(D)))
			warning("NA values in D")
		AGB <- (0.0673 * (WD * H * (D/10)^2)^0.976)/1000
	}
	else {
		INDEX <- match(tolower(site),site.info$site)
		if (is.na(INDEX)) {			stop("Site name should be one of the following: \n",paste(levels(factor(site.info$site)),collapse=" - ")) }
		E <- site.info$E[INDEX]
		AGB <- exp(-2.023977 - 0.89563505 * E + 0.92023559 *
				log(WD) + 2.79495823 * log(D/10) - 0.04606298 *
				(log(D/10)^2))/1000
	}
	return(AGB)
}


#' Load object
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description a wrapper to softly load R objects in the Global environment
#' @param saveFile the path to the object to be loaded
#' @return import the object
#' @export

LOAD <- function(saveFile) {
	env <- new.env()
	load(saveFile, envir=env)
	loadedObjects <- objects(env, all=TRUE)
	stopifnot(length(loadedObjects)==1)
	env[[loadedObjects]]
}

#' Assign status
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Assign status alive ("A"), alive with POM changed ("AC"), dead ("D"), recruited ("R") or resprout ("Rsp") to trees, and check for consistency over time (i.e. avoid resurrection)
#' @param DF a data.table
#' @return update the column 'status1' with consistent information.
#' @export

assign.status <- function(DT) {
	code <- rep("A",nrow(DT))
	code[code=="A" & is.na(DT$dbhc1) & DT$status1=="P"] <- "R"
	code[code=="A" & DT$code2=="R"] <- "R"
	code[code=="A" & DT$code2=="broken"] <- "B" # Broken
	code[code=="A" & is.na(DT$status1) & !is.na(DT$status2)] <- "Rsp"  # resprouted trees poses a problem only the first year, are alive/dead after
	code[code=="A" & grep("MAIN",DT$code1) & grep("SPROUT",DT$code2)] <- "Rsp"  # resprouted trees poses a problem only the first year, are alive/dead after
	code[DT$status2=="D"] <- "D"  # all dead trees have NA value for DBH
	if(any(code=="A")){
		loc <- max(which(code%in%c("A")))
		if (any(code[1:loc]=="D")) {
			code[which(code[1:loc]=="D")] <- "A"
			code[DT$code2=="broken"] <- "B" # Broken
		}}
	dHOM <- DT$hom2 - DT$hom1
	dHOM[is.na(dHOM)] <- 0
	code[code=="A" & dHOM!=0] <- "AC" # trees with POM changed are not accounted for in productivity
	return(list(code,dHOM))
}

assign.status2 <- function(DT) {
	code <- rep("A",nrow(DT))
	code[DT$status1=="P"] <- "R"
	code[is.na(DT$dbhc1) & DT$status2=="A" & code!="R"] <- "Rsp"  # resprouted trees poses a problem only the first year, are alive/dead after
	code[DT$status2=="D"] <- "D"
	if(any(code=="A")){
		loc <- max(which(code%in%c("A","Rsp")))
		if (any(code[1:loc]=="D")) {
			code[which(code[1:loc]=="D")] <- "A"
		}}
	dHOM <- DT$hom2 - DT$hom1
	dHOM[is.na(dHOM)] <- 0
	code[code=="A" & dHOM!=0] <- "AC" # trees with POM changed are not accounted for in productivity
	return(list(code,dHOM))
}


#' Create quadrats
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Creat a grid where all trees are allocated to a given quadrat of size (= grid size).
#' @param census a data.frame where trees have relative X-Y coordinates.
#' @param grid_size the size of the grid (in meter)
#' @param x the identifier of X coordinates (i.e. 'gx')
#' @param y the identifier of Y coordinates (i.e. 'gy')
#' @param make.squared does the coordinates needs to 
#' @return add three columns to the data.frame: quadrat's number, centroids X and Y.
#' @export

create.quadrats=function(census,grid_size,x="gx",y="gy",fit.in.plot) {
	X <- census[,grep(x,names(census)),with=F][[1]]
	Y <- census[,grep(y,names(census)),with=F][[1]]
	
	if (any(is.na(X))|any(is.na(Y))) {
		warning(paste(length(X[is.na(X)])," trees without coordinates were discarded."))
		census <- census[!is.na(X) & !is.na(Y)]
		X <- census[,grep(x,names(census)),with=F][[1]]
		Y <- census[,grep(y,names(census)),with=F][[1]]
	}
	minx=0
	miny=0
	maxx=max(X,na.rm=T)
	maxy=max(Y,na.rm=T)
	
	if(exists("fit.in.plot")) {
		maxx <- round(max(X)/10)*10
		maxy <- round(max(Y)/10)*10
	}
	x1=X
	x1[X<=0]=0.1
	x1[X>=maxx]=maxx-0.1
	y1=Y
	y1[Y<=0]=0.1
	y1[Y>=maxy]=maxy-0.1
	
	# specify grid size for division into quadrats
	w_grid = ceiling(maxx)/grid_size;
	h_grid = ceiling(maxy)/grid_size;
	n_quadrat = w_grid*h_grid;
	
	# for now, only allow grid sizes that fit neatly
	if ( max(x1,na.rm=T) > maxx | min(x1,na.rm=T) < 0 | max(y1,na.rm=T) > maxy | min(y1,na.rm=T) < 0 )
	{
		stop('Some trees are outside the plot boundaries. Consider fit.in.plot=T')
	}
	if ( round(w_grid) != w_grid | round(h_grid) != h_grid )
	{
		stop('Plot width and height must be divisible by grid_size');
	}
	
	census$quad<-(ceiling((maxy-miny)/grid_size))*(floor((x1-minx)/grid_size))+(floor((y1-miny)/grid_size)+1)   ## identifiant de cellule unique 1->100 en partant de 0,0
	if ( max(census$quad,na.rm=T) != n_quadrat )
	{
		stop(paste('Quadrat numbering error: expected ',n_quadrat,' quadrats; got ',max(census$quad,na.rm=T),sep=''))
	}
	census$centroX<-(floor(x1/grid_size)*grid_size)+(grid_size/2)
	census$centroY<-(floor(y1/grid_size)*grid_size)+(grid_size/2)
	# census[,.(max(gx)-min(gx),max(gy)-min(gy)),by=quad]
	return(census)
}

#' Unbiased recruitment flux
#' @author Helene Muller-Laudau (hmullerlandau@gmail.com )
#' @description Compute unbiased recruitment rate (i.e. account for unmeasured recruitment)
#' @param A0 biomass at initial census.
#' @param A1 biomass of alive trees at final census.
#' @param S1 initial biomass of individuals that survived to time t
#' @param time cenusus interval in year
#' @return absolute recruitment flux in % per year
#' @export

rec.flux <- function(A0,A1,S1,area,time) {
	rec <- log(A1/S1)*(A1-A0)/(area*time*log(A1/A0))
	ifelse(rec==Inf,0,rec)
	return(rec)
}

#' Unbiased loss flux
#' @author Helene Muller-Laudau (hmullerlandau@gmail.com )
#' @description Compute unbiased recruitment rate (i.e. account for unmeasured recruitment)
#' @param A0 biomass at initial census.
#' @param A1 biomass of alive trees at final census.
#' @param S1 initial biomass of individuals that survived to time t
#' @param time cenusus interval in year
#' @return absolute recruitment flux in % per year
#' @export

loss.flux <- function(A0,A1,S1,area,time) {
	LO <- (log(A0/S1)/log(A1/A0))*((A1-A0)/(area*time))
	ifelse(LO==Inf,0,LO)
	return(LO)
}
