#' CTFS-formated data preparation
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Main routine to format and correct CTFS-formated data
#' @param site provide the full name of your site (in lower case) i.e. 'barro colorado island'
#' @param stem TRUE or FALSE, are you using the stem data (stem=TRUE) or tree data (i.e. called 'full')?
#' @param taper.correction TRUE or FALSE, are you willing to apply Cushman et al (2014) taper correction?
#' @param fill_missing TRUE or FALSE, are you willing to extrapolate missing DBH from surrounding DBH?
#' @param palm TRUE or FALSE, if TRUE, biomass of palm trees is computed through a specific allometric model (Goodman et al. 2013)
#' @param strangler TRUE or FALSE, if TRUE, strangler figs tree are flagged (upon a list to published soon)
#' @param maxrel a numeric value setting the threshold over which relative productivity is assumed to be too high (usually set at 20 percents)
#' @param draw.graph TRUE or FALSE, if TRUE, draw graph for trees with 'erroneous' DBH measures/productivity
#' @param output.errors TRUE or FALSE, if TRUE, creat a CSV file with all trees with erroneous DBH measures/productivity
#' @param DATA_path allows to provide a different path where the data are located
#' @param exclude.interval a vector (i.e. c(1,2)) indicating if a set of census intervals must be discarded from computation due for instance to a change in  protocol of measurment
#' @return a data.table (data.frame) with all relevant variables.
#' @export

data.prep <- function(site,stem,taper.correction,fill_missing,palm,strangler,maxrel=20,draw.graph,output.errors,DATA_path,exclude.interval=exclude.interval) {
	site <- tolower(site)
	INDEX <- match(tolower(site),site.info$site)
	if (is.na(INDEX)) {			stop("Site name should be one of the following: \n",paste(levels(factor(site.info$site)),collapse=" - ")) }

	if(missing(DATA_path)){
		DATA_path <- paste0(path_folder,"/data/")
	}

	files <-list.files(DATA_path)
	ifelse(stem,files <- files[grep("stem",files)], files <- files [grep("full",files)])
	ifelse(!dir.exists(file.path(paste0(path_folder,"/output"))), dir.create(file.path(paste0(path_folder,"/output"))), FALSE)


	# Create the receiving data.frame
	df <- data.frame("treeID"=NA, "stemID"=NA, "tag"=NA, "StemTag"=NA, "sp"=NA, "quadrat"=NA, "gx"=NA, "gy"=NA,"dbh"=NA,"hom"=NA, "ExactDate"=NA, "DFstatus"=NA, "codes"=NA, "date"=NA, "status"=NA,"CensusID"=NA)

	for (i in 1:length(files)) {
		temp <- LOAD(paste(DATA_path,files[i],sep="/"))
		temp$CensusID <- i
		temp  <- temp[match(names(df),names(temp))]
		df <- rbind(df,temp)
	}
	rm(temp)
	df <- data.table(df[-1,])

	df <- data.correction(df,taper.correction,fill_missing)
	print("Step 1: data correction done.")

	df <- computeAGB(df,site,palm,DATA_path)
	print("Step 2: AGB calculation done.")

	DF <- format.interval(df,strangler)
	print("Step 3: data formating done.")

	DF <- flag.errors(DF,site,strangler=strangler,maxrel=maxrel,draw.graph=draw.graph,output.errors=output.errors,exclude.interval=exclude.interval)
	print("Step 4: errors flagged. Saving corrected data into /output folder.")

	save(DF,file=paste0(path_folder,"/output/",site,"_census_intervals.Rdata"))

	rm(list=setdiff(ls(), c("DF","path_folder","SITE", lsf.str())))
	return(DF)
}

#' Trim and correct data
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Stack all censuses together and correct DBH, if required
#' @param taper.correction TRUE or FALSE, are you willing to apply Cushman et al (2014) taper correction?
#' @param fill_missing TRUE or FALSE, are you willing to extrapolate missing DBH from surrounding DBH?
#' @return a data.table (data.frame) with all relevant variables.
#' @export

data.correction <- function(df,taper.correction,fill_missing) {
	df <- df[!status%in%c("P"),] # discard all priors
	df[,id :=paste(df$treeID,df$stemID,sep="-")] # creat a unique tree-stem ID
	df <- df[order(id,CensusID)]
	df[,status1:=normal.stat(.SD),by=id]

	df[status1=="M", nrow2 := seq_len(.N), by = id]
	df <- within(df,nrow2[is.na(nrow2)] <- 0)
	df <- within(df,status1[nrow2==1] <- "D") # multiple missing trees are considered as dead at first occurrence
	df <- df[nrow2<2]  # keeps only 1 line for dead trees

	df[status1=="D", nrow := seq_len(.N), by = id]
	df <- within(df,nrow[is.na(nrow)] <- 0)
	df <- df[nrow<2,]  # keeps only 1 line for dead trees
	df[,c("nrow","nrow2"):=NULL]

	df[,year:=round(mean(as.numeric(substr(ExactDate,1,4)),na.rm=T)),by=CensusID] # Assign 1 year per census
	NO.MEASURE <- df[,all(is.na(dbh)),by=id]
	df <- df[!id%in%NO.MEASURE$treeID[NO.MEASURE$V1]] # remove stems without any measure

	# Taper correction or missing values: -> fill gaps for missing values
	df[, c("dbh2","hom2") := corDBH(.SD,taper.correction=taper.correction,fill_missing=fill_missing), by=id] # might be time consuming (4 minutes for BCI)
	table(is.na(df$dbh2),df$status1)
	NO.MEASURE <- df[,all(is.na(dbh2)),by=treeID]
	table(NO.MEASURE$V1)
	df <- df[!treeID%in%NO.MEASURE$treeID[NO.MEASURE$V1]] # remove trees without any measurement
	return(df)
}

#' Biomass computation
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Allocate wood density and compute above-ground biomass using the updated model of Chave et al. (2014), given in Rejou-Mechain et al. (2017). Palm trees (palm=T) are computed using a different allometric model (Goodman et al. 2013).
#' @param site provide the full name of your site (in lower case) i.e. 'barro colorado island'
#' @param palm TRUE or FALSE, if TRUE, biomass of palm trees is computed through a specific allometric model (Goodman et al. 2013)
#' @param DATA_path allows to provide a different path where the data are located
#' @return a data.table (data.frame) with all relevant variables.
#' @export

computeAGB <- function(df,site,palm=T,DATA_path) {
	## Allocate wood density
	df$wsg <- density.ind(df=df,site,wsg=WSG)

	# Compute biomass
	df$agb <- AGB.comp(site,df$dbh2, df$wsg,H = NULL)

	# Compute biomass for palms
	if (palm) {
		SP <-  LOAD(paste(DATA_path,list.files(DATA_path)[grep("spptable",list.files(DATA_path))],sep="/"))
		if(is.na(match("genus",tolower(names(SP))))) {
			trim <- function (x) gsub("^\\s+|\\s+$", "", x)
			SP$genus <-  trim(substr(SP$Latin,1,regexpr(" ",SP$Latin)))
			SP$species <-  trim(substr(SP$Latin,regexpr(" ",SP$Latin),50))
			SP <- SP[,c("sp","genus","species","Family")]
			SP$name <- paste(SP$genus,SP$species,sep=" ")
			names(SP) <- c("sp","genus","species","Family","name")
			df <- merge(df,SP,by="sp",all.x=T)
			agbPalm <- function(D) { exp(-3.3488 + 2.7483*log(D/10) + ((0.588)^2)/2)/1000 }
			df['Family'=="Arecaceae",agb:= agbPalm(dbh2)]
		} else {
			SP <- SP[,c("sp","Genus","Species","Family")]
			SP$name <- paste(SP$Genus,SP$Species,sep=" ")
			names(SP) <- c("sp","genus","species","Family","name")
			df <- merge(df,SP,by="sp",all.x=T)
			agbPalm <- function(D) { exp(-3.3488 + 2.7483*log(D/10) + ((0.588)^2)/2)/1000 }
			df['Family'=="Arecaceae",agb:= agbPalm(dbh2)]
		}
	}
	return(df)
}

#' Format census intervals
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Create census intervals (i.e. put consecutive census side by side), assign status by tree (i.e. alive (A),dead (D), recruited (R) or resprout (Rsp))
#' @param df a data.table
#' @param strangler TRUE or FALSE, if TRUE, strangler figs tree are flagged (upon a list to published soon)
#' @return a formated data.table.
#' @export

format.interval <- function(df,strangler) {
	# Receiveing data set
	DF <- data.table("treeID"=NA,"dbh1"=NA,"dbhc1"=NA,"status1"=NA,"code1"=NA,"hom1"=NA,"agb1"=NA,"date1"=NA,"dbh2"=NA,"dbhc2"=NA,"status2"=NA,"code2"=NA,"hom2"=NA,"agb2"=NA,"agbl"=NA,"date2"=NA,"interval"=NA,"year"=NA)

	YEAR <- levels(factor(df$year))
	for (j in 1:(length(YEAR)-1)) {  # 4 minutes to run
		i1 <- df[year==YEAR[j] & status1 != "D", .I[which.max(dbh2)], by = treeID] # keep only information for the biggest stem per treeID
		A1 <- df[i1$V1,c("treeID","dbh","dbh2","status1","codes","hom")]
		names(A1) <- c("treeID","dbh1","dbhc1","status1","code1","hom1")
		B1 <- df[year==YEAR[j] & status1 != "D",list("agb1"=sum(agb,na.rm=T),"date1"=mean(date,na.rm=T)),by=treeID]
		BB <- merge(B1,A1,by="treeID",all.x=T)
		cens1 <- BB[,c("treeID","dbh1","dbhc1","status1","code1","hom1","agb1","date1")]


		i2 <- df[year==YEAR[j+1] & status1 != "D", .I[which.max(dbh2)], by = treeID]
		A2 <- df[i2$V1,c("treeID","dbh","dbh2","codes","hom","status1")]
		names(A2) <- c("treeID","dbh2","dbhc2","code2","hom2","status2")

		B2 <- df[year==YEAR[j+1],list("agb2"=sum(agb[status1!="D"],na.rm=T),"agbl"=sum(agb[status1=="D"],na.rm=T),"date2"=mean(date,na.rm=T)),by=treeID]
		BB <- merge(B2,A2,by="treeID",all.x=T)
		cens2 <- BB[,c("treeID","dbh2","dbhc2","status2","code2","hom2","agb2","agbl","date2")]
		cens2 <- within(cens2,status2[is.na(status2)] <- "D")

		ID <- data.table(treeID=unique(c(cens1$treeID,cens2$treeID)))
		ID <- merge(ID,cens1,by='treeID',all.x=T)
		ID <- merge(ID,cens2,by='treeID',all.x=T)
		ID$interval=j
		ID$year=YEAR[j+1]
		DF <- rbind(DF,ID)
	}
	DF <- DF[-1,]

	# Add coordinates & other mandatory information
	COORD <- df[,.(gx=round(10*gx[!is.na(gx)])/10,gy=round(10*gy[!is.na(gy)])/10,quadrat=quadrat,name=name,sp=sp),by=treeID]
	COORD <- unique(COORD[!duplicated(COORD$treeID)])
	DF <- merge(DF,COORD,by="treeID",all.x=T)

	# Add average date of census when missing
	DATE <- DF[,.(date1=mean(date1,na.rm=T),date2=mean(date2,na.rm=T)),by=year]
	DF <- within(DF,date1[is.na(date1)] <- DATE$date1[match(DF[is.na(date1),"year"]$year,DATE$year)])
	DF$int <- (DF$date2 - DF$date1)/365.5  # census interval in days


	# Update status for recruited trees
	DF[, nrow := seq_len(.N), by = treeID]
	DF <- within(DF,status1[is.na(status1) & !is.na(dbh2) & nrow==1] <- "P")  # recruited trees or Z code (=?)
	DF[,nrow:= NULL]

	# Assign status
	DF[, c("code","dHOM") := assign.status(.SD), by=treeID]

	# Remove unnecessary rows
	DF[code=="D", nrow := seq_len(.N), by = treeID]
	DF <- within(	DF,nrow[is.na(nrow)] <- 0)
	DF <- 	DF[nrow<2,]  # keeps only 1 line for dead trees
	DF[,"nrow":=NULL]

	# Compute annualized fluxes
	DF[code%in%c("A","AC"),prod.g := (agb2-agb1)/int,by=treeID] # annual prod for alive trees
	DF[code%in%c("R","Rsp"),prod.r:=agb2/int,by=treeID] # annual prod for resprouts and recruits
	DF[,loss:=agbl/int,by=treeID] # annualized loss for dead trees
	DF[code=="D",loss2:=agb1/int,by=treeID]
	# DF[,.(prod.g=mean(prod.g,na.rm=T),prod.r=mean(prod.r,na.rm=T),loss=mean(loss,na.rm=T),agbl=mean(agbl,na.rm=T)),by=code]
	# Flag strangler figs
	if(strangler) {
		DF$ficus <- 0
		ficus$name <- paste(ficus$Genus,ficus$Species,sep=" ")
		FIC <- match(DF$name,ficus$name)
		DF <- within(DF,ficus[!is.na(FIC)]<-1)
	}
	return(DF)
}

#' Flag major errors
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Identify trees with major errors in DBH measurments. A major error correspond to a relative individal productivity (or growth) is above a given percentage (set by 'maxrel') of the mean productivity computed at a site. Additionnaly, flagged trees that died at next census interval are also flagged. Option to see DBH measurement (=draw.graph) of flagged trees or print a csv (output.errors) are given.
#' @param DF a data.table
#' @param site provide the full name of your site (in lower case) i.e. 'barro colorado island'
#' @param strangler TRUE or FALSE, if TRUE, strangler figs tree are flagged (upon a list to published soon)
#' @param maxrel a numeric value setting the threshold over which relative productivity is assumed to be too high (usually set at 20 percents)
#' @param draw.graph TRUE or FALSE, if TRUE, draw graph for trees with 'erroneous' DBH measures/productivity
#' @param output.errors TRUE or FALSE, if TRUE, creat a CSV file with all trees with erroneous DBH measures/productivity
#' @param exclude.interval a vector (i.e. c(1,2)) indicating if a set of census intervals must be discarded from computation due for instance to a change in  protocol of measurment
#' @return a data.table (data.frame) with all relevant variables.
#' @export
#'
flag.errors <- function(DF,site,strangler,maxrel,draw.graph,output.errors,exclude.interval) {
	mean.prod <- determine.mean.prod(DF,site,strangler,exclude.interval)
	DF[,prod.rel:=as.numeric(NA),]
	DF[,prod.rel:=prod.g*100/mean.prod] # relative contribution to average total productivity
	DF[,error.prod:=0]
	DF <- within(DF,error.prod[prod.rel>maxrel & dHOM==0 & code!="D"] <- 1)
	DF <- within(DF,error.prod[prod.rel<(-maxrel) & dHOM==0 & code!="D"] <- -1)
	table(DF$error.prod,DF$year)

	# Filter for AGB loss: discard trees that have major errors at previous census of year of death
	DF[,error.loss:=0]
	ID.loss <- DF[,treeID[.I[code=="D"]][error.prod[.I[code=="D"]-1]!=0]]
	DF <- within(DF,error.loss[treeID%in%ID.loss & code=="D"] <- 1)
	table(DF$error.loss,DF$year)
	ID <- DF[error.prod!=0 & code!="D",nrow(.SD)>1,by=treeID]

	if(draw.graph) { # Plot trees with large major error
		if(nrow(ID[ID$V1])>0 & nrow(ID[ID$V1])<20) {
			YEAR <- levels(factor(DF$year))
			CX=2

			DF[,year:=as.numeric(year)]
			X <- DF[treeID%in%ID[V1==T,treeID] & !code%in%c("R") ][order(treeID,year)]
			X$point <- 0
			X$point[X$error.prod!=0] <- 1
			Y <- DF[treeID%in%ID[V1==T,treeID] & !code%in%c("D","R") ][order(treeID,year)]
			YY <- Y[,.(year=max(year),name=unique(name),d2=dbhc2[year==max(year)],d02=dbh2[year==max(year)],hom2=hom2[year==max(year)]),by=treeID]
			Y$line <- 0
			Y$line[Y$dHOM==0] <- 1
			Y$point <- 0
			Y$point[Y$error.prod!=0] <- 1

			A <- ggplot(X,aes(x=year,y=dbhc1)) + geom_point(size=2) + facet_wrap(~treeID+name,scale="free",ncol=3,labeller = label_wrap_gen(multi_line=FALSE,width=50))  + geom_segment(data=Y,aes(x=year,y=dbhc1,xend=year+5,yend=dbhc2,linetype=as.factor(line))) + geom_point(data=X[point==1],aes(x=year+5,y=dbhc2),col=2)	+ labs(x="Year",y="dbh (mm)") + geom_text(data=Y,aes(x=year,y=dbh1-(0.05*dbh1)),label=round(Y$hom1,2),cex=CX)	+ geom_text(data=YY,aes(x=year+5,y=d02-(0.05*d02)),label=round(YY$hom2,2),cex=CX) + geom_text(data=Y,aes(x=year,y=dbh1-(0.05*dbh1)),label=round(Y$hom1,2),cex=CX) + geom_text(data=Y,aes(x=year,y=0.2*max(dbh2)),label=Y$dbh1,cex=CX,angle=90,vjust=1) + geom_text(data=YY,aes(x=year+5,y=0.2*max(d2)),label=YY$d02,cex=CX,angle=90,vjust=1) + theme(plot.title = element_text(hjust = 0.5,size=2*CX,face="bold"),axis.text.y = element_text(size=4*CX),axis.text.x = element_text(size=4*CX,vjust=0,angle=30),panel.background = element_blank(),strip.text = element_text(size = 4*CX,face="bold"),strip.background = element_rect("lightgrey"),panel.spacing = unit(0.1, "lines")) + scale_linetype_manual(values=c('0'="dashed",'1'="solid")) + guides(linetype=F,colour=F) + scale_x_continuous(limits=c(min(as.numeric(YEAR))-3,max(as.numeric(YEAR))+3))
			print(A)
			ggsave(A,file=paste0(path_folder,"/output/trees_with_major_errors.pdf"),width = 21, height = 29.7, units = "cm")
		} else if (nrow(ID[ID$V1])>20)
		{print(paste0("Too many graphs (N=",nrow(ID[ID$V1]),") to be drawned"))
		} else
		{ print(paste0("No tree productivity above",maxrel, "% or below",-maxrel,"% of mean productivity at your plot. You may eventually try a lower threshold.")) }
	} # end of print
	if (output.errors & nrow(ID[ID$V1])>0){
		write.csv(DF[treeID%in%ID$treeID],file=paste0(path_folder,"/output/trees_with_major_errors.csv"))
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

determine.mean.prod <- function(DF,site,strangler,exclude.interval) {
	AREA <- site.info$size[match(site,site.info$site)]
	if (missing(exclude.interval)){
		ifelse(strangler,PRODA <- DF[ficus!=1, .(prod=sum(prod.g,na.rm=T)), by=interval],PRODA <- DF[, .(prod=sum(prod.g,na.rm=T)), by=interval])
	}	else {
		ifelse(strangler,PRODA <- DF[ficus!=1 & !interval%in%c(exclude.interval), .(prod=sum(prod.g,na.rm=T)), by=interval],PRODA <- DF[!interval%in%c(exclude.interval), .(prod=sum(prod.g,na.rm=T)), by=interval])
	}
	mPROD <- mean((PRODA$prod)/AREA)
	return(mPROD)
}

#' CTFS-formated data preparation
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Main routine to format and correct CTFS-formated data
#' @param site provide the full name of your site (in lower case) i.e. 'barro colorado island'
#' @param stem TRUE or FALSE, are you using the stem data (stem=TRUE) or tree data (i.e. called 'full')?
#' @param taper.correction TRUE or FALSE, are you willing to apply Cushman et al (2014) taper correction?
#' @param fill_missing TRUE or FALSE, are you willing to extrapolate missing DBH from surrounding DBH?
#' @param palm TRUE or FALSE, if TRUE, biomass of palm trees is computed through a specific allometric model (Goodman et al. 2013)
#' @param strangler TRUE or FALSE, if TRUE, strangler figs tree are flagged (upon a list to published soon)
#' @param maxrel a numeric value setting the threshold over which relative productivity is assumed to be too high (usually set at 20 percents)
#' @param draw.graph TRUE or FALSE, if TRUE, draw graph for trees with 'erroneous' DBH measures/productivity
#' @param output.errors TRUE or FALSE, if TRUE, creat a CSV file with all trees with erroneous DBH measures/productivity
#' @param DATA_path allows to provide a different path where the data are located
#' @param exclude.interval a vector (i.e. c(1,2)) indicating if a set of census intervals must be discarded from computation due for instance to a change in  protocol of measurment
#' @return a data.table (data.frame) with all relevant variables.
#' @export


#' Loess
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description a wrapper to get smoothed predictions of AGB fluxes using a loess function (library 'locfit')
#' @param x a data.table
#' @param var the name of the variable to be smoothed again intial AGB
#' @param range the range of initial AGB to be used for prediction (i.e. 5th and 95th percentiles of the whole distribution)
#' @return a smoothed prediction of the variable of interest
#' @export
loess.function <- function(x,var,range)  {
	fit <- locfit(var ~ lAGB, data=x)
	pred <- predict(fit,newdata=list(lAGB=range))
	return(as.numeric(pred))
}

#' Normalized tree status
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Check the consistency of stem/tree status over time (i.e. a tree that is 'alive' at last census can not be 'dead' inbetween)
#' @param x a data.table
#' @return a data.table (data.frame) with all relevant variables.
#' @export

normal.stat <-function(X) {
	STAT <- X$status
	if (any(is.na(X$dbh))|any(grep("D",X$status))) {
		locA <- which(X$status=="A" & !is.na(X$dbh))
		if (length(locA)!=0) {
			if(any(grep("D",STAT[1:min(locA)]))) {
				STAT[is.na(X$dbh)][1:min(locA)-1] <- "P" }
			if(any(grep("D",STAT[min(locA):max(locA)]))) { #
				STAT[min(locA):max(locA)] <- "A" }
			if (all(is.na(X$dbh[(max(locA)+1):nrow(X)]))) {
				STAT[(max(locA)+1):nrow(X)] <- "D" }
		}
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

density.ind <- function (df, site, wsgdata, denscol = "wsg") {
	wsgdatamatch = which(wsgdata$site %in% site.info$wsg.site.name[site.info$site==site])
	if (length(wsgdatamatch) == 0)
		stop("Site name doesn't match!")
	wsgdata = unique(wsgdata[wsgdatamatch, ])
	meanwsg = mean(subset(wsgdata, idlevel == "species")[, denscol],
		na.rm = TRUE)
	m = match(df$sp, wsgdata$sp)
	result = wsgdata[m, denscol]
	result[is.na(m)] = meanwsg
	result[is.na(result)] = meanwsg
	return(result)
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

#' Data correction
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Perform two mains tasks: (a) apply a taper correction when POM is > 130 cm, and (b) linear interpolation values when missing DBHs. Interpolation of missing values is done by averaging surrounding available DBH values.
#' @param DF a data.table
#' @param taper.correction TRUE or FALSE, are you willing to apply Cushman et al (2014) taper correction?
#' @param fill_missing TRUE or FALSE, are you willing to extrapolate missing DBH from surrounding DBH?
#' @return a data.table (data.frame) with all relevant variables.
#' @export

# Correction of DBH
corDBH <- function(DF,taper.correction,fill_missing) {
	hom2 <-  round(as.numeric(DF$hom)*100)/100
	dbh2 <-  DF$dbh

	if (!all(is.na(dbh2))) { # if all DBH are NA -> can't do much -> much be discarded
		if (fill_missing & any(is.na(dbh2)))  { # Feel gap for missing DBH with average between prior and neDFt census
			loc <- which(is.na(dbh2))
			if (any(grepl("R",DF$codes))) { # avoid resprouts
				RESP <- which(grepl("R",DF$codes))
				if (loc[1]!=min(RESP)) {
					if(any(loc==1)) {
						M <- matrix(c(NA,dbh2[loc-1],dbh2[loc+1]),nrow=length(loc),2)
						M2 <- matrix(c(NA,hom2[loc-1],hom2[loc+1]),nrow=length(loc),2)
					} else {
						M <- matrix(c(dbh2[loc-1],dbh2[loc+1]),nrow=length(loc),2)
						M2 <- matrix(c(hom2[loc-1],hom2[loc+1]),nrow=length(loc),2)}
					dbh2[is.na(dbh2)] <- apply(M,1,mean,na.rm=T)
					hom2[is.na(hom2)] <- apply(M2,1,mean,na.rm=T)

					while(any(is.na(dbh2))) {
						loc <- which(is.na(dbh2))
						if(any(loc==1)) {
							M <- matrix(c(NA,dbh2[loc-1],dbh2[loc+1]),nrow=length(loc),2)
							M2 <- matrix(c(NA,hom2[loc-1],hom2[loc+1]),nrow=length(loc),2)
						} else {
							M <- matrix(c(dbh2[loc-1],dbh2[loc+1]),nrow=length(loc),2)
							M2 <- matrix(c(hom2[loc-1],hom2[loc+1]),nrow=length(loc),2)}
						dbh2[is.na(dbh2)] <- apply(M,1,mean,na.rm=T)
						hom2[is.na(hom2)] <- apply(M2,1,mean,na.rm=T)
					}
				}}# end of resprout

			if(any(loc==1)) {
				M <- matrix(c(NA,dbh2[loc-1],dbh2[loc+1]),nrow=length(loc),2)
				M2 <- matrix(c(NA,hom2[loc-1],hom2[loc+1]),nrow=length(loc),2)
			} else {
				M <- matrix(c(dbh2[loc-1],dbh2[loc+1]),nrow=length(loc),2)
				M2 <- matrix(c(hom2[loc-1],hom2[loc+1]),nrow=length(loc),2)}
			dbh2[is.na(dbh2)] <- apply(M,1,mean,na.rm=T)
			hom2[is.na(hom2)] <- apply(M2,1,mean,na.rm=T)

			while(any(is.na(dbh2))) {
				loc <- which(is.na(dbh2))
				if(any(loc==1)) {
					M <- matrix(c(NA,dbh2[loc-1],dbh2[loc+1]),nrow=length(loc),2)
					M2 <- matrix(c(NA,hom2[loc-1],hom2[loc+1]),nrow=length(loc),2)
				} else {
					M <- matrix(c(dbh2[loc-1],dbh2[loc+1]),nrow=length(loc),2)
					M2 <- matrix(c(hom2[loc-1],hom2[loc+1]),nrow=length(loc),2)}
				dbh2[is.na(dbh2)] <- apply(M,1,mean,na.rm=T)
				hom2[is.na(hom2)] <- apply(M2,1,mean,na.rm=T)
			}}# end of missing loop
		if (any(is.na(hom2))) { hom2[is.na(hom2)] <- 1.3}

		# Apply Cushman's correction to trees with POM changed
		if (taper.correction & any(hom2 > 1.3)) {
			ifelse(any(grepl("R",DF$codes)),ind1 <- 1:(grep("R",DF$codes)[1]-1),ind1 <-  which(DF$status!="D"))
			dbh2[ind1] <- round(dbh2[ind1]*exp(0.0247*(hom2[ind1]-1.3)),1)
			# dbh2[DF$status=="D"] <- tail(dbh2[DF$status=="A" & !is.na(dbh2)],1)  # replicate last dbh to dead trees
			hom2 <- rep(1.3,nrow(DF))
		}}

	return(list(dbh2,hom2))
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
assign.status <- function(DF) {
	code <- rep("A",nrow(DF))
	code[is.na(DF$dbh1) & DF$status1=="P"] <- "R"
	code[is.na(DF$status1) & !is.na(DF$status2)] <- "Rsp"  # resprouted trees poses a problem only the first year, are alive/dead after
	code[is.na(DF$dbhc2)] <- "D"
	if(any(code=="A")){
		loc <- max(which(code%in%c("A","Rsp")))
		if (any(code[1:loc]=="D")) {
			code[which(code[1:loc]=="D")] <- "A"
		}}
	dHOM <- DF$hom2 - DF$hom1
	dHOM[is.na(dHOM)] <- 0
	code[code=="A" & dHOM!=0] <- "AC" # trees with POM changed are not accounted for in productivity
	return(list(code,dHOM))
}

#' Create quadrats
#' @author Ervan Rutishauser (er.rutishauser@gmail.com)
#' @description Creat a grid where all trees are allocated to a given quadrat of size (= grid size).
#' @param census a data.frame where trees have relative X-Y coordinates.
#' @param grid_size the size of the grid (in meter)
#' @return add three columns to the data.frame: quadrat's number, centroids X and Y.
#' @export

create.quadrats=function(census,grid_size) {
	X <- census[,grep("x",names(census)),with = FALSE][[1]]
	Y <- census[,grepl("y\\b",names(census)),with = FALSE][[1]]
	if (any(is.na(X))){
		warning(paste(length(X[is.na(X)])," trees without coordinates were discarded."))
		census <- census[!is.na(X)]
		X <- census[,grep("x",names(census)),with = FALSE][[1]]
		Y <- census[,grepl("y\\b",names(census)),with = FALSE][[1]]
	}
	minx=0
	miny=0
	maxx=max(X,na.rm=T)
	maxy=max(Y,na.rm=T)
	x1=X
	x1[X<=0]=0.1
	x1[X==maxx]=maxx-0.1
	y1=Y
	y1[Y<=0]=0.1
	y1[Y==maxy]=maxy-0.1

	# specify grid size for division into quadrats
	w_grid = ceiling(maxx)/grid_size;
	h_grid = ceiling(maxy)/grid_size;
	n_quadrat = w_grid*h_grid;

	# for now, only allow grid sizes that fit neatly
	if ( round(w_grid) != w_grid | round(h_grid) != h_grid )
	{
		stop('Plot width and height must be divisible by grid_size');
	}

	if ( max(X,na.rm=T) > maxx | min(X,na.rm=T) < 0 | max(Y,na.rm=T) > maxy | min(Y,na.rm=T) < 0 )
	{
		stop('Some trees are outside the plot boundaries')
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
