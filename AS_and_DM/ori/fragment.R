		#project-specific code; the bed.ids were prepared in 
		contrast<-logical(length(bed.ids))
		contrast[grep('HN',bed.ids)]<-TRUE
		norm.no<-length(which(!contrast))
		tumor.no<-length(which(contrast))
