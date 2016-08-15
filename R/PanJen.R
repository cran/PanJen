## read in packages
library(mgcv)
library(RColorBrewer)
library(Formula)
library(lasso2)

choose.fform <-function(data,variable,forms_to_test, distribution=gaussian, smoothing_splines=20){
  ########################################
  ###   convert input
  ########################################
  Forms=Formula(forms_to_test) ## turn formula-object into Formula-object
  data$x=data$variable=get(variable,data)
  nforms<-length(Forms)[2] ## number of functional forms+base
  dist<-distribution
  smo<-smoothing_splines

  ########################################
  ###   VALIDATION
  ########################################
  if (nforms<1) {
    stop("missing functional forms - please use 'fform()'or see ??PanJen() for how to defince user-functions ")
  }
  
  if (dim(data)[1]<dim(data)[1]) {
    stop("At least one functional form produces NaNs")
  }
  if (!smo>0) {
    print("Number of smoothing splines must be larger than 0 - Smoothing set to 20 splines by default")
  }
      
  ########################################
  ###   creating table and lists for output 
  ########################################
  frameNames=c("smoothing","base")
  for (j in 2:nforms){
    frameNames=c(frameNames, as.character(formula(Forms, rhs =j, lhs=0))[2])
    }
  tableOut=matrix(nrow=nforms+1,ncol=3)
  rownames(tableOut)=frameNames
  colnames(tableOut)= c("AIC", "BIC","ranking (AIC)")
  
  
  ########################################
  ###   Estimate models and save AIC + BIC
  ########################################  
  ## Estimate base
  baseForm=formula(Forms, rhs =1,collapse=T)
  base<-gam(baseForm, 
            family=dist,
            method="GCV.Cp",
            data=data)
  ## save for output
  tableOut[2,1]=AIC(base)
  tableOut[2,2]=BIC(base)
  models <- list(base)
  
  ## estimate smoothing
  equation_update<-
    update.formula(baseForm,
                   paste(paste(". ~ .+s(", all.vars(formula(Forms, rhs=2, lhs=0)),sep=""),
                         paste(paste("k=",smo,sep=),")"),sep=",") )
  smoothing<-gam(equation_update, 
                 family=dist,
                 method="GCV.Cp",
                 data=data)
  tableOut[1,1]=AIC(smoothing)
  tableOut[1,2]=BIC(smoothing)
  models <- c(models, list(smoothing))
  
  ## estimate one model for each transformation
  namesLL=c("model_base","model_smoothing") ## model names
  
  
  for (j in 1:(nforms-1)){    
    equation_update<-merge.formula(baseForm,formula(Forms, lhs=0, rhs = c(j+1)))
    model_c<-gam(equation_update, 
                 family=dist,
                 method="GCV.Cp",
                 data=data)
    tableOut[j+2,1]=AIC(model_c)
    tableOut[j+2,2]=BIC(model_c)
    models <- c(models, list(model_c))
    namesLL<-c(namesLL,paste("model",frameNames[j+2]
                             , sep="_"))
  }
  names(models)<-namesLL
  
  
  ## sort table after AIC
  tableOut[,3]=as.numeric(rank(tableOut[,1]))
  tableOut=tableOut[order(tableOut[,3]),]             
  
  ## sort models after BIC  
  order=paste("model",rownames(tableOut), sep="_")
  models<- models[order]
  #####################################################################
  ## return list  
  #####################################################################
  print(tableOut)   ## print table
  class(data)<-"data.frame"
  
  output=list("models"=models,
              "dataset"=data,
              "rank.table"=tableOut,
              "fforms"=Forms)
  class(output)<-"PJ"
  return(output) 
}


fform <-function(data,variable,base_form, distribution=gaussian, smoothing_splines=20){
  
  ## Predefined forms
  formsPack<-Formula(y~I(1/variable^2)|I(1/variable)|I(log(variable))|I(variable^0.5)|variable|I(variable^2)|variable+I(variable^2)) 
  ########################################
  ###   convert input
  ########################################
  Forms<-Formula(base_form) ## turn formula-object into Formula-object
  data$x=data$variable=get(variable,data)
  nforms<-length(Forms)[2] ## check that a variable is provided
  dist<-distribution
  smo <-smoothing_splines  
  
  ########################################
  ###   VALIDATION
  ########################################
  if (nforms<1) {
    stop("missing variabel to test - please see ??PanJen() for syntax")
  }
  if (!smo>0) {
    smo<-20
    print("Number of smoothing splines must be larger than 0 - Smoothing set to 20 splines by default")
  }
  
  
  ########################################
  ###   creating table and lists for output 
  ########################################
  frameNames<-c("smoothing","base", "1/x^2","1/x","log(x)","x^0.5","x", "x^2", "x+x^2")
  tableOut<-matrix(nrow=9,ncol=3)
  rownames(tableOut)<-frameNames
  colnames(tableOut)<-c("AIC", "BIC","ranking (AIC)")
  
  
  ########################################
  ###   Estimate models and save AIC + BIC
  ########################################  
  ## Estimate base
  baseForm<-formula(Forms, rhs =1,collapse=T)
  base<-gam(baseForm, 
            family=dist,
            method="GCV.Cp",
            data=data)
  ## save for output
  tableOut[2,1]<-AIC(base)
  tableOut[2,2]<-BIC(base)
  models <- list(base)
  
  ## estimate smoothing
  equation_update<-
  update.formula(baseForm,
                         paste(paste(". ~ .+s(variable, k=",smo,sep=""),")"))
  
  smoothing<-gam(equation_update, 
                 family=dist,
                 method="GCV.Cp",
                 data=data)
  tableOut[1,1]=AIC(smoothing)
  tableOut[1,2]=BIC(smoothing)
  models <- c(models, list(smoothing))
  
  ## estimate one model for each transformation
  namesLL=c("model_base","model_smoothing") ## model names
  for (j in 1:7){
    equation_update<-merge.formula(baseForm, formula(formsPack, rhs=j, lhs=0))
    model_c<-gam(equation_update, 
                 family=dist,
                 method="GCV.Cp",
                 data=data)
    tableOut[j+2,1]=AIC(model_c)
    tableOut[j+2,2]=BIC(model_c)
    models <- c(models, list(model_c))
    namesLL<-c(namesLL,paste("model",frameNames[j+2]
                             , sep="_"))
  }
  names(models)<-namesLL
  
  
  ## sort table after AIC
  tableOut[,3]=as.numeric(rank(tableOut[,1]))
  tableOut=tableOut[order(tableOut[,3]),]             
  
  ## sort models after AIC  
  order=paste("model",rownames(tableOut), sep="_")
  models<- models[order]
  
  
  ## formula for output
  formsAll=paste(as.character(Forms[2]),as.character(Forms[3]), sep="~")
  FormsAll=Formula(formula(paste(formsAll,as.character(formsPack[3]),sep="|")))
    
  #####################################################################
  ## return list  
  #####################################################################
  print(tableOut)   ## print table
  class(data)<-"data.frame"
  output=list("models"=models,
              "dataset"=data,
              "rank.table"=tableOut,                
              "fforms"=FormsAll) 
  class(output)<-"PJ"
  return(output)
}

plotff<-function(input){
   nModels= length(input$models) ## number of model
   namesLL<-names(input$models)
   ms<-grep("model_smoothing",namesLL)
   varNam= all.vars(formula(input$fforms, lhs=0, rhs=2))[1] ## variable to test
   contNam= all.vars(formula(input$fforms, lhs=0, rhs=1))## control variables
   depNam=all.vars(formula(input$fforms, lhs=1, rhs=0)) ## dependent variable
   input$dataset$variable=get(all.vars(formula(input$fforms, lhs=0, rhs=1)),input$dataset)  
   
   ######################################################################
   ### creating prediction frame 
   ######################################################################
   ## 100 points from min to max of variable
   ## median of all control variables
   names<-all.vars(input$fforms)
   pred_frame<-data.frame(matrix(NA, nrow=100, ncol=(length(names))))
   names(pred_frame)<-c(names)
   pred_frame$variable<-seq(as.numeric(quantile(input$dataset$x,0.05)),as.numeric(quantile(input$dataset$x,0.95)),(as.numeric(quantile(input$dataset$x,0.95))-as.numeric(quantile(input$dataset$x,0.05)))/100)[1:100]
   pred_frame[varNam]<-pred_frame$variable

   for (i in 1:length(contNam)){ ## predict y for each model
     pred_frame[contNam[i]]<-
     median(get(contNam[i],input$dataset))
  }
   
  ## predict dependent variable for each model
  for (i in namesLL){ ## predict y for each model
    pred_frame[paste("y",i,sep="_")]=predict(input$models[[i]],newdata=pred_frame, type="response")
  }   
  firstM=length(contNam) + 4 ## first model in pred_frame
  lastM=dim(pred_frame)[2] ## last model in pred_frame
  
  ### Plotting
  limx=c(min(pred_frame$variable),max(pred_frame$variable)) ## limits for x-axis
  limy=c(max(c(0,min(pred_frame[,firstM:lastM]))), max(pred_frame[,firstM:lastM])) ## limits for y-axis   
   
  #### plotting the the fit 
  par(mar=c(4, 4, 8.1, 12), xpd=TRUE)
  opt <- options("scipen" = 20)
   
  ## Start plot
  plot(y_model_smoothing~variable, data=pred_frame, type="l",main="",sub="",xlab="x",ylab=depNam, lwd=3, col="black", xlim=limx, ylim=limy) ## name variable of interest
  
  ## adding base and  functional form lines 
  color<-c("darkred","aquamarine4","darkgrey","blue3","orange3","brown3","chartreuse4","darkgreen","springgreen","gold","steelblue2","hotpink","blueviolet")
  color[ms]<-"black"
  n<-0
  colL<-vector()
   for (i in namesLL){
     n<-n+1
     lines(get(paste("y",i,sep="_"),pred_frame)~get("variable",pred_frame),
           col= color[n], lwd=2)
     colL=cbind(colL,as.character(color[n]))
   }
  namesL<-gsub("model_","",namesLL)
  legend(limx[2],limy[2],
         cex=1,bty="n",
         namesL, 
         fill=colL, 
         horiz=FALSE)
}
