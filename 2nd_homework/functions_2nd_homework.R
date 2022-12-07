#first of all set as a working directory the folder where there are yours transcr. 
data
MvAplot <- function(exprData, pdffilename){
  #exprData <- read.table(exprData)
  exprData[exprData==0] <- 1 #replace all 0 to 1
  destination = pdffilename
  pdf(file=destination) #open PDF
  for (i in (2:ncol(exprData))) {
    M <- log(exprData[,1], base = exp(2)) - log(exprData[,i], base = exp(2))
    A <- (log(exprData[,1], base = exp(2)) + log(exprData[,i], base = exp(2)))/2
    plot(A, M, main = paste("sample 1 vs sample ",i)) #assign names to the differents plots
  }
  dev.off() #turn off pdf plotting
  
}

MvAplot_norm <- function(exprData, Atrim = c(0,8), Mtrim=0.02, pdffilename){
  exprData <- read.table(exprData)
  exprData[exprData==0] <- 1 #replace all 0 to 1
  destination = pdffilename
  pdf(file=destination) #open PDF
  #scale the data by their sequencing depth and multiply by 10^6
  for(j in (1:ncol(exprData))){
    exprData[,j] = (exprData[,j]/sum(exprData[,j]))*10^6
  }
  for (i in (2:ncol(exprData))) {
    #create a list where to save the selected M from selected A values
    Msel_tot <- c()
    #M and A will have the same length of the number of rows of the starting exprData matrix
    M <- log(exprData[,1], base = exp(2)) - log(exprData[,i], base = exp(2)) 
    A <- (log(exprData[,1], base = exp(2)) + log(exprData[,i], base = exp(2)))/2
    for (el in (1:length(A))){
      if(A[el]>Atrim[1] & A[el]<Atrim[2]){
        Msel <- log(exprData[el,1], base = exp(2)) - log(exprData[el,i], base = exp(2))
        Msel_tot <- c(Msel_tot, Msel)
      }
    }
    #compute the scaling factor
    sf = mean(Msel_tot, trim = Mtrim)
    M_normalized = M + sf
    plot(A, M_normalized, main = paste("sample 1 vs sample ",i)) #assign names to the different plots
  }
    dev.off() #turn off pdf plotting
}






TMMnorm <- function(exprData, annot, Atrim = c(0,8), Mtrim=0.02){
  #in that function I've not already considered the annotation data set
  exprData <- read.table(exprData)
  exprData[exprData==0] <- 1 #replace all 0 to 1
  #create the normalized matrix with the same shape of exprData matrix 
  normalized_matrix <- matrix(nrow = nrow(exprData), ncol = ncol(exprData))
  vector_of_sf <- c()
  
  #scale the data by their sequencing depth and multiply by 10^6
  for(j in (1:ncol(exprData))){
    exprData[,j] = (exprData[,j]/sum(exprData[,j]))*10^6
  }
  
  #define the value of reference sample(first column) of the normalized_matrix 
  for (i in (2:ncol(exprData))) {
    #create a list where to save the selected M from selected A values
    Msel_tot <- c()
    #M and A will have the same length of the number of rows of the starting exprData matrix
    M <- log(exprData[,1], base = exp(2)) - log(exprData[,i], base = exp(2)) 
    A <- (log(exprData[,1], base = exp(2)) + log(exprData[,i], base = exp(2)))/2
    for (el in (1:length(A))){
      if(A[el]>Atrim[1] & A[el]<Atrim[2]){
        Msel <- log(exprData[el,1], base = exp(2)) - log(exprData[el,i], base = exp(2))
        Msel_tot <- c(Msel_tot, Msel)
      }
    }
    #compute the sf and store it in the vector
    sf = mean(Msel_tot, trim = Mtrim)
    vector_of_sf <- c(vector_of_sf, sf)
    
    #this last part is not finished
    #scale the genes in the original scale by their length(not already implemented) and multiply by 10^3
    normalized_matrix[,i] = M + sqrt(2^sf)
  }
  return(normalized_matrix)
}


