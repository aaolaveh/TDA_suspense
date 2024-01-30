# ---
# jupyter:
#   jupytext:
#     formats: ipynb,R:light
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# ## Labeling MNI
#
# <ol>
#  <li> set_dictionaries_rois 
#  <li> set_dictionaries_attr
#  <li> find_structure 
#  <li> mni2cor
#  <li> my_round  

# +
library(rjson)

#------------------ set_dictionaries_rois --------------
#Set the dictionaries and the list of attributes for 
#function r_info
#Input: - file_info: file with json dictionaries
#       - file_nodes: file with json  attributes
#Output: - r_info: with those file sets 
#----------------- r_info -----------------
#Gives information of nodes of the given Atlas (Shen for me)
#Input: i: number of the region
#Output result: dictionarie with the MNI coordinates,
#      lobe, network and Brodmann area of region i

set_dictionaries_rois <- function(file_info,file_nodes){
    dict_info = fromJSON(file = file_info)
    dict_nodes = fromJSON(file = file_nodes)
    
    r_info <- function(i){
        info = dict_nodes$rois[[i]]
        
        l=toString(info$attr[1])
        n1=toString(info$attr[3])
        b=toString(info$attr[4])
        n2=toString(info$attr[5])
        
        result=list()
        result$'MNI'=c(info$x,info$y,info$z)
        result$'Lobe'= dict_info$'gui_Lobes'[[l]]
        result$'Network1'=dict_info$'gui_Networks'[[n1]]
        result$'Network2'=dict_info$'gui_Networks2'[[n2]]
        result$'Brodmann Area'=dict_info$'gui_BrodLabels'[[b]]
        
        return(result)
    }
    return(r_info)
}

# +
library(rjson)

#------------------ set_dictionaries_attr --------------
#Set the dictionaries and the list of attributes for 
#function attr_info
#Input: - file_info: file with json dictionaries
#       - file_nodes: file with json  attributes
#Output: - attr_info: with those file sets 
#----------------- attr_info -----------------
#Gives list of nodes for the given attributes 
#Input: i: list of attributes (Lobe, Network, Brodmann,Network2)
#Output result: list with the index of regions that posses the attributes

set_dictionaries_attr <- function(file_info,file_nodes){
    dict_info = fromJSON(file = file_info)
    dict_nodes = fromJSON(file = file_nodes)
    
    attr_info <- function(input_lista=list()){
        lista=list(Lobe=NULL, Network=NULL,Brodmann=NULL,Network2=NULL)
        lista=modifyList(lista,input_lista)
        n=length(lista)
        
        lista2=c('gui_Lobes','gui_Networks','gui_BrodLabels','gui_Networks2')
        
        wheres=list() #wheres store all numbers of the matches i.e Network=Def matches Default  in gui_Networks
                      #wheres[[3]]= NA NA NA '5' NA NA NA NA NA NA NA
        
        positions = c(1,3,4,5) #indices where information is stored in attributes
        for (j in 1:n){
            dict=dict_info[[lista2[j]]]
            if(is.null(lista[[j]])){
                y=rep(TRUE,length(dict))
            }else{
                x=lapply(dict, function(ch) grep(lista[[j]], ch)) #Find all matchings in dictionary
                y=x==1  #x to TRUE and FALSE
            }
            z=names(dict)[y] #keys of TRUES
            wheres[[j]]=z 
        }    
  
        #Find the indices for elements in wheres    
        indices=c()                 
        for (i in 1:length(dict_nodes$rois)){
            j=0
            repeat{
                j=j+1
                info = dict_nodes$rois[[i]]
                attr = toString(info$attr[positions[j]]) # OJO indices[j]
                z = wheres[[j]]
                if(sum(is.element(z,attr))==0){
                    break
                }
                if(j==n){ #If it meets all requirements it reachs n 
                    indices=c(indices,i)
                    break
                }
            }
        } 
        
        return(indices)
    }
    return(attr_info)
}
# -

#----------------------------  find_structure ---------------------------
#This function converts MNI coordinate to a description of brain structure in aal
#Input: - mni : the coordinates (MNI) of some points, in mm.  It is Mx3 matrix
#where each row is the coordinate for one point.
#        -DB (optional): The database. If is omit, make sure TDdatabase.mat is in the 
#        same folder
#Output: -one_line_result: A list of M elements, each describing each point.
#        -table_result:  A  MxN matrix being N the size of the database (DB)
#
find_structure<- function(mni){
    #Vectorize this functions
    
    mat=readMat('../Data/TDdatabase.mat') 
            
    if (is.vector(mni)){
        mni = matrix(mni,1,length(mni))
    }else{
        if(!is.matrix(mni)){stop("mni is not a matrix")}
    } 
    
    #round coordinates
    mni=apply(mni,c(1,2),round)
    
    T=matrix(c(2 ,0 ,0 ,-92,0,2,0,-128,
             0,0,2,-74,0,0,0,1),4,4,byrow=TRUE)
    
    index=mni2cor(mni,T)
    M=dim(index)[1] # Number of mni coordinates
    
    N=length(mat$DB) # Number of attributes
    list_result=list()
    one_line_result=rep("", M)
    
    for (i in 1:M){
        list_result[[i]]=list()
        for (j in 1:N){
            #mat$DB[[j]][[1]][[1]] is the j-th 3D-matrix 
            graylevel=mat$DB[[j]][[1]][[1]][index[i,1],index[i,2],index[i,3]] 
            if (graylevel == 0){
                 label = 'undefined'
            }else{    
                #mat$DB[[j]][[1]][[2]] is the list with regions
                label= mat$DB[[j]][[1]][[2]][[graylevel]][[1]][1,1]
            }
            #mat$DB[[j]][[1]][[3]] is the name of the attribute 
            name = mat$DB[[j]][[1]][[3]][1,1]  
            list_result[[i]][[name]]=label
            one_line_result[i] = paste(one_line_result[i], ' // ', label,sep="")
        }
    }
    return(list(one_line=one_line_result,lista=list_result))
}

#---------------------------- mni2cor --------------------------------
# convert mni coordinate to matrix coordinate
#Input: - mni : the coordinates (MNI) of some points, in mm.  It is Mx3 matrix
#        where each row is the coordinate for one point.
#        -T (optional): transform matrix coordinate is the returned coordinate in matrix.
#Output: -coords : Coordinate matrix
#
mni2cor <- function(mni,T=matrix(c(-4 ,0 ,0 ,84,0,4,0,-116,
                    0,0,4,-56,0,0,0,1),4,4,byrow=TRUE)){
    
    if (is.vector(mni)){
        mni = matrix(mni,1,length(mni))
    }else{
        if(!is.matrix(mni)){stop("mni is not a matrix")}
    } 
    
    if (dim(mni)[2] != 3){
        stop('are not 3-length coordinates')
        return(c())
    }
        
    a=cbind(mni,c(rep(1,dim(mni)[1])))
    b=t(solve(T))
    coords=a%*%b
    coords=coords[,1:3]
    
    if (is.vector(coords)){
        coords = matrix(coords,1,length(coords))
    }
        
    coords = apply(coords,c(1,2),my_round)
    
    return(coords)
}    

#----------------- my_round -------------------
#Integer rounding that behaves like round from MATLAB
my_round <- function(x){
    r=x-floor(x)
    if(r == 0.5){
        if (x<0){ return(x-0.5)}
        else{return(x+0.5)}
    }else{
        return(round(x))
    }
}

# ## Mapper
#
# <ol>
#   <li>cart2pol
#   <li>color_bar
#   <li>id 
#   <li>make_image
#   <li>make_integer_fun
#   <li>make_integer_fun2
#   <li>make_layout
#   <li>make_pie_fun 
#   <li>multiplot
#   <li>my_multi_with_colorbars 
#   <li>my_names
#   <li>my_plot_params
#   <li>plot_mapper_pie 
#   <li>vertex_size 

#=================================================================
#cart2pol
# Translate cartesian coordinates into polar coordinates
# based on coordinates x and y
# degrees indiccates if the result is in degrees (T) or radians (F)
cart2pol<-function(x,y,degrees=TRUE){
    # calculate r with sqrt of x and y squared
    r = sqrt(x^2 + y^2)
    # calculate theta with arctan
    theta =atan2(y, x)
    
    ## adjust angle for appropriate quadrant
    ## quadrants I and II need no adjustment
    ## quadrant III and IV, add 360
    theta = theta + (y < 0)*2*pi
    
    # return as degrees if requested
    if(degrees)
    {
       theta = theta*180/pi
    }
    
    result=matrix(c(r,theta),ncol=2)
    colnames(result) = c("r","theta")
    return(result)
}

#=================================================================
#color_bar
#Make a colorbar
color_bar<-function(zlim,breaks=NULL,lab.breaks=NULL,col=NULL,axis.args=NULL,
                    width=.2,heigth=.75,pos=NULL,title=NULL,
                    title.cex=NULL,title.line=NULL,title.side=NULL,title.at=NULL,
                    horizontal=FALSE){
    # set defaults for color scale 
    if( is.null(col))  {
        col=  tim.colors(nlevel)}
    else{
        nlevel = length(col)
    }
    
    #  Now set up the breaks
    if( is.null(breaks)){
        midpoints<- seq( zlim[1], zlim[2],,nlevel)
        delta<- (midpoints[2]- midpoints[1])/2
        # nlevel +1 breaks with the min and max as midpoints 
        # of the first and last bins.
    
        breaks <- c( midpoints[1]- delta, midpoints + delta)
        breaks=round(breaks,5)
    }  
    
    old.par <- par(no.readonly = TRUE)
    
    #plotting region of the colorbar
    plt=old.par$plt
    if(horizontal){
        if(is.null(pos)){
            pos=c(.2,.1)
        }
        plt[3]=pos[1]
        plt[4]=plt[3]+width
        plt[1]=pos[2]
        plt[2]=plt[1]+heigth
    }else{
        if(is.null(pos)){
            pos=c(.1,.1)
        }
        plt[1]=pos[1]
        plt[2]=plt[1]+width
        plt[3]=pos[2]
        plt[4]=plt[3]+heigth
    }
    
    ix = 1:2
    iy = breaks
    nBreaks = length(breaks)
    midpoints<- (breaks[1:(nBreaks-1)] +  breaks[2:nBreaks] )/2
    iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints)) 
     
    pin=old.par$pin
           
    if (!horizontal) {
        pin[1]=.2
        par(mar=c(.1,2,.1,2),pin=pin,pty='m',plt=plt)
        image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
              ylab = "", col = col, breaks=breaks)
    }else {
        pin[2]=.2
        par(mar=c(2,.1,2,.1),pin=pin,pty='m',plt=plt)
        image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
              ylab = "", col = col, breaks=breaks)
    }
     
    if (!is.null(lab.breaks)) {
        # axis with labels at break points
        axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2), 
            at = breaks, labels = lab.breaks,
            cex.axis=ifelse(horizontal,5*par()$pin[2],5*par()$pin[1])), axis.args)
        # now add the axis to the legend strip.
        # notice how all the information is in the list axis.args
        do.call("axis", axis.args)
    }
    else {
        if ((!is.null(axis.args))){
            # If lab.breaks is not specified ( with or without breaks), pretty
            # tick mark locations and labels are computed internally,
            # or as specified in axis.args at the function call
            axis.args <- c(list(side = ifelse(horizontal, 3, 4), 
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)), 
                        cex.axis=ifelse(horizontal,6*par()$pin[2],6*par()$pin[1]), axis.args)
            do.call("axis", axis.args)
        } 
    }  
  
    # add a title to the axis if information has been  supplied
    # using the mtext function. The arguments to mtext are
    # passed as a list like the drill for axis (see above)
    #
    if (!is.null(title)) {
        if(is.null(title.line)){title.line=ifelse((!is.null(axis.args)),2.3,0.5)}
        if(is.null(title.side)){title.side=ifelse(horizontal,3,4)}
        if(is.null(title.cex)){title.cex=ifelse(horizontal,5*par()$pin[2],5*par()$pin[1])}
        if(is.null(title.at)){title.at=(zlim[2]+zlim[1])/2}
        legend.args <- list(text = title, cex=title.cex, side = title.side, 
                            at=title.at,line = title.line)
        # just guessing at a good default for line argument!
    }
    else{legend.args=NULL}
    
    # add the label using mtext function
    if (!is.null(legend.args)) {
        do.call(mtext, legend.args)
    }
     # add a box around legend strip
    box()
}


#==================================================================
id<-function(x){x}
#===================================================================
#make_image.
#Input: -a real vector vec
#Outout: A function g where g(i)=v[i]
make_image<-function(vec){
    g<-function(x){vec[x]}
    return(g)
}


#====================================================================
#make_integer_fun.
#Input: -f is a real function 
#       -dom is the domain (vector) where f is going to be applied
#       -min and max of f(dom)
#       -n_int is the number of intervals f(dom) is divided
# Output: -g a vectorixed function where g(x) is the number of the interval 
#         where f(x) belongs
make_integer_fun<-function(f,dom=NULL,min_max=NULL,n_int){
    if(missing(dom) || is.null(dom)){
        if(missing(min_max) || is.null(min_max)){
        stop("Missing domain or maximum and minimum")
        }  
    }else{ #Calculate minimum and maximum from dom
        dom=as.array(dom)
        v_f=Vectorize(f)
        image=v_f(dom)
        min_max=c(1,2)
        min_max[1]=min(image)
        min_max[2]=max(image)
    }
    #length of each interval 
    l=(min_max[2]-min_max[1])/(n_int)
    #Calculates the interval where f(x) belongs
    g<-function(x){
        if(f(x)==min_max[1]){return(1)}
        else{ceiling((f(x)-min_max[1])/l)}
    }
    return(Vectorize(g))
}


#====================================================================
#make_integer_fun2.
#Input: -f is a real function 
#       -inter: vector with limits of intervals of the image of f
# Output: -g a vectorixed function where the integer part g(x) is the number of 
#         the interval where f(x) belongs and decimal part an approximation of highness
make_integer_fun2<-function(f,inter,continuous=FALSE){
    h<-function(x,min,max){return((x-min)/(max-min))}
    n_int=length(inter)-1
    print(paste0("Calculating on " ,toString(n_int), " intervals"))
    #Calculates the interval where f(x) belongs
    g<-function(x){
        int=findInterval(f(x),inter,rightmost.closed=TRUE)
        if (int == 0 | int == length(inter)){
            stop("Value outside limits")
        }else{
            if (continuous){
                value=h(f(x),inter[int],inter[int+1])+int #linear value inside interval
                if (value==length(inter)){
                    return(value-0.00001)
                }else{
                    return(value)  
                }
            }else{return(int)}
        }
    }
    return(Vectorize(g))
}

#========================================================
#make_layout.
#Intercalated layout for my mapper output
#Input -total: Total of mappers
#      -cols: number of columns
#      -Intercalado
#Output: Layout matrix
make_layout<-function(total,cols,intercalado=FALSE){
    ntotal=total/2
    rows=ntotal/cols
    if (intercalado) {
        vector=c()
        for (i in 1:rows){
            vector=c(vector,(1+(i-1)*cols):(i*cols),(1+(i-1)*cols+ntotal):(i*cols+ntotal))
        }
    }else{
       vector=1:total
    }
    return(matrix(vector,ncol=cols,byrow=TRUE))
}

#=====================================================================
#make_pie_function. 
#Input: -f is a real function
#       -vert a list which is element is a vector. 
#Output: values is a list where the element u is a vector and u[i] stores the
#        number of entries with value f_i (the i-th posible value of f)  
make_pie_fun<-function(f,vert,n_int){
    #Calculates integer function from f and vert
    g=make_integer_fun(f,unique(unlist(vert)),n_int=n_int)
    #Make the pie function 
    make_pie<-function(vert){
        g_vert=lapply(vert,g) #Image of vert using g
        types=1:n_int
        values=list()
        i=1
        for (ver in g_vert){ 
            a=sapply(types,function(x)sum(ver==x)) #Total on each interval
            values[[i]]=a
            i=i+1
        }
        return(values) 
    }    
    return(make_pie)
}

#=====================================================================
# Multiple plot function
#
# - plots is a list of parameters for my_plot function
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(fun_plot, plots=NULL, title, cols=1, layout=NULL){
  library(grid)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots<2) {
    do.call(fun_plot,plots[[1]])
  }else {
    layout(layout)
    for (plt in plots){
        do.call(fun_plot,plt)
    }
    mtext(title,side=3, line=16, at=-4, cex=.8 , font=2)
  }  
} 


# +
#==================================================================================
#my_multi_with_colorbars.
#Input: -fun_plot: A plotting function 
#       -plots: Arguments used by fun_plot
#         -layout: A matrix specifying the layout. 
#         -colors: List of colors for bars
#         -bars_per_line: Numbers of bars in same line
#         -title: Title of plot
#         -title_bars: List of title for the bars
#         -ticks: Labels and positions in the bars
#         -mins,maxs: List of mins and maxs for the bars   
#         -Horizontal: TRUE plots horizontal bars, else, vertical
#         -width, heigth: Dimensions of bars
#Output: Multiplot function with colorbars 

my_multi_with_colorbars<-function(fun_plot,plots=NULL,layout,colors,bars_per_line=1,
                                  title=NULL,title_bars,axis.args=NULL,mins,maxs,
                                  horizontal=FALSE,width=.2,heigth=.7){
  numPlots = length(plots)
  
  rows=nrow(layout)
  cols=ncol(layout) 
  #cocient and reminder of rows and cols divided by bars_per_line  
  c1=ncol(layout) %/% bars_per_line
  r1=ncol(layout) %% bars_per_line
  c2=nrow(layout) %/% bars_per_line
  r2=nrow(layout) %% bars_per_line
  
  #reminders must be zero
  if(horizontal){
      if(r1 != 0){stop("bars_per_line must divide the number of columns")}
  }else
      if(r2 != 0){stop("bars_per_line must divide the number of rows")}
    
  total=length(colors)
  #lines with bars 
  lines=ceiling(total/bars_per_line)
  m=max(layout)
  #Computes a new layout to put the bars  
  if(horizontal){
      N=matrix(0:(lines*ncol(layout)-1),nrow=lines,byrow=TRUE)
      N=apply(N,c(1,2),function(x){x%/%c1+1}) #N is te matrix to be added to layout
      N=N+m
      layout=rbind(layout,N)
      if(!is.null(title)){
          N=matrix(1,ncol=ncol(layout))
          layout=rbind(N,layout+1) #N is te matrix to be added to layout
          layout(layout,heights=c(.2,rep(.5,rows),rep(.2,lines))) #Dimensions of plots in layout
          par(mar=c(.1,.1,.1,.1),pin=c(5,.5))
          plot.new()
          text(0.5,0.5,title,cex=1.3,font=2)
      }else{
          layout(layout,heights=c(rep(.5,rows),rep(.2,lines))) #Dimensions of plots in layout
      }
  }else{
      N=matrix(0:(lines*nrow(layout)-1),ncol=lines)
      N=apply(N,c(1,2),function(x){x%/%c2+1}) #N is te matrix to be added to layout
      N=N+m
      layout=cbind(layout,N)
      if(!is.null(title)){
          N=matrix(1,ncol=ncol(layout))
          layout=rbind(N,layout+1)#N is te matrix to be added to layout
          layout(layout,widths=c(rep(.5,cols),rep(.2,lines)),heights=c(.2,rep(.5,rows))) #Dimensions of plots in layout
          par(mar=c(.1,.1,.1,.1),pin=c(5,.5))
          plot.new()
          text(0.5,0.5,title,cex=1.3,font=2)
      }else{
          layout(layout,widths=c(rep(.5,cols),rep(.2,lines))) #Dimensions of plots in layout
      }
  }
  
  #plots using fun_plot  
  for (plt in plots){
      do.call(fun_plot,plt)
  }
  cat('Ploting graphs done \t')
  
  #plots colorbars
  for (i in 1:total){
      if(!is.null(axis.args[[i]])){#Ticks converted to labels
        args=axis.args[[i]]
      }else{args=NULL}
      
      col=na.omit(unlist(colors[[i]]))
      
      #Plotting color_bar
      color_bar(zlim=c(mins[[i]],maxs[[i]]),axis.args=args,
                width=width,heigth=heigth,col=col,title=title_bars[[i]],
                title.side=NULL,title.line=NULL,title.cex=.6,horizontal=horizontal)
  }
  cat('Plotting colorbars done \t')
}
# -
#=====================================================================         
#my_names is.Input:file. Output is a list
#with 7 elements:  part, metric, filter, dist_file, file, pdf and title 
my_names<-function(file){
    lista=list()
    file_sp=strsplit(file , '/',fixed=TRUE)[[1]]
    path=paste0(file_sp[1:(length(file_sp)-2)],collapse='/')
    name=file_sp[length(file_sp)]
    name_sp=strsplit(name , '_',fixed=TRUE)[[1]] #split string 
    lista[["part"]]=name_sp[1]
    lista[["metric"]]=name_sp[2]
    lista[["filter"]]=str_remove(name_sp[3],".npy")
    name_file=paste0(name_sp[1:2],collapse='_')
    lista[["dist_file"]]=paste(path,'/dist/',name_file,'.npy',sep="")
    
    name=str_replace(name_file,"../Data/dist","mapper")
    name=paste(name,"_",lista[["filter"]],sep="")
    
    lista[["file"]]=paste("../Data/mapper/",name,".rds",sep="")
    
    lista[["pdf"]]=paste("../Imagenes/mapper_graphs/mapper_",name,".pdf",sep="")
    
    lista[["title"]]=paste("Mapper for part = ",lista[["part"]],", metric=", lista[["metric"]],", filter =", lista[["filter"]])
    return(lista)
}

#=====================================================================         
#my_names2 is.Input:file. Output is a list
#with 5 elements:   index, metric, filter, dist_file and file
my_names2<-function(file){
    lista=list()
    file_sp=strsplit(file , '/',fixed=TRUE)[[1]]
    path=paste0(file_sp[1:(length(file_sp)-2)],collapse='/')
    name=file_sp[length(file_sp)]
    name_sp=strsplit(name , '_',fixed=TRUE)[[1]] #split string 
    lista[["index"]]=name_sp[1]
    lista[["metric"]]=name_sp[2]
    lista[["filter"]]=str_remove(name_sp[3],".npy")
    name_file=paste0(name_sp[1:2],collapse='_')
    lista[["dist_file"]]=paste(path,'/dists/',name_file,'.npy',sep="")
    
    name=str_replace(name_file,"Data/dists","mapper")
    name=paste(name,"_",lista[["filter"]],sep="")
    
    lista[["file"]]=paste("../../../../Documentos/Mapper/Data/mapper/",name,".rds",sep="")
    lista[["title"]]=paste("Mapper for index = ",lista[["index"]],", metric=", lista[["metric"]],", filter =", lista[["filter"]])

    return(lista)
}

#=====================================================================
#params to be used for plot_mapper_pie
my_plot_params<-function(m,layout,pie_function,palette,pars=NULL,labels=NULL,margin=c(1,1,1,1),pin=c(1,1)){
    vert=m$vertices
    params=list() #Parametes for plotting
    params[["graph"]]=m$as_igraph()
    params[["layout"]]=layout
    params[["pies"]]=pie_function(vert) #Computing pies
    params[["size"]]=vertex_size(vert)
    params[["colors"]]=palette
    if (is.null(pars)){
        params[["title"]]=NULL     
    }else{
        params[["title"]]=paste("No. of Intervals = ", pars[1], 
        "| % of overlap = ",pars[2], "| No. of vertices = ", length(vert))
    }
    if (is.null(labels)){
        params[['labels']]=NA
    }else{
        params[['labels']]=labels
    }
    params[["margin"]]=margin
    params[["pin"]]=pin  
    return(params)
} 

#======================================================================
#plot mapper graph with vertex shape of pies 
plot_mapper_pie<-function(graph,layout,pies,size,colors,title=NULL,labels=NULL,pin,margin=c(1,1,1,1)){
    par(mar=margin,pin=pin)
    plot.igraph(graph,layout=layout,vertex.label=labels,   #NA;sapply(pies,sum), 
                vertex.label.color="black",
                vertex.shape="pie",vertex.pie = pies, 
                vertex.pie.color = colors,vertex.size=size,
                vertex.pie.lty = 0, vertex.pie.density=-1, 
                edge.color="black",edge.width=.5,pin=pin)
    if(!is.null(title)){
        mtext(title,side=3,cex=.3*pin[1],font=1,line=0)
    }
}

#=====================================================================
#vertex_size gives the size of vertex
# It assings size=10 to the smallest vertex
# and 25 to the biggest. Anything else is proportional                
vertex_size<-function(vert,m=10,M=25){
    lens=sapply(vert,length)
    mini=min(lens)
    maxi=max(lens)
    g<-function(x){
       r=(M-m)*(x-mini)/(maxi-mini)+m
       return(r)
    }          
    sizes=sapply(lens,g)
    return(sizes)
}

# ## Graph
#
# <ol>
#   <li>bapply
#   <li>connectors
#   <li>recog_block
#   <li>rotate
#   <li>tcm_matrix


#================================================================
#bapply applies a matrix function on blocks of the matrix M
#Input:   -FUN: Function to be applied
#         -M: Matrix
#         -X,Y: List or matrix of intervals (n x 2). A element 
#         of the list or row would be c(1,5) representing the interval between indices 1 and 5
#         -comb: Use all combinations of intervals in X and Y, Otherwise the intervals will be matched
#         one to one
# Output: List of results of FUN applied to defined blocks
bapply<-function(FUN,M,X,Y=NULL,comb=FALSE){
    block_from_index <-function(A,B,index){
        j1=index[1]
        j2=index[2]
        int1=A[j1,]
        int2=B[j2,]
        block=M[int1[1]:int1[2],int2[1]:int2[2]]
        return(block)
    }
    
    if(is.list(X)){
        X=matrix(unlist(X,use.names=FALSE), nrow = length(X), byrow = TRUE)
    }
    if(ncol(X)!=2){stop("X must have 2 columns")}
    
    if(!is.null(Y)){
        if(is.list(Y)){
            Y=matrix(unlist(Y,use.names=FALSE), nrow = length(Y), byrow = TRUE)
        }
        if(ncol(Y)!=2){stop("Y must have 2 columns")}
        rY=nrow(Y)
    }
    
    rX=nrow(X)
    if (comb){ #comb is TRUE
        if (is.null(Y)){
            index=combinations(rX,2,repeats.allowed=TRUE)
            Y=X
        }
        else{
            index=expand.grid(1:rX,1:rY)
        }
    }else{
        index=matrix(1:rX,nrow=rX,ncol=2)
        if  (!is.null(Y)){
            if (rX != rY){stop("X and Y must have same dimensions")}
        }else{
            Y=X
        }
    }
    result=list()
    for (i in 1:nrow(index)){
        block=block_from_index(X,Y,index[i,])
        result[[i]]=FUN(block)  
    }
    return(result)   
}

#===============================================================
#connectors returns the connectors of the graph
#connector is a node when remove change the number of connected components of the graph
#Input - graph: igraph element
#Outpout - list of connector nodes
connectors <- function(graph){
    comps=components(graph)$no
    v=vcount(graph)
    conns=c()
    for (i in 1:v){
        graph1=delete_vertices(graph,i)
        new_comps=components(graph1)$no
        if (new_comps > comps){conns=c(conns,i)}
    }
    return(conns)
}

#====================================================================
#recog_block returns the continous intervals on the vector
#Input: -vec: vector 
#       -output: matrix where each row has the firts ad last element
#       of the interval
recog_block <-function(vec){
    ord_vec=sort(vec)
    l=length(ord_vec)
    blocks=matrix(0,nrow=0,ncol=2)
    vec=c(0,0)
    vec[1]=ord_vec[1]
    for (j in 2:l){
        a=ord_vec[j-1]
        b=ord_vec[j]
        if (a+1<b){
            vec[2]=a
            blocks=rbind(blocks,vec)
            vec[1]=b
        }
    }
    vec[2]=ord_vec[l]
    blocks=rbind(blocks,vec)
    return(blocks)
}

rotate <- function(coords,center,angle){
    coords = sweep(coords, 2, center)
    M <- matrix( c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2 )
    coords = coords  %*%  M 
    coords = sweep(coords, 2, -center)
    return(coords)
}

#=====================================================================
#tcm_matrix returns temporal connectivity matrix from graph 
#(vertices are composed, if not it is the same as adjacency matrix)
#Input: -vert: list of vertices of the graph
#       -adj: adjacency list of graph
#Output: -M: temporal connectivity matrix
tcm_matrix <-function(vert,adj){
    v=length(vert) #total of vertices
    p=length(unique(unlist(vert,use.names = FALSE))) #total of points
    M=matrix(0,nrow=p,ncol=p) #initialize matrix
    for (ver in names(vert)){
        points=vert[[ver]] #points inside vertex
        v_neig=c(adj[[ver]],ver) #neighbor vertices of vertex + itself
        for (vn in v_neig){ # for each neighbor
            ver=toString(vn) 
            neig=vert[[ver]] #recover points inside neighbor vertex
            for (j in points){   
                for (k in neig){
                    if (j != k){M[j,k]=M[j,k]+1}  #+1 when nodes are in neighbor vertices}
                }
            }
        }
        
    }
    return(M)         
}

