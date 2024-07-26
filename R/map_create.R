# R file for functions relating to the creation of maps visually showing region-specific parameters
#-------------------------------------------------------------------------------
#' @title map_shapes_load
#'
#' @description Create a set of shape data to make into one or more maps
#'
#' @details Takes in one or more shapefiles (.shp) and extracts data for selected regions of a specified type
#'
#' @param regions Vector of names of the regions for which to extract data
#' @param shapefiles Vector of names of shapefiles from which to extract data
#' @param region_label_type Type of region ID used in vector of regions, corresponding to a data type appearing in the
#'   shapefiles (e.g. "GID_1" for first subnational region IDs in the form of the three-letter country code plus a
#'   number, e.g. "AGO.1_1", "AGO.2_1", etc.)
#' '
#' @export
#'
map_shapes_load <- function(regions=c(),shapefiles=c(),region_label_type=""){

  assert_that(is.character(regions))
  assert_that(is.character(shapefiles))
  assert_that(is.character(region_label_type))

  n_regions=length(regions)

  for(i in 1:length(shapefiles)){
    shape_data=read_sf(shapefiles[i])
    assert_that(region_label_type %in% names(shape_data))
    if(i==1){
      shape_data_all=list(regions=regions,shapes=shape_data$geometry[1],
                          lat_min=Inf,long_min=Inf,lat_max=-Inf,long_max=-Inf)
    }
    bbox=st_bbox(shape_data)
    shape_data_all$lat_min=min(shape_data_all$lat_min,bbox[2])
    shape_data_all$lat_max=max(shape_data_all$lat_max,bbox[4])
    shape_data_all$long_min=min(shape_data_all$long_min,bbox[1])
    shape_data_all$long_max=max(shape_data_all$long_max,bbox[3])
    j=match(region_label_type,names(shape_data))
    file_regions=shape_data[[j]]
    for(n_region in 1:n_regions){
      k=match(regions[n_region],file_regions)
      if(is.na(k)==FALSE){
        shape_data_all$shapes[[n_region]]=shape_data$geometry[[k]]
      }
    }
  }
  assert_that(length(shape_data_all$shapes)==n_regions,msg="Region data missing")
  for(n_region in 1:n_regions){assert_that(is.null(shape_data_all$shapes[[n_region]])==FALSE,
                                                    msg=paste("No shape data found for region",n_region))}

  return(shape_data_all)
}
#-------------------------------------------------------------------------------
#' @title create_map
#'
#' @description Create a map of one or more regions with colours denoting parameter values
#'
#' @details Takes in region shape data generated using map_shapes_load() and parameter values for each region, plots map
#'   of the regions and fills regions with colour based on parameter values
#'
#' @param shape_data Region shape data generated using map_shapes_load()
#' @param param_values Vector of parameter values for regions in shape_data
#' @param scale Vector of scale intervals to use for param_values
#' @param colour_scale Vector of colours with size greater than or equal to scale - used to convert scale to colours
#' @param pixels_max Number of pixels to use for largest dimension of map
#' @param text_size Size of text to appear in legend and titles
#' @param display_axes TRUE/FALSE flag indicating whether to frame map and display latitude/longitude axes
#' @param border_colour_regions Colour to use for borders of regions. Set to NA if borders to be invisible.
#' @param ... Additional optional parameters: \cr
#'    lat_max,lat_min,long_max,long_min: borders if not to be set default \cr
#'    additional_border_shapes: Shape data for optional additional borders \cr
#'    border_colour_additional: colour to use for additional borders if any. \cr
#     map_title: Title to show above map \cr
#     legend_title: Title to show above legend \cr
#     legend_position: Position to place map legend if to be used \cr
#     legend_format: Number format to use for scale values in legend if used \cr
#     legend_dp: Number of decimal places to use in scale values in legend \cr
#     legend_columns: Number of columns in which to display legend values \cr
#     output_file: Name of file to which to output map \cr
#' '
#' @export
#'
create_map <- function(shape_data=list(),param_values=c(),scale=c(),colour_scale=c(),pixels_max=720,text_size=1,
                       display_axes=FALSE,border_colour_regions="grey",...){

  assert_that(is.list(shape_data))
  assert_that(is.numeric(param_values))
  assert_that(is.numeric(scale))
  assert_that(is.logical(display_axes))
  n_regions=length(param_values)
  assert_that(n_regions==length(shape_data$shapes))
  ap=list(...) #Get additional optional parameters

  #Set map dimensions
  if(is.null(ap$lat_max)){ap$lat_max=shape_data$lat_max}
  if(is.null(ap$lat_min)){ap$lat_min=shape_data$lat_min}
  if(is.null(ap$long_max)){ap$long_max=shape_data$long_max}
  if(is.null(ap$long_min)){ap$long_min=shape_data$long_min}
  height_ll=ap$lat_max-ap$lat_min
  width_ll=ap$long_max-ap$long_min
  pixel_scale=pixels_max/max(height_ll,width_ll)
  width_px=width_ll*pixel_scale
  height_px=height_ll*pixel_scale

  #Assign parameter values within scale
  assert_that(min(param_values,na.rm=TRUE)>=min(scale))
  assert_that(max(param_values,na.rm=TRUE)<=max(scale))
  scale_values=rep(NA,length(param_values))
  for(i in 1:length(param_values)){
    scale_values[i]=findInterval(param_values[i],scale)
  }
  if(is.null(ap$legend_format)==FALSE && ap$legend_format=="integer"){n_intervals=length(scale)} else {n_intervals=length(scale)-1}

  #Create legend labels
  if(is.null(ap$legend_position)==FALSE){
    assert_that(ap$legend_position %in% c("bottomright","bottom","bottomleft","left","topleft","top","topright","right","center"))
    assert_that(is.null(ap$legend_format)==FALSE)
    assert_that(ap$legend_format %in% c("f","e","pc","integer"))
    if(ap$legend_format=="integer"){assert_that(is.integer(param_values) && is.integer(scale))}
    legend_labels=rep("",n_intervals)
    if(ap$legend_format=="integer"){
      for(i in 1:n_intervals){
        legend_labels[i]=paste0(scale[i])
      }
    }
    if(ap$legend_format=="pc"){
      for(i in 1:n_intervals){legend_labels[i]=paste0(formatC(scale[i]*100,format="f",digits=ap$legend_dp)," - ",
                                                     formatC(scale[i+1]*100,format="f",digits=ap$legend_dp))}
    }
    if(ap$legend_format=="f"){
      for(i in 1:n_intervals){legend_labels[i]=paste0(formatC(scale[i],format="f",digits=ap$legend_dp)," - ",
                                                     formatC(scale[i+1],format="f",digits=ap$legend_dp))}
    }
    if(ap$legend_format=="e"){
      for(i in 1:n_intervals){legend_labels[i]=paste0(formatC(scale[i],format="e",digits=ap$legend_dp)," - ",
                                                     formatC(scale[i+1],format="e",digits=ap$legend_dp))}
    }
  }

  #Set colours
  ratio=length(colour_scale)/n_intervals
  values=ratio*c(1:length(colour_scale))[c(1:n_intervals)]
  for(i in 1:n_intervals){values[i]=max(1,floor(values[i]))}
  colour_scale2 <- colour_scale[values]

  #Create graph
  par(mar=c(1,1,1,1))
  if(is.null(ap$output_file)==FALSE){png(filename=ap$output_file,width=width_px,height=height_px)}

  matplot(x=c(ap$long_min,ap$long_max),y=c(ap$lat_min,ap$lat_max),col=0,xlab="",ylab="",
          axes=display_axes,frame.plot=display_axes)
  for(n_region in 1:n_regions){
    plot(st_geometry(shape_data$shapes[[n_region]]),col=colour_scale2[scale_values[n_region]],border=border_colour_regions,
         add=TRUE)
  }
  if(is.null(ap$additional_border_shapes)==FALSE){
    n_shapes_additional=length(ap$additional_border_shapes$shapes)
    for(j in 1:n_shapes_additional){
      plot(st_geometry(ap$additional_border_shapes$shapes[[j]]),col=NA,border=ap$border_colour_additional,add=TRUE)
    }
  }
  if(is.null(ap$legend_position)==FALSE){
    if(is.null(ap$legend_columns)){ap$legend_columns=1}
    legend(ap$legend_position,legend=legend_labels,fill=colour_scale2,cex=text_size,title=ap$legend_title,
           ncol=ap$legend_columns)
  }
  title(main=ap$map_title,cex.main=text_size)

  if(is.null(ap$output_file)==FALSE){dev.off()}
  par(mar=c(4,4,4,4))

  return(NULL)
}
