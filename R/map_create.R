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
  shape_data_all=list(regions=regions,shapes=rep(NULL,n_regions),
                      lat_min=Inf,long_min=Inf,lat_max=-Inf,long_max=-Inf)

  for(i in 1:length(shapefiles)){
    shape_data=read_sf(shapefiles[i])
    assert_that(region_label_type %in% names(shape_data))
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
#' @param lat_max Maximum latitude (set to NULL to default to maximum latitude from shape_data)
#' @param lat_min Minimum latitude (set to NULL to default to minimum latitude from shape_data)
#' @param long_max Minimum longitude (set to NULL to default to minimum longitude from shape_data)
#' @param long_min Minimum longitude (set to NULL to default to minimum longitude from shape_data)
#' @param display_axes TRUE/FALSE flag indicating whether to frame map and display latitude/longitude axes
#' @param map_title Title to show above map
#' @param legend_title Title to show above legend
#' @param legend_position Position to place map legend (select from"bottomright", "bottom", "bottomleft", "left",
#'   "topleft", "top", "topright", "right" and "center")
#' @param legend_format Number format to use for scale values in legend - set to "f" (basic number), "e" (scientific
#'   notation) or "pc" (percentage)
#' @param legend_dp Number of decimal places to use in scale values in legend (e.g. if legend_format is "f" and
#'   legend_dp is 1, numbers will appear as e.g. 10.0)
#' @param legend_columns Number of columns in which to display legend values
#' @param output_file Name of file to which to output map; if set to NULL, map is displayed without being saved to file
#' '
#' @export
#'
create_map <- function(shape_data=list(),param_values=c(),scale=c(),colour_scale=c(),pixels_max=720,text_size=1,
                       lat_max=NULL,lat_min=NULL,long_max=NULL,long_min=NULL,display_axes=FALSE,map_title="",
                       legend_title="",legend_position="topleft",legend_format="f",legend_dp=1,legend_columns=1,
                       output_file=NULL){

  assert_that(is.list(shape_data))
  assert_that(is.numeric(param_values))
  assert_that(is.numeric(scale))
  assert_that(is.logical(display_axes))
  assert_that(legend_position %in% c("bottomright","bottom","bottomleft","left",
                                     "topleft","top","topright","right","center"))
  assert_that(legend_format %in% c("f","e","dp"))
  n_regions=length(param_values)
  assert_that(n_regions==length(shape_data$shapes))

  #Set map dimensions
  if(is.null(lat_max)){lat_max=shape_data$lat_max}
  if(is.null(lat_min)){lat_min=shape_data$lat_min}
  if(is.null(long_max)){long_max=shape_data$long_max}
  if(is.null(long_min)){long_min=shape_data$long_min}
  height_ll=lat_max-lat_min
  width_ll=long_max-long_min
  pixel_scale=pixels_max/max(height_ll,width_ll)
  width_px=width_ll*pixel_scale
  height_px=height_ll*pixel_scale

  #Assign parameter values within scale
  assert_that(min(param_values)>=min(scale))
  assert_that(max(param_values)<=max(scale))
  scale_values=rep(NA,length(param_values))
  for(i in 1:length(param_values)){
    scale_values[i]=findInterval(param_values[i],scale)
  }

  #Set colours
  n_scale_values=length(scale)
  ratio=length(colour_scale)/n_scale_values
  values=ratio*c(1:length(colour_scale))[c(1:n_scale_values)]
  for(i in 1:n_scale_values){values[i]=max(1,floor(values[i]))}
  colour_scale <- colour_scale[values]

  #Create legend labels
  legend_labels=rep("",n_scale_values-1)
  if(legend_format=="pc"){
    for(i in 1:(n_scale_values-1)){legend_labels[i]=paste(formatC(scale[i]*100,format="f",digits=legend_dp)," - ",
                                                          formatC(scale[i+1]*100,format="f",digits=legend_dp),sep="")}

  }
  if(legend_format=="f"){
    for(i in 1:(n_scale_values-1)){legend_labels[i]=paste(formatC(scale[i],format="f",digits=legend_dp)," - ",
                                                          formatC(scale[i+1],format="f",digits=legend_dp),sep="")}
  }
  if(legend_format=="e"){
    for(i in 1:(n_scale_values-1)){legend_labels[i]=paste(formatC(scale[i],format="e",digits=legend_dp)," - ",
                                                          formatC(scale[i+1],format="e",digits=legend_dp),sep="")}
  }

  #Create graph
  par(mar=c(1,1,1,1))
  if(is.null(output_file)==FALSE){png(filename=output_file,width=width_px,height=height_px)}
  {
    matplot(x=c(shape_data$long_min,shape_data$long_max),y=c(shape_data$lat_min,shape_data$lat_max),
            col=0,xlab="",ylab="",axes=display_axes,frame.plot=display_axes)
    for(n_region in 1:n_regions){
      plot(st_geometry(shape_data$shapes[[n_region]]),col=colour_scale[scale_values[n_region]],border="grey",add=TRUE)
    }
    legend(legend_position,legend=legend_labels,fill=colour_scale,cex=text_size,title=legend_title,ncol=legend_columns)
    title(main=map_title,cex=text_size)
  }
  if(is.null(output_file)==FALSE){dev.off()}
  par(mar=c(4,4,4,4))

  return(NULL)
}
