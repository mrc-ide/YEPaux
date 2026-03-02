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
map_shapes_load <- function(regions=c(), shapefiles=c(), region_label_type=""){

  #TODO - Change function to deal with inputs where "geometry" called something else

  assert_that(is.character(regions))
  assert_that(is.character(shapefiles))
  assert_that(is.character(region_label_type))

  n_regions=length(regions)
  for(i in 1:length(shapefiles)){
    shape_data=read_sf(shapefiles[i])
    assert_that(region_label_type %in% names(shape_data), msg=paste0("Region label not found in ",
                                                                    shapefiles[i]))
    if(i==1){shape_data_all=st_sf(data.frame(region=rep(NA, n_regions),
                                             geometry=rep(shape_data$geometry[1], n_regions)))}
    file_regions=shape_data[[match(region_label_type, names(shape_data))]]

    for(n_region in 1:n_regions){
      region=regions[n_region]
      k=match(region, file_regions)
      if(is.na(k)==FALSE){
        shape_data_all$region[n_region]=region
        shape_data_all$geometry[n_region]=shape_data$geometry[k]
      }
    }
  }
  assert_that(all(shape_data_all$region==regions), msg="Missing region data")

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
#' @param text_size Size of text to appear in legend and titles
#' @param display_axes TRUE/FALSE flag indicating whether to frame map and display latitude/longitude axes
#' @param border_colour_regions Colour to use for borders of regions. Set to NA if borders to be invisible.
#' @param ... Additional optional parameters: \cr
#'    scale_manual Vector of scale intervals to use for param_values \cr
#'    colour_scale_manual Vector of colours with size greater than or equal to scale_manual - used to convert scale_manual to colours\cr
#'    lat_max, lat_min, long_max, long_min: borders if not to be set default \cr
#'    additional_border_shapes: Shape data for optional additional borders \cr
#'    border_colour_additional: colour to use for additional borders if any. \cr
#     map_title: Title to show above map \cr
#     legend_title: Title to show above legend \cr
#     legend_format: Number format to use for scale values in legend if used \cr
#     legend_dp: Number of decimal places to use in scale values in legend \cr
#     legend_position: TBA \cr
#     legend_columns: Number of columns in which to display legend values (TBA) \cr
#' '
#' @export
#'
create_map <- function(shape_data=list(), param_values=c(), text_size=1,
                       display_axes=FALSE, border_colour_regions="grey", ...){

  #TODO - fix text size

  assert_that(is.list(shape_data))
  assert_that(is.numeric(param_values))
  assert_that(is.logical(display_axes))
  n_regions=length(param_values)
  assert_that(n_regions==length(shape_data$geometry))
  ap=list(...) #Get additional optional parameters

  #Set map dimensions
  bbox=st_bbox(shape_data)
  if(is.null(ap$lat_max)){ap$lat_max=bbox$ymax}
  if(is.null(ap$lat_min)){ap$lat_min=bbox$ymin}
  if(is.null(ap$long_max)){ap$long_max=bbox$xmax}
  if(is.null(ap$long_min)){ap$long_min=bbox$xmin}

  #Assign parameter values within scale
  if(is.null(ap$scale_manual)==FALSE){
    assert_that(is.numeric(ap$scale_manual))
    assert_that(min(param_values, na.rm=TRUE)>=min(ap$scale_manual))
    scale_values=rep(NA, length(param_values))
    for(i in 1:length(param_values)){
      scale_values[i]=findInterval(param_values[i], ap$scale_manual)
    }
    map_values=as.character(ap$scale_manual[scale_values])
    n_intervals=length(ap$scale_manual)
    ratio=length(ap$colour_scale_manual)/n_intervals
    values=ratio*c(1:length(ap$colour_scale_manual))[c(1:n_intervals)]
    for(i in 1:n_intervals){values[i]=max(1, floor(values[i]))}
    palette_vector=ap$colour_scale_manual[values]
  } else {
    map_values = param_values
  }

  #Create graph (ggplot, new)
  map_output <- ggplot()
  map_output <- map_output + geom_sf(data = shape_data,
                                     mapping = aes(fill=map_values),
                                     colour = border_colour_regions,
                                     show.legend=TRUE)
  map_output <- map_output + xlim(ap$long_min, ap$long_max) +  ylim(ap$lat_min, ap$lat_max)
  map_output <- map_output + labs(fill = ap$legend_title, title = ap$map_title)
  if(display_axes){map_output <- map_output+theme_linedraw()}else{map_output <- map_output+theme_void()}
  if(is.null(ap$legend_position)==FALSE){
    map_output <- map_output + theme(legend.position = ap$legend_position)
  }
  if(is.null(ap$scale_manual)==FALSE){
    map_output <- map_output + scale_fill_manual(limits=as.character(ap$scale_manual[c(1:n_intervals)]),
                                                 aesthetics="fill",
                                                 values=palette_vector,
                                                 breaks=ap$scale_manual[c(1:n_intervals)],
                                                 na.value = "grey50")
  } else {
    map_output <- map_output + scale_fill_viridis_c(option="magma",
                                                    na.value = "grey50")
  }
  if(is.null(ap$additional_border_shapes)==FALSE){
    map_output <- map_output + geom_sf(data=ap$additional_border_shapes$geometry,
                                       fill = NA,
                                       colour=ap$border_colour_additional,
                                       show.legend=NA)
  }
  map_output <- map_output + theme(text = element_text(size = text_size))

  return(map_output)
}
