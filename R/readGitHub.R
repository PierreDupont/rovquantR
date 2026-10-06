#' Read files from GitHub 
#'
#' The \code{readGitHub} function identifies and loads specific files hosted on GitHub.
#' 
#' @name readGitHub
#'
#' @param data.dir A character string denoting the data directory path.
#' 
#' @return The data loaded
#'
#' @author Pierre Dupont
#'
#' @importFrom gh gh
#' @importFrom terra rast 
#' @importFrom sf read_sf st_length st_coordinates st_centroid
#' @importFrom dplyr mutate filter
#' 
NULL
#' @rdname readGitHub
#' @export
readGitHub <- function(
   report = "fagrapport112_Wolverine2016-2025",
   subdirectory = "output/rasters/RasterForRovbase",
   extension = ".tif"
   ){
  
## 1. Identify files

# Fetch the recursive tree of the repository
repo_tree <- gh::gh("GET /repos/:owner/:repo/git/trees/:tree_sha",
                    owner = "richbi",
                    repo = "RovQuantPublic",
                    tree_sha = "main",      # Can be a branch name or commit SHA
                    recursive = 1)          # 1 tells GitHub to return all nested subfolders

# Extract the file paths into a vector
gh_files <- sapply(repo_tree$tree, function(x) x$path) 

# Identify the paths to the correct report
path <- ifelse(is.null(subdirectory),
               report,
               file.path(report,subdirectory))
gh_files <- gh_files[grep(path, gh_files)]

# Identify the paths to the correct report
gh_files <- gh_files[grep(extension, gh_files)]

# list raw GitHub URLs
base_url <- "https://github.com/richbi/RovQuantPublic/raw/refs/heads/main" 
raw_url <- file.path(base_url, gh_files)



## 2. load .tiff

if( extension = ".tif"){
  
  
  }
## Function to read .tif from GitHub
readGitHub_tif <- function(raw_url){
  # Create a temporary file path
  temp_tif <- tempfile(fileext = ".tif")
  
  # Download the binary file from GitHub
  download.file(raw_url,
                destfile = temp_tif,
                mode = "wb")
  
  # Read and return the GeoTIFF data
  spatial_raster <- rast(temp_tif)
  return(spatial_raster)
}
out <- lapply(raw_url, readGitHub_tif)



## 3. load .csv
out <- lapply(raw_url, read.csv)


# Inspect the data
print(spatial_raster)
plot(spatial_raster)


##-- Output
return(data)
}



