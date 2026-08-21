#' setArgs
#' @param script name of current script
#' @returns optparse class object
#' @export
setArgs <- function(script=NULL){

  if(is.null(script)){
    stop("no script name provided")
  }

  argparse <- optparse::OptionParser()

  if(script %in% c("downsampleBams","gridsearch","qdnaseqMod","qdnaseqRelToAbs","splitInput","computeSignatures")){
    argparse <- optparse::add_option(object = argparse,"--sampleName",action = "store",
                                     type = "character",default = NULL,
                                     help = "sample name string")
  }

  if(script %in% c("bamCheck","downsampleBams","qdnaseqMod")){
    argparse <- optparse::add_option(object = argparse,"--bam",action = "store",
                                     type = "character",default = NULL,
                                     help = "Input path for sequence alignment file")
  }

  if(script %in% c("bamCheck","downsampleBams")){
    argparse <- optparse::add_option(object = argparse,"--fileType",action = "store",
                                     type = "character",default = NULL,
                                     help = "file type of aligned file - 'BAM' or 'CRAM'")
  }

  if(script %in% c("combineResults","downsampleBams","gridsearch","qdnaseqRelToAbs")){
    argparse <- optparse::add_option(object = argparse,"--rds",action = "store",
                                     type = "character",default = NULL,
                                     help = "Path for rds input file")
  }

  if(script %in% c("autofit","combineResults","splitInput","computeSignatures")){
    argparse <- optparse::add_option(object = argparse,"--tsv",action = "store",
                                     type = "character",default = NULL,
                                     help = "Input path for TSV file")
  }

  if(script %in% c("downsampleBams","gridsearch","qdnaseqMod","qdnaseqRelToAbs")){
    argparse <- optparse::add_option(object = argparse,"--meta",action = "store",
                                     type = "character",default = NULL,
                                     help = "Path for metadata file")
  }

  if(script == "combineResults"){
    argparse <- optparse::add_option(object = argparse,"--tab",action = "store",
                                     type = "character",default = NULL,
                                     help = "Path for tab input file")
  }

  if(script %in% c("autofit","combineResults","gridsearch","qdnaseqRelToAbs","splitInput")){
    argparse <- optparse::add_option(object = argparse,"--tsvOut",action = "store",
                                     type = "character",default = NULL,
                                     help = "Output path for TSV file")
  }

  if(script %in% c("combineResults","qdnaseqMod","qdnaseqRelToAbs")){
    argparse <- optparse::add_option(object = argparse,"--rdsOut",action = "store",
                                     type = "character",default = NULL,
                                     help = "Path for rds output file")
  }

  if(script %in% c("combineResults","qdnaseqRelToAbs")){
    argparse <- optparse::add_option(object = argparse,"--tabOut",action = "store",
                                     type = "character",default = NULL,
                                     help = "Path for tab output file")
  }

  if(script %in% c("gridsearch","qdnaseqRelToAbs","computeSignatures")){
    argparse <- optparse::add_option(object = argparse,"--plotOut",action = "store",
                                     type = "character",default = NULL,
                                     help = "Path for plotting output file") |>
      optparse::add_option("--project",action = "store",type = "character",
                           default = NULL,help = "charater string for project name")
  }

  if(script == "bamCheck"){
    argparse <- optparse::add_option(object = argparse,"--okOut",action = "store",
                                     type = "character",default = NULL,
                                     help = "Path for bamCheck output file")
  }

  if(script %in% c("gridsearch","qdnaseqMod","qdnaseqRelToAbs","computeSignatures")){
    argparse <- optparse::add_option(object = argparse,"--genome",action = "store",
                                     type = "character",default = NULL,
                                     help = "string of genome build") |>
      optparse::add_option("--bin",action = "store",type = "double",
                           default = NULL,help = "bin size in kilobases")
  }

  if(script == "autofit"){
    argparse <- optparse::add_option(object = argparse,"--flagThreshold",action = "store",
                                     type = "double",default = NULL,
                                     help = "Threshold for flagging fits") |>
      optparse::add_option("--fitMethod",action = "store",type = "character",
                           default = NULL,help = "Fit method for autofit") |>
      optparse::add_option("--errorMetric",action = "store",type = "character",
                           default = NULL,help = "Error metric for autofit")
  }

  if(script == "downsampleBams"){
    argparse <- optparse::add_option(object = argparse,"--prplpu",action = "store",
                                     type = "logical",default = NULL,
                                     help = "preset ploidy-purity boolean") |>
      optparse::add_option("--reference",action = "store",type = "character",
                           default = NULL,help = "Path for genome reference fasta file (CRAM)") |>
      optparse::add_option("--bamOut",action = "store",type = "character",
                           default = NULL,help = "Path for sequence alignment output file")
  }

  if(script == "gridsearch"){
    argparse <- optparse::add_option(object = argparse,"--filtOut",action = "store",
                                     type = "character",default = NULL,
                                     help = "Path for filtered output file") |>
      optparse::add_option("--pdfOut",action = "store",type = "character",
                           default = NULL,help = "Path for sunrise output file") |>
      optparse::add_option("--fitOut",action = "store",type = "character",
                         default = NULL,help = "Path for fitting output file") |>
      optparse::add_option("--pl_min",action = "store",type = "double",
                         default = NULL,help = "Minimum ploidy") |>
      optparse::add_option("--pl_max",action = "store",type = "double",
                         default = NULL,help = "Maximum ploidy") |>
      optparse::add_option("--pu_min",action = "store",type = "double",
                         default = NULL,help = "Minimum purity") |>
      optparse::add_option("--pu_max",action = "store",type = "double",
                         default = NULL,help = "Maximum purity") |>
      optparse::add_option("--af_cutoff",action = "store",type = "double",
                         default = NULL,help = "Allele frequency difference cutoff") |>
      optparse::add_option("--filter_underpowered",action = "store",type = "logical",
                         default = NULL,help = "Filter underpowered fits boolean") |>
      optparse::add_option("--filter_homozygous",action = "store",type = "logical",
                         default = NULL,help = "Filter homozygous loss proportion") |>
      optparse::add_option("--homozygous_prop",action = "store",type = "double",
                         default = NULL,help = "Genome proportion in Basepairs") |>
      optparse::add_option("--homozygous_thrsh",action = "store",type = "double",
                         default = NULL,help = "Threshold for homozygous loss")
  }

  if(script == "qdnaseqMod"){
    argparse <- optparse::add_option(object = argparse,"--autoSmooth",action = "store",
                                     type = "logical",default = NULL,
                                     help = "autosmooth boolean") |>
      optparse::add_option("--smoothThreshold",action = "store",type = "double",
                           default = NULL,help = "Threshold for smoothing") |>
      optparse::add_option("--use_seed",action = "store",type = "logical",
                           default = NULL,help = "use seed boolean") |>
      optparse::add_option("--seed_val",action = "store",type = "character",
                           default = NULL,help = "seed value") |>
      optparse::add_option("--pairedEnd",action = "store",type = "logical",
                           default = NULL,help = "pairedEnd boolean")
  }

  if(script == "computeSignatures"){
    argparse <- optparse::add_option(object = argparse,"--tsvDrews",action = "store",
                                     type = "character",default = NULL,
                                     help = "Drews signature output path") |>
      optparse::add_option("--tsvMac",action = "store",type = "character",
                           default = NULL,help = "Mac signature output path")
  }
  return(argparse)
}
