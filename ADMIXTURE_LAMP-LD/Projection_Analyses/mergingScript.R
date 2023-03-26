#' Write SNP or sample IDs to file
#' 
#' @param ids Vector of SNP or Sample IDs
#' @param file File name to write to
#' @param update Logical, Default is FALSE. Skip writing file if it exists
#'   
#' @param sample Logical. Should be TRUE for samples (the default) and FALSE for
#'   SNPs
#'   
#' @return NULL
#' @import utils
#' @export
plink_write_ids <- function(ids, file, update, sample = TRUE) {
  if (sample) {
    out <- cbind(ids, ids)
  } else {
    out <- ids
  }
  if (update | !file.exists(file)) {
    utils::write.table(
      out,
      file,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  }
}

#' Subset samples and variants using PLINK
#'
#' @param bfile              [\code{string}]\cr
#'                           The basename of the binary PLINK files.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param remove             [\code{string}]\cr
#'                           File name of a list of sample FIDs and IIDs to remove.
#' @param keep               [\code{string}]\cr
#'                           File name of a list of sample FIDs and IIDs to keep.
#' @param exclude            [\code{string}]\cr
#'                           File name of a list of variant names to remove.
#' @param extract            [\code{string}]\cr
#'                           File name of a list of variant names to keep.
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK.
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @details See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate
#'
plink_subset <- function(bfile, output.prefix, remove, keep, exclude, extract, ...,
                         bed.file = NULL, bim.file = NULL, fam.file = NULL,
                         exec = "plink2",
                         num.threads,
                         memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(bfile)) {
    checkmate::assert_file(sprintf("%s.bed", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.bim", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.fam", bfile), add = assertions)
    input <- sprintf("--bfile %s", bfile)
  } else {
    checkmate::assert_file(bed.file, add = assertions)
    checkmate::assert_file(bim.file, add = assertions)
    checkmate::assert_file(fam.file, add = assertions)
    input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file)
  }
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  if (missing(remove)) {
    remove <- ""
  } else {
    checkmate::assert_file(remove, add = assertions)
    remove <- sprintf("--remove %s", remove)
  }
  if (missing(keep)) {
    keep <- ""
  } else {
    checkmate::assert_file(keep, add = assertions)
    keep <- sprintf("--keep %s", keep)
  }
  
  if (missing(exclude)) {
    exclude <- ""
  } else {
    checkmate::assert_file(exclude, add = assertions)
    exclude <- sprintf("--exclude %s", exclude)
  }
  if (missing(extract)) {
    extract <- ""
  } else {
    checkmate::assert_file(extract, add = assertions)
    extract <- sprintf("--extract %s", extract)
  }
  
  assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  }  
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {  
    memory = max(5000,     
                 -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000),           
                      -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE),                 
                 na.rm = TRUE)  
  }  
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  system_call(
    bin = exec,
    args = c(input,
             "--threads", num.threads,
             "--memory", memory,
             remove,
             keep,
             extract,
             exclude,
             "--make-bed",
             "--keep-allele-order",
             "--allow-extra-chr",
             "--out", output.prefix, ...)
  )
  
}

#' Normalize and convert VCF files to PLINK binary files
#'
#' @param vcf.file           [\code{string}]\cr
#'                           The VCF file path.
#' @param ref.file           [\code{string}]\cr
#'                           A human reference genome \code{fasta} file to normalize indels against.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK.
#' @param build              [\code{string}]\cr
#'                           Human genome build code to split the X chromosome into pseudo-autosomal region and pure X. Default is 'b37' (alias 'hg19').
#' @param bcftools.exec      [\code{string}]\cr
#'                           Path of bcftools executable.
#' @param plink.exec         [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @details Based on the best practices on \url{http://apol1.blogspot.de/2014/11/best-practice-for-converting-vcf-files.html}. The procedure will take the VCF file, strip the variant IDs, split multi-allelic sites into bi-allelic sites, assign names to make sure indels will not become ambiguous, and finally convert to PLINK format. See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate
#'
plink_normalize_vcf_conversion <- function(vcf.file, ref.file, output.prefix, ..., 
                                           build = "b37",
                                           bcftools.exec = "bcftools", plink.exec = "plink2",
                                           num.threads,
                                           memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(vcf.file, add = assertions)
  checkmate::assert_file(ref.file, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  checkmate::assertChoice(build, c("b36", "hg18", "b37", "hg19", "b38", "hg38"), add = assertions)
  
  assert_command(bcftools.exec, add = assertions)
  assert_command(plink.exec, add = assertions)
  
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  } 
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {  
    memory = max(5000,       
                 -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000),              
                      -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE),                   na.rm = TRUE)  
  } 
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  # Run PLINK
  system_call(
    bin = bcftools.exec,
    args = c("norm", "-Ou", "-m", "-any", vcf.file,
             "|",
             bcftools.exec, "norm", "-Ou", "-f", ref.file,
             "|",
             bcftools.exec, "annotate", "-Ob", "-x", "ID", "-I", "+'%CHROM:%POS:%REF:%ALT'",
             "|",
             plink.exec, "--bcf", "/dev/stdin",
             "--vcf-idspace-to", "_",
             "--keep-allele-order",
             "--allow-extra-chr",
             "--split-x", build, "no-fail",
             "--threads", num.threads,
             "--memory", memory,
             "--make-bed",
             "--out", output.prefix, ...)
  )
  
}

#' Convert VCF files to PLINK binary files
#' 
#' Be sure the VCF file is well prepared!
#'
#' @param vcf.file           [\code{string}]\cr
#'                           The VCF file path.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK.
#' @param build              [\code{string}]\cr
#'                           Human genome build code to split the X chromosome into pseudo-autosomal region and pure X. Default is 'b37' (alias 'hg19').
#' @param plink.exec         [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @details Based on the best practices on \url{http://apol1.blogspot.de/2014/11/best-practice-for-converting-vcf-files.html}. The procedure will take the VCF file and convert it to PLINK format. Be sure the VCF file is well prepared! See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate
#' 
#' @seealso plink_normalize_vcf_conversion
#' @seealso vcf_normalization
#'
plink_vcf_conversion <- function(vcf.file, output.prefix, ..., 
                                 build = "b37",
                                 plink.exec = "plink2",
                                 num.threads,
                                 memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(vcf.file, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  checkmate::assertChoice(build, c("b36", "hg18", "b37", "hg19", "b38", "hg38"), add = assertions)
  
  assert_command(plink.exec, add = assertions)
  
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  } 
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {  
    memory = max(5000,       
                 -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000),              
                      -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE),                   na.rm = TRUE)  
  } 
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  # Run PLINK
  system_call(
    bin = plink.exec,
    args = c("--vcf", vcf.file,
             "--vcf-idspace-to", "_",
             "--keep-allele-order",
             "--allow-extra-chr",
             "--split-x", build, "no-fail",
             "--threads", num.threads,
             "--memory", memory,
             "--make-bed",
             "--out", output.prefix, ...)
  )
  
}

#' Merge two PLINK datasets
#'
#' @param first.prefix     [\code{string}]\cr
#'                         The basename of the first binary PLINK file set.
#' @param second.prefix    [\code{string}]\cr
#'                         The basename of the second binary PLINK file set.
#' @param merge.mode       [\code{int}]\cr
#'                         Merge mode.
#' @param output.prefix    [\code{string}]\cr
#'                         The basename of the output binary PLINK file set.
#' @param ...              [\code{character}]\cr
#'                         Additional arguments passed to PLINK.
#' @param first.bed.file   [\code{string}]\cr
#'                         Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param first.bim.file   [\code{string}]\cr
#'                         Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param first.fam.file   [\code{string}]\cr
#'                         Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param second.bed.file  [\code{string}]\cr
#'                         Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param second.bim.file  [\code{string}]\cr
#'                         Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param second.fam.file  [\code{string}]\cr
#'                         Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec             [\code{string}]\cr
#'                         Path of PLINK executable.
#' @param num.threads      [\code{int}]\cr
#'                         Number of CPUs usable by PLINK.
#'                         Default is determined by SLURM environment variables and at least 1.
#' @param memory           [\code{int}]\cr
#'                         Memory for PLINK in Mb.
#'                         Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @details See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate
#'
plink_merge <- function(first.prefix, second.prefix, merge.mode,
                        output.prefix,
                        ...,
                        first.bed.file = NULL, first.bim.file = NULL, first.fam.file = NULL,
                        second.bed.file = NULL, second.bim.file = NULL, second.fam.file = NULL,
                        exec = "plink2",
                        num.threads,
                        memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(first.prefix)) {
    checkmate::assert_file(sprintf("%s.bed", first.prefix), add = assertions)
    checkmate::assert_file(sprintf("%s.bim", first.prefix), add = assertions)
    checkmate::assert_file(sprintf("%s.fam", first.prefix), add = assertions)
    first <- sprintf("--bfile %s", first.prefix)
  } else {
    checkmate::assert_file(first.bed.file, add = assertions)
    checkmate::assert_file(first.bim.file, add = assertions)
    checkmate::assert_file(first.fam.file, add = assertions)
    first <- sprintf("--bed %s --bim %s --fam %s", first.bed.file, first.bim.file, first.fam.file)
  }
  
  if (!missing(second.prefix)) {
    checkmate::assert_file(sprintf("%s.bed", second.prefix), add = assertions)
    checkmate::assert_file(sprintf("%s.bim", second.prefix), add = assertions)
    checkmate::assert_file(sprintf("%s.fam", second.prefix), add = assertions)
    second <- sprintf("--bmerge %s", second.prefix)
  } else {
    checkmate::assert_file(second.bed.file, add = assertions)
    checkmate::assert_file(second.bim.file, add = assertions)
    checkmate::assert_file(second.fam.file, add = assertions)
    second <- sprintf("--bmerge %s %s %s", second.bed.file, second.bim.file, second.fam.file)
  }
  
  checkmate::assert_int(merge.mode, lower = 1, upper = 7, add = assertions)
  
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE) 
  }   
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {   
    memory = max(5000,           
                 -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000),                    
                      -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE),      
                 na.rm = TRUE) 
  }  
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  # Make a merge with PLINK
  system_call(
    bin = exec,
    args = c(first,
             second,
             "--threads", num.threads,
             "--memory", memory,
             "--keep-allele-order",
             "--allow-extra-chr",
             "--merge-mode ", merge.mode,
             "--out", output.prefix,
             ...)
  )
  
}

#' LD pruning with PLINK
#'
#' @param bfile              [\code{string}]\cr
#'                           The basename of the binary PLINK files.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param window.size        [\code{int}]\cr
#'                           Window size.
#' @param kb.window          [\code{boolean}]\cr
#'                           Is window size in kb units?
#' @param step.size          [\code{int}]\cr
#'                           Step size in variant counts.
#' @param chr                [\code{number}]\cr
#'                           Restrict LD pruning to specified chromosome.
#' @param threshold          [\code{number}]\cr
#'                           Pairwise \eqn{r^2} threshold.
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK.
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @details See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate
#'
plink_ld_pruning <- function(bfile, output.prefix, 
                             window.size, kb.window = FALSE, step.size, threshold, chr, ...,
                             bed.file = NULL, bim.file = NULL, fam.file = NULL,
                             exec = "plink2",
                             num.threads,
                             memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(bfile)) {
    checkmate::assert_file(sprintf("%s.bed", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.bim", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.fam", bfile), add = assertions)
    input <- sprintf("--bfile %s", bfile)
  } else {
    checkmate::assert_file(bed.file, add = assertions)
    checkmate::assert_file(bim.file, add = assertions)
    checkmate::assert_file(fam.file, add = assertions)
    input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file)
  }
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  checkmate::assert_int(window.size, lower = 2, add = assertions)
  checkmate::assert_logical(kb.window, add = assertions)
  if (kb.window) {
    window_size <- sprintf("%dkb", window.size)
    message("When window size is in kb units, step size is set to 1.")
    step_size <- 1
  } else {
    window_size <- window.size
    step_size <- step.size
  }
  checkmate::assert_int(step_size, lower = 1, add = assertions)
  
  checkmate::assert_number(threshold, lower = 0, upper = 1, add = assertions)
  
  if (!missing(chr)) {
    checkmate::assert_int(chr, lower = 1, add = assertions)
    chr <- sprintf("--chr %d", chr)
  } else {
    chr <- ""
  }
  
  assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {  
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE) 
  }  
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) { 
    memory = max(5000,       
                 -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000),                   
                      -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE),                
                 na.rm = TRUE)  
  }   
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  system_call(
    bin = exec,
    args = c(input,
             "--threads", num.threads,
             "--memory", memory,
             "--keep-allele-order",
             "--allow-extra-chr",
             chr,
             "--indep-pairwise", window_size, step_size, threshold,
             "--out", output.prefix, ...)
  )
  
}

#' Sex imputation with PLINK
#'
#' @param bfile              [\code{string}]\cr
#'                           The basename of the binary PLINK files.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param f.values           [\code{numeric(2)}]\cr
#'                           \code{numeric} of length \code{2} with minimum F coefficients to impute females (\code{[1]}) or males (\code{[2]}).
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK.
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @details See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate tools
#'
plink_sex_imputation <- function(bfile, output.prefix, 
                                 f.values, ...,
                                 bed.file = NULL, bim.file = NULL, fam.file = NULL,
                                 exec = "plink2",
                                 num.threads,
                                 memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(bfile)) {
    checkmate::assert_file(sprintf("%s.bed", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.bim", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.fam", bfile), add = assertions)
    input <- sprintf("--bfile %s", bfile)
  } else {
    checkmate::assert_file(bed.file, add = assertions)
    checkmate::assert_file(bim.file, add = assertions)
    checkmate::assert_file(fam.file, add = assertions)
    input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file)
  }
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  checkmate::assert_numeric(f.values, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, all.missing = FALSE, len = 2, null.ok = FALSE, add = assertions)
  
  assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {     
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)   
  }   
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {     
    memory = max(5000,                   
                 -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000),                       -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE),                   
                 na.rm = TRUE)   
  }   
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  system_call(
    bin = exec,
    args = c(input,
             "--threads", num.threads,
             "--memory", memory,
             "--impute-sex", f.values,
             "--keep-allele-order",
             "--allow-extra-chr",
             "--make-bed",
             "--out", output.prefix, ...)
  )
  
}

#' Merge multiple PLINK datasets
#'
#' @param bfile              [\code{string}]\cr
#'                           The basename of the binary PLINK files.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param merge.list         [\code{string}]\cr
#'                           File listing other PLINK datasets to be merged.
#' @param merge.mode         [\code{int}]\cr
#'                           Merge mode.
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK.
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_NODE} and \code{num.threads * SLURM_MEM_PER_CPU} and at least 5000.
#'
#' @details See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate tools
#'
plink_merge_list <- function(bfile, output.prefix, 
                             merge.list, merge.mode, ...,
                             bed.file, bim.file, fam.file,
                             exec = "plink2",
                             num.threads,
                             memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(bfile)) {
    checkmate::assert_file(sprintf("%s.bed", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.bim", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.fam", bfile), add = assertions)
    input <- sprintf("--bfile %s", bfile)
  } else if (!any(c(missing(bed.file), missing(bim.file), missing(fam.file)))) {
    checkmate::assert_file(bed.file, add = assertions)
    checkmate::assert_file(bim.file, add = assertions)
    checkmate::assert_file(fam.file, add = assertions)
    input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file)
  } else {
    input <- ""
  }
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  checkmate::assert_file(merge.list, add = assertions)
  
  if (!missing(merge.mode)) {
    checkmate::assert_int(merge.mode, lower = 1, upper = 7, add = assertions)
    merge.mode <- sprintf("--merge-mode %d", merge.mode)
  } else {
    merge.mode <- ""
  }
  
  assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")),
                       na.rm = TRUE)
  }
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  
  if (missing(memory)) {
    memory = max(5000, 
                 -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000), 
                      -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE), 
                 na.rm = TRUE)
  }
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  system_call(
    bin = exec,
    args = c(input,
             "--threads", num.threads,
             "--memory", memory,
             "--merge-list", merge.list,
             merge.mode,
             "--keep-allele-order",
             "--allow-extra-chr",
             "--make-bed",
             "--out", output.prefix, ...)
  )
  
}

#' Remove duplicated variants
#' 
#' Uses \code{cut}, \code{sort}, \code{uniq} and \code{awk} to find duplicated markers excludes them using PLINK.
#'
#' @param bfile              [\code{string}]\cr
#'                           The basename of the binary PLINK files.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK.
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_NODE} and \code{num.threads * SLURM_MEM_PER_CPU} and at least 5000.
#'
#' @details See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate tools
#'
plink_dedup <- function(bfile, output.prefix, 
                        ...,
                        bed.file = NULL, bim.file = NULL, fam.file = NULL,
                        exec = "plink2",
                        num.threads,
                        memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(bfile)) {
    checkmate::assert_file(sprintf("%s.bed", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.bim", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.fam", bfile), add = assertions)
    input <- sprintf("--bfile %s", bfile)
    input_prefix <- bfile
    bim_file <- sprintf("%s.bim", bfile)
  } else {
    checkmate::assert_file(bed.file, add = assertions)
    checkmate::assert_file(bim.file, add = assertions)
    checkmate::assert_file(fam.file, add = assertions)
    input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file)
    input_prefix <- sub("\\.bim", "", bed.file)
    bim_file <- bim.file
  }
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")),
                       na.rm = TRUE)
  }
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  
  if (missing(memory)) {
    memory = max(5000, 
                 -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000), 
                      -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE), 
                 na.rm = TRUE)
  }
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  dups_file <- sprintf("%s.dups", input_prefix)
  
  # Find duplicated markers
  system_call(
    bin = "cut",
    args = c("-f2",
             bim_file,
             "|", "sort", 
             "|", "uniq -c", 
             "|", "awk", "'$1>=2 {print $2}'", 
             ">", dups_file)
  )
  
  # Remove duplicated markers
  system_call(
    bin = exec,
    args = c(input, ...,
             "--keep-allele-order",
             "--allow-extra-chr",
             "--exclude", dups_file,
             "--make-bed",
             "--out", output.prefix)
  )
  
}

#' Remove very long INDELS
#' 
#' Uses \code{cut}, \code{sort}, \code{uniq} and \code{awk} to find very long INDELS causing PLINK to be very memory hungry during analyses and excludes them using PLINK.
#'
#' @param bfile              [\code{string}]\cr
#'                           The basename of the binary PLINK files.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param max_length         [\code{integer}]\cr
#'                           The maximum length of any allele.
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK.
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_NODE} and \code{num.threads * SLURM_MEM_PER_CPU} and at least 5000.
#'
#' @details See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate tools
#'
plink_rm_long_indels <- function(bfile, output.prefix, 
                                 max_length,
                                 ...,
                                 bed.file = NULL, bim.file = NULL, fam.file = NULL,
                                 exec = "plink2",
                                 num.threads,
                                 memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(bfile)) {
    checkmate::assert_file(sprintf("%s.bed", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.bim", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.fam", bfile), add = assertions)
    input <- sprintf("--bfile %s", bfile)
    input_prefix <- bfile
    bim_file <- sprintf("%s.bim", bfile)
  } else {
    checkmate::assert_file(bed.file, add = assertions)
    checkmate::assert_file(bim.file, add = assertions)
    checkmate::assert_file(fam.file, add = assertions)
    input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file)
    input_prefix <- sub("\\.bim", "", bed.file)
    bim_file <- bim.file
  }
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  checkmate::assert_int(max_length, add = assertions)
  
  assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")),
                       na.rm = TRUE)
  }
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  
  if (missing(memory)) {
    memory = max(5000, 
                 -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000), 
                      -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE), 
                 na.rm = TRUE)
  }
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  long_indels_file <- sprintf("%s.longindels", input_prefix)
  
  # Find very long indels
  system_call(
    bin = "awk",
    args = c(sprintf("'{if (length($5) > %d || length($6) > %d) print $2}'", max_length, max_length),
             bim_file,
             "|", "sort", 
             "|", "uniq", 
             ">", long_indels_file)
  )
  
  # Remove very long indels
  system_call(
    bin = exec,
    args = c(input, ...,
             "--keep-allele-order",
             "--allow-extra-chr",
             "--exclude", long_indels_file,
             "--make-bed",
             "--out", output.prefix)
  )
  
}

#' Remove regions of high LD
#' 
#' There are regions of long-range, high linkage diequilibrium in the human genome. 
#' These regions should be excluded when performing certain analyses such as principal component analysis on genotype data.
#'
#' @param bfile              [\code{string}]\cr
#'                           The basename of the binary PLINK files.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param build              [\code{integer}]\cr
#'                           The genome build version of PLINK file. Default b37.
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK.
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_NODE} and \code{num.threads * SLURM_MEM_PER_CPU} and at least 5000.
#'
#' @details See \url{https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate tools
#'
plink_rm_high_ld <- function(bfile, output.prefix, 
                             build = "b37",
                             ...,
                             bed.file = NULL, bim.file = NULL, fam.file = NULL,
                             exec = "plink2",
                             num.threads,
                             memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(bfile)) {
    checkmate::assert_file(sprintf("%s.bed", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.bim", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.fam", bfile), add = assertions)
    input <- sprintf("--bfile %s", bfile)
    input_prefix <- bfile
    bim_file <- sprintf("%s.bim", bfile)
  } else {
    checkmate::assert_file(bed.file, add = assertions)
    checkmate::assert_file(bim.file, add = assertions)
    checkmate::assert_file(fam.file, add = assertions)
    input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file)
    input_prefix <- sub("\\.bim", "", bed.file)
    bim_file <- bim.file
  }
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  checkmate::assert_choice(build, choices = c("b37", "b36"), add = assertions)
  
  assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")),
                       na.rm = TRUE)
  }
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  
  if (missing(memory)) {
    memory = max(5000, 
                 -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000), 
                      -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE), 
                 na.rm = TRUE)
  }
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  if (build == "b37") {
    high_ld_regions <- fread(
      "Chr	Start	Stop
      1	48000000	52000000
      2	86000000	100500000
      2	134500000	138000000
      2	183000000	190000000
      3	47500000	50000000
      3	83500000	87000000
      3	89000000	97500000
      5	44500000	50500000
      5	98000000	100500000
      5	129000000	132000000
      5	135500000	138500000
      6	25000000	35000000
      6	57000000	64000000
      6	140000000	142500000
      7	55000000	66000000
      8	7000000	13000000
      8	43000000	50000000
      8	112000000	115000000
      10	37000000	43000000
      11	46000000	57000000
      11	87500000	90500000
      12	33000000	40000000
      12	109500000	112000000
      20	32000000	34500000")
  } else if (build == "b36") {
    high_ld_regions <- fread(
      "Chr	Start	Stop
      1	48060567	52060567
      2	85941853	100407914
      2	134382738	137882738
      2	182882739	189882739
      3	47500000	50000000
      3	83500000	87000000
      3	89000000	97500000
      5	44500000	50500000
      5	98000000	100500000
      5	129000000	132000000
      5	135500000	138500000
      6	25500000	33500000
      6	57000000	64000000
      6	140000000	142500000
      7	55193285	66193285
      8	8000000	12000000
      8	43000000	50000000
      8	112000000	115000000
      10	37000000	43000000
      11	46000000	57000000
      11	87500000	90500000
      12	33000000	40000000
      12	109521663	112021663
      20	32000000	34500000
      23	14150264	16650264
      23	25650264	28650264
      23	33150264	35650264
      23	55133704	60500000
      23	65133704	67633704
      23	71633704	77580511
      23	80080511	86080511
      23	100580511	103080511
      23	125602146	128102146
      23	129102146	131602146")
  }
  high_ld_regions[, SetID := sprintf("HiLD%d", .I)]
  
  high_ld_regions_file <- tempfile()
  on.exit(unlink(high_ld_regions_file, recursive = TRUE, force = TRUE), add = TRUE)
  
  fwrite(high_ld_regions, high_ld_regions_file, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # Remove very long indels
  system_call(
    bin = exec,
    args = c(input, ...,
             "--keep-allele-order",
             "--allow-extra-chr",
             "--exclude", "range", high_ld_regions_file,
             "--make-bed",
             "--out", output.prefix)
  )
  
}

#' Quality control on marker level with PLINK
#'
#' @param bfile              [\code{string}]\cr
#'                           The basename of the binary PLINK files.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param call.rate          [\code{double}]\cr
#'                           Filter out all variants with call rates falling below the provided threshold.
#' @param min.maf            [\code{double}]\cr
#'                           Filter out all variants with minor allele frequency below the provided threshold.
#' @param max.maf            [\code{double}]\cr                          
#'                           Filter out all variants with minor allele frequency above the provided threshold.
#' @param hwe.pval           [\code{double}]\cr
#'                           Filter out all variants which have Hardy-Weinberg equilibrium exact test p-value below the provided threshold.
#' @param mid.p              [\code{double}]\cr                          
#'                           Apply the mid-p adjustment? See Details.
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK.
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @details The 'mid.p' modifier applies the mid-p adjustment described in Graffelman J, Moreno V (2013). 
#' The mid-p adjustment tends to bring the null rejection rate in line with the nominal p-value, and also reduces the filter's tendency to favor retention of variants with missing data. 
#' PLINK recommends its use.
#' See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector and number of markers excluded per criteria.
#' @export
#'
#' @import checkmate tools
#'
plink_marker_qc <- function(bfile, output.prefix, 
                            call.rate, min.maf, max.maf, hwe.pval, mid.p = TRUE, ...,
                            bed.file = NULL, bim.file = NULL, fam.file = NULL,
                            exec = "plink2",
                            num.threads,
                            memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(bfile)) {
    checkmate::assert_file(sprintf("%s.bed", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.bim", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.fam", bfile), add = assertions)
    input <- sprintf("--bfile %s", bfile)
  } else {
    checkmate::assert_file(bed.file, add = assertions)
    checkmate::assert_file(bim.file, add = assertions)
    checkmate::assert_file(fam.file, add = assertions)
    input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file)
  }
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  if (!missing(call.rate)) {
    checkmate::assert_number(call.rate, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    geno <- sprintf("--geno %f", 1 - call.rate)
  } else {
    geno <- NULL
  }
  if (!missing(min.maf)) {
    checkmate::assert_number(min.maf, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    min_maf <- sprintf("--maf %f", min.maf)
  } else {
    min_maf <- NULL
  }
  if (!missing(max.maf)) {
    checkmate::assert_number(max.maf, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    max_maf <- sprintf("--max-maf %f", max.maf)
  } else {
    max_maf <- NULL
  }
  if (!missing(hwe.pval)) {
    checkmate::assert_number(hwe.pval, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    hwe_pval <- sprintf("--hwe %g", hwe.pval)
    if (mid.p) {
      hwe_pval <- sprintf("%s midp", hwe_pval)
    }
  } else {
    hwe_pval <- NULL
  }
  
  assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {     
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)   
  }   
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {     
    memory = max(5000,                   
                 -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000),                       
                      -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE),                   
                 na.rm = TRUE)   
  }   
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  log <- system_call(
    bin = exec,
    args = c(input,
             "--threads", num.threads,
             "--memory", memory,
             geno, 
             min_maf, 
             max_maf, 
             hwe_pval,
             "--keep-allele-order",
             "--allow-extra-chr",
             "--make-bed",
             "--out", output.prefix, ...)
  )
  
  num_geno_rm <- max(0, as.integer(gsub("^(\\d+) variants? removed due to missing genotype data.*", "\\1", grep("^\\d+ variants? removed due to missing genotype data.*", log, value = TRUE))))
  num_hwe_rm <- max(0, as.integer(gsub("^--hwe: (\\d+) variants?.*", "\\1", grep("^--hwe: \\d", log, value = TRUE))))
  num_maf_rm <- max(0, as.integer(gsub("^(\\d+) variants? removed due to minor allele threshold.*", "\\1", grep("^\\d+ variants? removed due to minor allele threshold.*", log, value = TRUE))))
  
  return(
    list(
      num_marker_rm = num_geno_rm + num_hwe_rm + num_maf_rm,
      num_geno_rm = num_geno_rm,
      num_hwe_rm = num_hwe_rm,
      num_maf_rm = num_maf_rm,
      log = log
    )
  )
  
}


#' Quality control on sample level with PLINK
#'
#' @param bfile              [\code{string}]\cr
#'                           The basename of the binary PLINK files.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param call.rate          [\code{double}]\cr
#'                           Filter out all samples with call rates falling below the provided threshold.
#' @param het.sigma          [\code{double}]\cr
#'                           Filter out all samples with mean heterozygosity rate above or below \code{het.sigma} times estimated heterozygosity standard deviation.
#' @param ld.pruning.params  [\code{list}]\cr
#'                           List with function arguments passed to \code{\link{plink_ld_pruning}}.
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to ALL PLINK calls (LD pruning, heterozygosity estimation, exclusion of samples).
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param tmp.dir            [\code{string}]\cr
#'                           Path where to save temporary files. If not set by user, defaults to \code{tempdir()}.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @details Heterozygosity estimation is not LD-sensitive, thus LD pruning is performend first using \code{\link{plink_ld_pruning}}. Heterozygosity estimate is based on SNPs only. See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system outputs as \code{list} of \code{character} vectors and number of samples excluded per criteria.
#' @export
#'
#' @import checkmate tools data.table stats
#' @importFrom batchtools makeRegistry reduceResultsList batchMap submitJobs waitForJobs
#'
plink_sample_qc <- function(bfile, output.prefix, 
                            call.rate, het.sigma, 
                            ld.pruning.params, ...,
                            bed.file = NULL, bim.file = NULL, fam.file = NULL,
                            exec = "plink2", tmp.dir = NULL,
                            num.threads,
                            memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(bfile)) {
    checkmate::assert_file(bed_file <- sprintf("%s.bed", bfile), add = assertions)
    checkmate::assert_file(bim_file <- sprintf("%s.bim", bfile), add = assertions)
    checkmate::assert_file(fam_file <- sprintf("%s.fam", bfile), add = assertions)
    input <- sprintf("--bfile %s", bfile)
  } else {
    checkmate::assert_file(bed_file <- bed.file, add = assertions)
    checkmate::assert_file(bim_file <- bim.file, add = assertions)
    checkmate::assert_file(fam_file <- fam.file, add = assertions)
    input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file)
  }
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  if (!missing(call.rate)) {
    checkmate::assert_number(call.rate, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
    mind <- sprintf("--mind %f", 1 - call.rate)
  } else {
    mind <- NULL
  }
  if (!missing(het.sigma)) {
    checkmate::assert_number(het.sigma, lower = 0, finite = TRUE, null.ok = FALSE, add = assertions)
  }
  
  checkmate::assert_list(ld.pruning.params, min.len = 3, names = "unique", any.missing = FALSE, all.missing = FALSE, null.ok = FALSE, add = assertions)
  checkmate::assert_subset(c("window.size", "step.size", "threshold"), names(ld.pruning.params), add = assertions)
  
  assert_command(exec, add = assertions)
  
  if (!missing(tmp.dir)) {
    checkmate::assert_string(tmp.dir, add = assertions)
    tmp_dir <- tmp.dir
    dir.create(tmp.dir, recursive = TRUE)
  } else {
    tmp_dir <- tempdir()
  }
  
  if (missing(num.threads)) {     
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)   
  }   
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {     
    memory = max(5000,                   
                 -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000),                       
                      -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE),                   
                 na.rm = TRUE)   
  }   
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  reg_dir <- tempfile(pattern = "reg", tmpdir = tmp_dir)
  on.exit(unlink(reg_dir, recursive = TRUE, force = TRUE), add = TRUE)
  file.create(conf_file <- tempfile(tmpdir = tmp_dir))
  on.exit(unlink(conf_file, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines(sprintf("cluster.functions = batchtools::makeClusterFunctionsSocket(ncpus = %d)", num.threads), con = conf_file)
  ld_reg <- batchtools::makeRegistry(
    file.dir = reg_dir,
    work.dir = dirname(output.prefix),
    packages = c("imbs"),
    conf.file = conf_file
  )
  
  chromosomes <- intersect(
    data.table::fread(cmd = sprintf("awk '{chr[$1]++} END {for (key in chr) print key}' %s", bim_file), col.names = "CHR")$CHR, # chromosomes in dataset
    1:22 # autosomes
  )
  
  if (length(chromosomes) == 0) {
    stop("No autosomes in dataset!")
  }
  
  batchtools::batchMap(
    fun = plink_ld_pruning, 
    chr = chromosomes,
    output.prefix = sprintf("%s_chr%d", output.prefix, chromosomes),
    more.args = c(
      ld.pruning.params, 
      list(bed.file = bed_file, 
           bim.file = bim_file, 
           fam.file = fam_file, 
           snps.only = "--snps-only just-acgt",
           num.threads = 1,
           memory = floor(memory/num.threads),
           exec = exec
      )
    ),
    reg = ld_reg
  )
  
  batchtools::submitJobs(reg = ld_reg)
  
  if (!batchtools::waitForJobs(reg = ld_reg)) {
    stop(sprintf("LD pruning failed! Check registry at %s", reg_dir))
  }
  
  ld_log <- batchtools::reduceResultsList(reg = ld_reg)
  
  file.remove(sprintf("%s.prune.in", output.prefix))
  file.create(sprintf("%s.prune.in", output.prefix))
  lapply(
    X = 1:22, 
    FUN =  function(chr) file.append(sprintf("%s.prune.in", output.prefix), sprintf("%s_chr%d.prune.in", output.prefix, chr))
  )
  
  het_log <- system_call(
    bin = exec,
    args = c(input,
             "--threads", num.threads,
             "--memory", memory,
             sprintf("--extract %s.prune.in", output.prefix),
             "--keep-allele-order",
             "--allow-extra-chr",
             "--het",
             "--out", output.prefix, ...)
  )
  
  het <- data.table::fread(sprintf("%s.het", output.prefix))
  het[, HET_RATE := (`N(NM)` - `O(HOM)`)/`N(NM)`]
  het[, MEAN_HET_RATE := mean(HET_RATE, na.rm = TRUE)]
  het[, SD_HET_RATE := sd(HET_RATE, na.rm = TRUE)]
  data.table::fwrite(
    x = het[HET_RATE < MEAN_HET_RATE - het.sigma*SD_HET_RATE | HET_RATE > MEAN_HET_RATE + het.sigma*SD_HET_RATE],
    file = sprintf("%s.het_remove", output.prefix),
    sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE
  )
  
  qc_log <- system_call(
    bin = exec,
    args = c(input,
             "--threads", num.threads,
             "--memory", memory,
             sprintf("--remove %s.het_remove", output.prefix),
             mind,
             "--keep-allele-order",
             "--allow-extra-chr",
             "--make-bed",
             "--out", output.prefix, ...)
  )
  
  num_mind_rm <- max(0, as.integer(gsub("^(\\d+) people removed due to missing genotype data.*", "\\1", grep("^\\d+ people removed due to missing genotype data.*", qc_log, value = TRUE))))
  num_het_rm <- het[HET_RATE < MEAN_HET_RATE - het.sigma*SD_HET_RATE | HET_RATE > MEAN_HET_RATE + het.sigma*SD_HET_RATE, .N]
  
  return(
    list(
      num_sample_rm = num_mind_rm + num_het_rm,
      num_mind_rm = num_mind_rm,
      num_het_rm = num_het_rm,
      ld_log = ld_log,
      het_log = het_log, 
      qc_log = qc_log
    )
  )
  
}

#' Internal PCA computation with PLINK
#'
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the output files.
#' @param num.evec           [\code{int}]\cr
#'                           Number of principal componentns to calculate.
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK call.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @return Captured system outputs as \code{list} of \code{character} vectors and number of samples excluded per criteria.
#' 
.plink_pca <- function(bed.file, bim.file, fam.file, output.prefix, num.evec, ..., exec, num.threads, memory) {
  
  system_call(
    bin = exec,
    args = c(
      sprintf("--bed %s", bed.file),
      sprintf("--bim %s", bim.file),
      sprintf("--fam %s", fam.file),
      sprintf("--out %s", output.prefix),
      sprintf("--pca %d header", num.evec),
      "--threads", num.threads,
      "--memory", memory,
      "--keep-allele-order",
      "--allow-extra-chr", ...
    )
  )
  
}

#' Internal Fst computation with PLINK
#'
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the output files.
#' @param pop.file           [\code{string}]\cr
#'                           File defining disjoint clusters/strata of samples. FIDs in the first column, IIDs in the second column, and cluster names in the third column.
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK call.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @return Captured system outputs as \code{list} of \code{character} vectors and number of samples excluded per criteria.
#' 
.plink_fst <- function(bed.file, bim.file, fam.file, output.prefix, pop.file, ..., exec, num.threads, memory) {
  
  system_call(
    bin = exec,
    args = c(
      sprintf("--bed %s", bed.file),
      sprintf("--bim %s", bim.file),
      sprintf("--fam %s", fam.file),
      sprintf("--out %s", output.prefix),
      sprintf("--fst"),
      sprintf("--within %s", pop.file),
      "--threads", num.threads,
      "--memory", memory,
      "--keep-allele-order",
      "--allow-extra-chr", ...
    )
  )
  
}

#' PCA with PLINK
#'
#' @param bfile              [\code{string}]\cr
#'                           The basename of the binary PLINK files.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param num.evec           [\code{int}]\cr
#'                           Number of principal componentns to calculate.
#' @param outlier.removal    [\code{flag}]\cr
#'                           Indicate if outliers shall be removed before final PCA computation.
#' @param outlier.sigma      [\code{int}]\cr
#'                           How many standard deviations from the mean of a single PC is not an outlier?
#' @param num.outlier.evec   [\code{int}]\cr
#'                           How many PCs shall be used for outlier detection?
#' @param num.outlier.iter   [\code{int}]\cr
#'                           Number of outlier iterations.
#' @param ld.pruning.params  [\code{list}]\cr
#'                           List with function arguments passed to \code{\link{plink_ld_pruning}}.
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param tmp.dir            [\code{string}]\cr
#'                           Path where to save temporary files. If not set by user, defaults to \code{tempdir()}.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @details First independent SNPs are extracted using PLINK LD pruning (\code{plink_ld_pruning}). If requested, outliers are detected and removed by calculating a PCA, calculating the mean (\eqn{\hat{\mu}_{j}}) and standard deviation (\eqn{\hat{\sigma}_{j}}) of each principal component, and excluding those samples who's PC value \eqn{v_{j}} exceeds \eqn{\hat{\mu}_{j}\pm\nu\,\hat{\sigma}_{j}} for at least one \eqn{j}, where \eqn{\nu} is \code{outlier.sigma}.
#'
#' @return A \code{list} of logs, Eigenvectors and Eigenvalues of both the outlier removal and the final PCA.
#' @export
#' 
#' @import checkmate tools data.table
#' @importFrom batchtools makeRegistry reduceResultsList batchMap submitJobs waitForJobs
#'
plink_pca <- function(bfile, output.prefix, 
                      num.evec, 
                      outlier.removal = FALSE, outlier.sigma, num.outlier.evec, num.outlier.iter,
                      ld.pruning.params,
                      bed.file = NULL, bim.file = NULL, fam.file = NULL,
                      exec = "plink2", tmp.dir = NULL,
                      num.threads,
                      memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(bfile)) {
    checkmate::assert_file(bed_file <- sprintf("%s.bed", bfile), add = assertions)
    checkmate::assert_file(bim_file <- sprintf("%s.bim", bfile), add = assertions)
    checkmate::assert_file(fam_file <- sprintf("%s.fam", bfile), add = assertions)
    input <- sprintf("--bfile %s", bfile)
  } else {
    checkmate::assert_file(bed_file <- bed.file, add = assertions)
    checkmate::assert_file(bim_file <- bim.file, add = assertions)
    checkmate::assert_file(fam_file <- fam.file, add = assertions)
    input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file)
  }
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  checkmate::assert_int(num.evec, na.ok = FALSE, lower = 1, upper = Inf, null.ok = FALSE, add = assertions)
  
  checkmate::assert_flag(outlier.removal, add = assertions)
  if (outlier.removal) {
    checkmate::assert_number(outlier.sigma, na.ok = FALSE, lower = 0, upper = Inf, finite = TRUE, null.ok = FALSE, add = assertions)
    checkmate::assert_int(num.outlier.evec, na.ok = FALSE, lower = 1, upper = Inf, null.ok = FALSE, add = assertions)
    checkmate::assert_int(num.outlier.iter, na.ok = FALSE, lower = 1, upper = Inf, null.ok = FALSE, add = assertions)
  }  else {
    num.outlier.iter = 0
  }
  
  checkmate::assert_list(ld.pruning.params, min.len = 3, names = "unique", any.missing = FALSE, all.missing = FALSE, null.ok = FALSE, add = assertions)
  checkmate::assert_subset(c("window.size", "step.size", "threshold"), names(ld.pruning.params), add = assertions)
  
  assert_command(exec, add = assertions)
  
  if (!missing(tmp.dir)) {
    checkmate::assert_string(tmp.dir, add = assertions)
    tmp_dir <- tmp.dir
    dir.create(tmp.dir, recursive = TRUE)
  } else {
    tmp_dir <- tempdir()
  }
  
  if (missing(num.threads)) {     
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)   
  }   
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {     
    memory = max(5000,                   
                 -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000),                       
                      -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE),                   
                 na.rm = TRUE)   
  }   
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  reg_dir <- tempfile(pattern = "reg", tmpdir = tmp_dir)
  on.exit(unlink(reg_dir, recursive = TRUE, force = TRUE), add = TRUE)
  file.create(conf_file <- tempfile(tmpdir = tmp_dir))
  on.exit(unlink(conf_file, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines(sprintf("cluster.functions = batchtools::makeClusterFunctionsSocket(ncpus = %d)", num.threads), con = conf_file)
  ld_reg <- batchtools::makeRegistry(
    file.dir = reg_dir,
    work.dir = dirname(output.prefix),
    packages = c("imbs"),
    conf.file = conf_file
  )
  
  chromosomes <- intersect(
    data.table::fread(cmd = sprintf("awk '{chr[$1]++} END {for (key in chr) print key}' %s", bim_file), col.names = "CHR")$CHR, # chromosomes in dataset
    1:22 # autosomes
  )
  
  if (length(chromosomes) == 0) {
    stop("No autosomes in dataset!")
  }
  
  batchtools::batchMap(
    fun = plink_ld_pruning, 
    chr = chromosomes,
    output.prefix = sprintf("%s_chr%d", output.prefix, chromosomes),
    more.args = c(
      ld.pruning.params, 
      list(bed.file = bed_file, 
           bim.file = bim_file, 
           fam.file = fam_file, 
           snps.only = "--snps-only just-acgt",
           num.threads = 1,
           memory = floor(memory/num.threads),
           exec = exec
      )
    ),
    reg = ld_reg
  )
  
  batchtools::submitJobs(reg = ld_reg)
  
  if (!batchtools::waitForJobs(reg = ld_reg)) {
    stop(sprintf("LD pruning failed! Check registry at %s", reg_dir))
  }
  
  ld_log <- batchtools::reduceResultsList(reg = ld_reg)
  
  prune_in_file <- sprintf("%s.prune.in", output.prefix)
  file.remove(prune_in_file)
  file.create(prune_in_file)
  lapply(
    X = 1:22, 
    FUN =  function(chr) file.append(prune_in_file, sprintf("%s_chr%d.prune.in", output.prefix, chr))
  )
  
  rm_samples_file <- sprintf("%s.rm_samples", output.prefix)
  file.create(rm_samples_file)
  
  outlier_evec_list <- vector("list", num.outlier.iter)
  outlier_eval_list <- vector("list", num.outlier.iter)
  outlier_samples_list <- vector("list", num.outlier.iter)
  outlier_log_list <- vector("list", num.outlier.iter)
  
  if (outlier.removal) {
    
    for (i in seq_len(num.outlier.iter)) {
      
      outlier_prefix <- sprintf("%s_outlier_removal_iteration%03d", output.prefix, i)
      rm_samples_iteration_file <- sprintf("%s.rm_samples", outlier_prefix)
      
      outlier_log_list[[i]] <- .plink_pca(
        bed.file = bed_file, 
        bim.file = bim_file, 
        fam.file = fam_file, 
        output.prefix = outlier_prefix, 
        num.evec = num.outlier.evec,
        rm.samples = sprintf("--remove %s", rm_samples_file), 
        pruned.snps = sprintf("--extract %s", prune_in_file),
        exec = exec,
        memory = memory,
        num.threads = num.threads
      )
      
      evecs <- data.table::fread(sprintf("%s.eigenvec", outlier_prefix))
      evecs <- data.table::melt(evecs, id.vars = c("FID", "IID"))
      
      evals <- cbind(
        "variable" = sprintf("PC%d", 1:num.outlier.evec), 
        "value" = data.table::fread(sprintf("%s.eigenval", outlier_prefix))
      )
      
      outlier_evec_list[[i]] <- data.table::copy(evecs)
      outlier_eval_list[[i]] <- data.table::copy(evals)
      
      evecs[, MEAN := mean(value, na.rm = TRUE), by = "variable"]
      evecs[, SD := sd(value, na.rm = TRUE), by = "variable"]
      evecs[, RM := value < MEAN - outlier.sigma*SD | value > MEAN + outlier.sigma*SD]
      
      outlier_samples <- unique(evecs[which(RM), .(FID, IID)])
      outlier_samples_list[[i]] <- outlier_samples
      
      if (length(outlier_samples) <= 1) {
        # no more outliers, exit loop early
        break
      }
      
      data.table::fwrite(
        outlier_samples, 
        file = rm_samples_file, 
        quote = FALSE, sep = "\t",
        row.names = FALSE, col.names = FALSE, append = TRUE
      )
      data.table::fwrite(
        outlier_samples, 
        file = rm_samples_iteration_file, 
        quote = FALSE, sep = "\t",
        row.names = FALSE, col.names = FALSE, append = FALSE
      )
    }
    
  }
  
  pca_log <- .plink_pca(
    bed.file = bed_file, bim.file = bim_file, fam.file = fam_file,
    output.prefix = output.prefix, num.evec = num.evec,
    rm.samples = sprintf("--remove %s", rm_samples_file), 
    pruned.snps = sprintf("--extract %s", prune_in_file),
    exec = exec, num.threads = num.threads, memory = memory
  )
  
  evecs <- data.table::fread(sprintf("%s.eigenvec", output.prefix))
  evals <- cbind(
    "variable" = sprintf("PC%d", 1:num.evec), 
    "value" = data.table::fread(sprintf("%s.eigenval", output.prefix))
  )
  evecs <- data.table::melt(evecs, id.vars = c("FID", "IID"))
  
  return(
    list(
      outlier_samples = outlier_samples_list,
      outlier_evecs = outlier_evec_list,
      outlier_evals = outlier_eval_list,
      outlier_logs = outlier_log_list,
      evecs = evecs,
      evals = evals,
      ld_log = ld_log,
      pca_log = pca_log
    )
  )
  
}

#' Fst with PLINK
#'
#' @param bfile              [\code{string}]\cr
#'                           The basename of the binary PLINK files.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param pop.file           [\code{string}]\cr
#'                           File defining disjoint clusters/strata of samples. FIDs in the first column, IIDs in the second column, and cluster names in the third column.
#' @param all.comb           [\code{flag}]\cr
#'                           Indicate if all population combinations should be considered for Fst calculation.   
#' @param with.na            [\code{flag}]\cr
#'                           Indicate if samples with missing population identifier ('NA') should be considered as a valid population.
#' @param outlier.removal    [\code{flag}]\cr
#'                           Indicate if outliers shall be removed before final PCA computation.
#' @param outlier.sigma      [\code{int}]\cr
#'                           How many standard deviations from the mean of a single PC is not an outlier?
#' @param num.outlier.evec   [\code{int}]\cr
#'                           How many PCs shall be used for outlier detection?
#' @param num.outlier.iter   [\code{int}]\cr
#'                           Number of outlier iterations.
#' @param ld.pruning.params  [\code{list}]\cr
#'                           List with function arguments passed to \code{\link{plink_ld_pruning}}.
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param tmp.dir            [\code{string}]\cr
#'                           Path where to save temporary files. If not set by user, defaults to \code{tempdir()}.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @details First independent SNPs are extracted using PLINK LD pruning (\code{plink_ld_pruning}). If requested, outliers are detected and removed by calculating a PCA, calculating the mean (\eqn{\hat{\mu}_{j}}) and standard deviation (\eqn{\hat{\sigma}_{j}}) of each principal component, and excluding those samples who's PC value \eqn{v_{j}} exceeds \eqn{\hat{\mu}_{j}\pm\nu\,\hat{\sigma}_{j}} for at least one \eqn{j}, where \eqn{\nu} is \code{outlier.sigma}.
#'
#' @return A \code{list} of logs, Eigenvectors and Eigenvalues of the outlier removal and the final Fst values.
#' @export
#' 
#' @import checkmate tools data.table
#' @importFrom batchtools makeRegistry reduceResultsList batchMap submitJobs waitForJobs
#'
plink_fst <- function(bfile, output.prefix, 
                      pop.file, all.comb = TRUE, with.na = FALSE,
                      outlier.removal = FALSE, outlier.sigma, num.outlier.evec, num.outlier.iter,
                      ld.pruning.params,
                      bed.file = NULL, bim.file = NULL, fam.file = NULL,
                      exec = "plink2", tmp.dir = NULL,
                      num.threads,
                      memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(bfile)) {
    checkmate::assert_file(bed_file <- sprintf("%s.bed", bfile), add = assertions)
    checkmate::assert_file(bim_file <- sprintf("%s.bim", bfile), add = assertions)
    checkmate::assert_file(fam_file <- sprintf("%s.fam", bfile), add = assertions)
    input <- sprintf("--bfile %s", bfile)
  } else {
    checkmate::assert_file(bed_file <- bed.file, add = assertions)
    checkmate::assert_file(bim_file <- bim.file, add = assertions)
    checkmate::assert_file(fam_file <- fam.file, add = assertions)
    input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file)
  }
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  checkmate::assert_string(pop.file, add = assertions)
  checkmate::assert_file(pop.file, add = assertions)
  
  checkmate::assert_flag(all.comb, add = assertions)
  
  checkmate::assert_flag(with.na, add = assertions)
  
  checkmate::assert_flag(outlier.removal, add = assertions)
  if (outlier.removal) {
    checkmate::assert_number(outlier.sigma, na.ok = FALSE, lower = 0, upper = Inf, finite = TRUE, null.ok = FALSE, add = assertions)
    checkmate::assert_int(num.outlier.evec, na.ok = FALSE, lower = 1, upper = Inf, null.ok = FALSE, add = assertions)
    checkmate::assert_int(num.outlier.iter, na.ok = FALSE, lower = 1, upper = Inf, null.ok = FALSE, add = assertions)
  }  else {
    num.outlier.iter = 0
  }
  
  checkmate::assert_list(ld.pruning.params, min.len = 3, names = "unique", any.missing = FALSE, all.missing = FALSE, null.ok = FALSE, add = assertions)
  checkmate::assert_subset(c("window.size", "step.size", "threshold"), names(ld.pruning.params), add = assertions)
  
  assert_command(exec, add = assertions)
  
  if (!missing(tmp.dir)) {
    checkmate::assert_string(tmp.dir, add = assertions)
    tmp_dir <- tmp.dir
    dir.create(tmp.dir, recursive = TRUE)
  } else {
    tmp_dir <- tempdir()
  }
  
  if (missing(num.threads)) {     
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)   
  }   
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {     
    memory = max(5000,                   
                 -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000),                       
                      -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE),                   
                 na.rm = TRUE)   
  }   
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  ld_reg_dir <- tempfile(pattern = "reg", tmpdir = tmp_dir)
  on.exit(unlink(ld_reg_dir, recursive = TRUE, force = TRUE), add = TRUE)
  file.create(conf_file <- tempfile(tmpdir = tmp_dir))
  on.exit(unlink(conf_file, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines(sprintf("cluster.functions = batchtools::makeClusterFunctionsSocket(ncpus = %d)", num.threads), con = conf_file)
  ld_reg <- batchtools::makeRegistry(
    file.dir = ld_reg_dir,
    work.dir = dirname(output.prefix),
    packages = c("imbs"),
    conf.file = conf_file
  )
  
  chromosomes <- intersect(
    data.table::fread(cmd = sprintf("awk '{chr[$1]++} END {for (key in chr) print key}' %s", bim_file), col.names = "CHR")$CHR, # chromosomes in dataset
    1:22 # autosomes
  )
  
  if (length(chromosomes) == 0) {
    stop("No autosomes in dataset!")
  }
  
  batchtools::batchMap(
    fun = plink_ld_pruning, 
    chr = chromosomes,
    output.prefix = sprintf("%s_chr%d", output.prefix, chromosomes),
    more.args = c(
      ld.pruning.params, 
      list(bed.file = bed_file, 
           bim.file = bim_file, 
           fam.file = fam_file, 
           snps.only = "--snps-only just-acgt",
           num.threads = 1,
           memory = floor(memory/num.threads),
           exec = exec
      )
    ),
    reg = ld_reg
  )
  
  batchtools::submitJobs(reg = ld_reg)
  
  if (!batchtools::waitForJobs(reg = ld_reg)) {
    stop(sprintf("LD pruning failed! Check registry at %s", ld_reg_dir))
  }
  
  ld_log <- batchtools::reduceResultsList(reg = ld_reg)
  
  prune_in_file <- sprintf("%s.prune.in", output.prefix)
  file.remove(prune_in_file)
  file.create(prune_in_file)
  lapply(
    X = 1:22, 
    FUN =  function(chr) file.append(prune_in_file, sprintf("%s_chr%d.prune.in", output.prefix, chr))
  )
  
  rm_samples_file <- sprintf("%s.rm_samples", output.prefix)
  file.create(rm_samples_file)
  
  outlier_evec_list <- vector("list", num.outlier.iter)
  outlier_eval_list <- vector("list", num.outlier.iter)
  outlier_samples_list <- vector("list", num.outlier.iter)
  outlier_log_list <- vector("list", num.outlier.iter)
  
  if (outlier.removal) {
    
    for (i in seq_len(num.outlier.iter)) {
      
      outlier_prefix <- sprintf("%s_outlier_removal_iteration%03d", output.prefix, i)
      rm_samples_iteration_file <- sprintf("%s.rm_samples", outlier_prefix)
      
      outlier_log_list[[i]] <- .plink_pca(
        bed.file = bed_file, 
        bim.file = bim_file, 
        fam.file = fam_file, 
        output.prefix = outlier_prefix, 
        num.evec = num.outlier.evec,
        rm.samples = sprintf("--remove %s", rm_samples_file), 
        pruned.snps = sprintf("--extract %s", prune_in_file),
        exec = exec,
        memory = memory,
        num.threads = num.threads
      )
      
      evecs <- data.table::fread(sprintf("%s.eigenvec", outlier_prefix))
      evecs <- data.table::melt(evecs, id.vars = c("FID", "IID"))
      
      evals <- cbind(
        "variable" = sprintf("PC%d", 1:num.outlier.evec), 
        "value" = data.table::fread(sprintf("%s.eigenval", outlier_prefix))
      )
      
      outlier_evec_list[[i]] <- data.table::copy(evecs)
      outlier_eval_list[[i]] <- data.table::copy(evals)
      
      evecs[, MEAN := mean(value, na.rm = TRUE), by = "variable"]
      evecs[, SD := sd(value, na.rm = TRUE), by = "variable"]
      evecs[, RM := value < MEAN - outlier.sigma*SD | value > MEAN + outlier.sigma*SD]
      
      outlier_samples <- unique(evecs[which(RM), .(FID, IID)])
      outlier_samples_list[[i]] <- outlier_samples
      
      if (length(outlier_samples) <= 1) {
        # no more outliers, exit loop early
        break
      }
      
      data.table::fwrite(
        outlier_samples, 
        file = rm_samples_file, 
        quote = FALSE, sep = "\t",
        row.names = FALSE, col.names = FALSE, append = TRUE
      )
      data.table::fwrite(
        outlier_samples, 
        file = rm_samples_iteration_file, 
        quote = FALSE, sep = "\t",
        row.names = FALSE, col.names = FALSE, append = FALSE
      )
    }
    
  }
  
  populations <- data.table::fread(pop.file, col.names = c("FID", "IID", "POP"))
  
  if (!with.na) {
    populations <- populations[!is.na(POP)]
  }
  
  num_populations <- populations[, length(unique(POP))]
  
  if (all.comb) {
    
    population_combinations <- unlist(
      lapply(
        X = 2:num_populations,
        FUN = function(m) combn(unique(populations[, POP]), m = m, simplify = FALSE)
      ),
      recursive = FALSE
    )
    
  } else {
    
    population_combinations <- list(populations[, (unique(POP))])
    
  }
  
  fst_reg_dir <- tempfile(pattern = "reg", tmpdir = tmp_dir)
  on.exit(unlink(fst_reg_dir, recursive = TRUE, force = TRUE), add = TRUE)
  fst_reg <- batchtools::makeRegistry(
    file.dir = fst_reg_dir,
    work.dir = dirname(output.prefix),
    packages = c("imbs"),
    conf.file = conf_file
  )
  
  batchtools::batchExport(
    export = list(
      tmp_dir = tmp_dir,
      populations = populations,
      with.na = with.na,
      .plink_fst = .plink_fst
    ),
    reg = fst_reg
  )
  
  batchtools::batchMap(
    fun = function(pops, ...) {
      pop_file <- tempfile(tmpdir = tmp_dir)
      on.exit(unlink(pop_file), add = TRUE)
      
      output_prefix <- tempfile(tmpdir = tmp_dir)
      on.exit(file.remove(list.files(tmp_dir, pattern = basename(output_prefix), full.names = TRUE)), add = TRUE) 
      
      data.table::fwrite(
        x = populations[POP %in% pops],
        pop_file, 
        quote = FALSE, 
        sep = "\t",
        row.names = FALSE, 
        col.names = FALSE
      )
      
      if (with.na) {
        
        pop_file <- sprintf("%s keep-NA", pop_file)
        
      }
      
      fst_log <- .plink_fst(
        pop.file = pop_file, output.prefix = output_prefix, ...
      )
      
      fst <- data.table::fread(sprintf("%s.fst", output_prefix))
      
      return(list(
        log = fst_log,
        fst = fst
      ))
    }, 
    pops = population_combinations,
    more.args = list(
      bed.file = bed_file, 
      bim.file = bim_file, 
      fam.file = fam_file,
      rm.samples = sprintf("--remove %s", rm_samples_file), 
      pruned.snps = sprintf("--extract %s", prune_in_file),
      exec = exec, 
      num.threads = 1, 
      memory = floor(memory/num.threads)
    ),
    reg = fst_reg
  )
  
  batchtools::submitJobs(reg = fst_reg)
  
  if (!batchtools::waitForJobs(reg = fst_reg)) {
    stop(sprintf("Fst calculation failed! Check registry at %s", fst_reg_dir))
  }
  
  fst <- batchtools::reduceResultsList(reg = fst_reg)
  
  return(
    list(
      outlier_samples = outlier_samples_list,
      outlier_evecs = outlier_evec_list,
      outlier_evals = outlier_eval_list,
      outlier_logs = outlier_log_list,
      fst = fst,
      population_combinations = population_combinations,
      ld_log = ld_log
    )
  )
  
}

#' GWA (single variant logistic regression) with PLINK
#'
#' @param bfile              [\code{string}]\cr
#'                           The basename of the binary PLINK files.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param test               [\code{string}]\cr
#'                           Genetic model to use. Choices are 'genotypic', 'hethom', 'dominant', and 'recessive'.
#'                           The 'genotypic' modifier adds an additive effect/dominance deviation 2df joint test (with two genotype-dependent variables in the regression, one with 0/1/2 coding and the second with 0/1/0 coding), while 'hethom' uses 0/0/1 and 0/1/0 coding instead. 
#'                           The 'dominant' and 'recessive' modifiers specify a model assuming full dominance or recessiveness, respectively, for the A1 allele.
#' @param sex                [\code{flag}]\cr
#'                           Include sex as covariate?
#' @param beta               [\code{flag}]\cr
#'                           Output regression coefficients instead of OR?
#' @param intercept          [\code{flag}]\cr
#'                           Include intercept in output?
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK call.
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param tmp.dir            [\code{string}]\cr
#'                           Path where to save temporary files. If not set by user, defaults to \code{tempdir()}.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @details Performs logistic regression given a case/control phenotype. For details see PLINK manual.
#'
#' @return A \code{list} of log (\code{gwa_log}, as \code{character}) and regression results (\code{gwa_result}, as \code{data.table}).
#' @export
#' 
#' @import checkmate tools data.table
#'
plink_logistic <- function(bfile, output.prefix,
                           test, sex, beta, intercept,
                           ...,
                           bed.file = NULL, bim.file = NULL, fam.file = NULL,
                           exec = "plink2",
                           num.threads,
                           memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(bfile)) {
    checkmate::assert_file(bed_file <- sprintf("%s.bed", bfile), add = assertions)
    checkmate::assert_file(bim_file <- sprintf("%s.bim", bfile), add = assertions)
    checkmate::assert_file(fam_file <- sprintf("%s.fam", bfile), add = assertions)
    input <- sprintf("--bfile %s", bfile)
  } else {
    checkmate::assert_file(bed_file <- bed.file, add = assertions)
    checkmate::assert_file(bim_file <- bim.file, add = assertions)
    checkmate::assert_file(fam_file <- fam.file, add = assertions)
    input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file)
  }
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  if (!missing(test)) {
    checkmate::assert_choice(test, choices = c("genotypic", "hethom", "dominant", "recessive"), null.ok = FALSE, add = assertions)
    test <- sprintf("%s", test)
  } else {
    test <- ""
  }
  if (!missing(sex)) {
    checkmate::assert_flag(sex, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (sex) {
      sex <- "sex"
    } else {
      sex <- "no-x-sex"
    }
  } else {
    sex <- ""
  }
  if (!missing(beta)) {
    checkmate::assert_flag(beta, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (beta) {
      beta <- "beta"
    } else {
      beta <- ""
    }
  } else {
    beta <- ""
  }
  if (!missing(intercept)) {
    checkmate::assert_flag(intercept, na.ok = FALSE, null.ok = FALSE, add = assertions)
    if (intercept) {
      intercept <- "intercept"
    } else {
      intercept <- ""
    }
  } else {
    intercept <- ""
  }
  
    assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {     
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)   
  }   
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {     
    memory = max(5000,                   
                 -min(-(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000),                       
                      -(num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000), na.rm = TRUE),                   
                 na.rm = TRUE)   
  }   
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  gwa_log <- system_call(
    bin = exec,
    args = c(input,
             "--threads", num.threads,
             "--memory", memory,
             "--logistic", test, sex, beta, intercept,
             "--out", output.prefix, ...)
  )
  
  gwa_result <- data.table::fread(file = sprintf("%s.assoc.logistic", output.prefix))
  
  return(
    list(
      gwa_log = gwa_log,
      gwa_result = gwa_result
    )
  )
  
}