process parquet2gr {
    memory "30 GB"
    cpus 1
    input:
	tuple val(meta),
	      val(parquets)
    output:
	path "all_concatemers.rds", 
	emit: "concatemers"
    """
    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/parquet2gr_task.R")
    Rscript \$RSCRIPT_PATH ${parquets}

    """
}   


process initial_diagnostics {
    memory "30 GB"
    cpus 1
    input:
	tuple val(concatemers)
    output:
	path "concatemer_analysis.rds"
    """
    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/analyze_concats.R")
    Rscript \$RSCRIPT_PATH ${concatemers}

    """
}
