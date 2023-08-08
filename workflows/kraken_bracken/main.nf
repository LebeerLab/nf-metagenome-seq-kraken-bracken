include { KRAKEN } from './modules/kraken'
include { BRACKEN; DETERMINE_MIN_LENGTH } from './modules/bracken' //addParams(OUTPUT: "${params.outdir}")

workflow KRACKEN_BRACKEN {
    take: reads

    main:
        ch_KrakenDB = Channel.value(file ("${params.krakendb}"))

        KRAKEN(filteredReads, ch_KrakenDB)
            .branch {
                success: it[1].countLines() > 1
                failed: it[1].countLines() == 1
            }
            .set { kraken_reports }

        kraken_reports.failed.subscribe { println "Sample ${it[0]} only contained unclassified reads and is excluded for further bracken processing."}
        // Run bracken on appropriate readlength sizes   

        // Determine % reads of sizes 200, 100, 50
        DETERMINE_MIN_LENGTH(reads.success)
            // Filter out empty reads (readlength=0)
            .branch { 
                success : it[2] as int > 0
                failed : it[2] as int == 0
            }
            .set { max_length }
        kraken_reports.success
            .join(max_length.success)
            .set {kraken_reportsLength}

        kraken_reportsLength
	        .map{ tuple(it[0], it[1], it[3] )}
    	    .set{ kreports }

        BRACKEN(kraken_reportsLength, ch_KrakenDB)

    emit: BRACKEN.out
}