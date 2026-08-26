process COUNT {
    tag { lane }

    publishDir path: "${params.outputDir}/counts/${saf.baseName}/${lane}",
               mode: 'copy',
               overwrite: true

    input:
    tuple val(lane), path(sams), path(saf)

    output:
    path("${lane}.txt"), emit: countedFiles
    path("${lane}.txt.summary"), emit: featureCountsResults

    script:
    """
    GOOD=()
    BAD=()
    for f in ${sams}; do
        n=\$(samtools view -c -F 0x904 "\$f")
        if [ "\$n" -gt 0 ]; then GOOD+=("\$f"); else BAD+=("\$f"); fi
    done

    if [ \${#GOOD[@]} -gt 0 ]; then
        featureCounts \
            -T ${task.cpus} \
            -a ${saf} \
            -p \
            --countReadPairs \
            -B \
            -C \
            -F SAF \
            -o ${lane}.txt \
            "\${GOOD[@]}"
    else
        echo "# no aligned reads in any sample of this lane" > ${lane}.txt
        awk 'BEGIN{OFS="\t"} {print \$1,\$2,\$3,\$4,\$5,(\$4-\$3+1)}' ${saf} >> ${lane}.txt
        echo -e "Status" > ${lane}.txt.summary
    fi

    for f in "\${BAD[@]}"; do
        awk -v name="\$f" 'NR==1{print; next} NR==2{print \$0"\t"name; next} {print \$0"\t0"}' ${lane}.txt > ${lane}.txt.tmp && mv ${lane}.txt.tmp ${lane}.txt
        awk -v name="\$f" 'NR==1{print \$0"\t"name; next} {print \$0"\t0"}' ${lane}.txt.summary > ${lane}.txt.summary.tmp && mv ${lane}.txt.summary.tmp ${lane}.txt.summary
    done
    """
}