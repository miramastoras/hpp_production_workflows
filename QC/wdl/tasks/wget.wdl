version 1.0

workflow runWget{
    call wget
}
task wget {
    input {
        Array[File] files
        String outDir
        # runtime configurations
        Int memSize=16
        Int threadCount=4
        Int diskSize=128
        String dockerImage="mobinasri/bio_base:v0.1"
        Int preemptible=2
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        mkdir ~{tarGzName}
        cp ~{sep=" " files} ~{tarGzName}
        tar -cf ~{tarGzName}.tar ~{tarGzName}
        gzip ~{tarGzName}.tar
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File fileTarGz = "~{tarGzName}.tar.gz"
    }
}
